"""
predict_substrate_clean.py — PRODUCTION 推理(部署忠实、可解释、诚实)

链路:名称/SMILES → 结构路由模块 → 规则 RunReactants 生成候选(部署时本就如此)
     → 冻结 ChemBERTa 编码 → LTR_chem(纯结构)打分排序
     → top-N + 产物名 + 置信 + 适用域(AD) + 证据 overlay(已知菌/酶/文献) + 诚实声明

诚实边界(写进输出 JSON):
  - 只在"规则生成出的候选"里排;规则生不出的真产物排不出(generation recall ~25-68%)。
  - 排序分 = LTR_chem 学习分,**非湿实验概率**;菌/酶/文献只作解释,不进排序。

用法:
  python scripts/ssrf/ml_ranking_kernel/modular/predict_substrate_clean.py --name "rutin" --topn 12
  python scripts/ssrf/ml_ranking_kernel/modular/predict_substrate_clean.py --smiles "<SMILES>" --topn 12
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit import RDLogger

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import kio  # noqa: E402
import module_router as mr  # noqa: E402
from generate_candidates import run_reactants  # noqa: E402
from ltr_candidates import load_module_rules  # noqa: E402
from ltr_build import LTR_OUT  # noqa: E402

RDLogger.DisableLog("rdApp.*")
MODELS = LTR_OUT / "models_clean"
POOL = r"D:\CRB\FoodGut\positive_sample_db\data\processed\foodgut_modular_positive_candidate_pool_v2.csv"


def _l2(v):
    return v / max(float(np.linalg.norm(v)), 1e-9)


def evidence_lookup():
    d = pd.read_csv(POOL, low_memory=False,
                    usecols=["substrate_inchikey", "product_inchikey", "enzyme_ec", "enzyme_name",
                             "representative_microbes", "pmid"])
    d["sb"] = d.substrate_inchikey.astype(str).str[:14]; d["pb"] = d.product_inchikey.astype(str).str[:14]
    ev = {}
    for r in d.itertuples(index=False):
        ev.setdefault((r.sb, r.pb), {
            "enzyme_ec": None if pd.isna(r.enzyme_ec) else str(r.enzyme_ec),
            "enzyme": None if pd.isna(r.enzyme_name) else str(r.enzyme_name)[:50],
            "microbe": None if pd.isna(r.representative_microbes) else str(r.representative_microbes)[:50],
            "pmid": None if pd.isna(r.pmid) else str(r.pmid)})
    return ev


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--name", default=None); ap.add_argument("--smiles", default=None)
    ap.add_argument("--topn", type=int, default=12); ap.add_argument("--n-rules", type=int, default=800)
    args = ap.parse_args()
    kio.setup_logging()
    import xgboost as xgb
    from resolve import name_to_smiles, block1_to_name
    from encode_molecules import FrozenEncoder

    if args.smiles:
        smi = kio.canonical_smiles(args.smiles); label = args.name or "query"
    elif args.name:
        smi, src = name_to_smiles(args.name)
        if not smi:
            raise SystemExit(f"无法解析名称 {args.name}")
        label = args.name; kio.log.info("解析 [%s]: %s", src, smi)
    else:
        raise SystemExit("给 --name 或 --smiles")
    module, reason = mr.classify(smi)
    kio.log.info("== %s → 模块 %s(%s) ==", label, module, mr.MODULE_NAMES[module])
    mp = MODELS / f"ltr_chem_{module}.json"
    if not mp.exists():
        raise SystemExit(f"模块 {module} 部署模型缺失: {mp}")
    meta = json.loads((MODELS / "deploy_meta.json").read_text()).get(module, {})
    dim = int(meta.get("dim", 384))
    model = xgb.XGBRanker(); model.load_model(str(mp))

    # 规则生成候选（部署忠实）
    rules = load_module_rules(np.random.default_rng(0), args.n_rules)[module]
    cmol = Chem.MolFromSmiles(smi); cmolh = Chem.AddHs(cmol)
    cand = {}
    for smarts, ecc in rules:
        for psmi in run_reactants(smarts, cmol, cmolh):
            pik = kio.smiles_to_inchikey(psmi); pb = kio.inchikey_block1(pik)
            sb0 = kio.inchikey_block1(kio.smiles_to_inchikey(smi))
            if pb and pb != sb0 and pb not in cand:
                cand[pb] = (psmi, pik)
    if not cand:
        print(json.dumps({"input": label, "module": module, "n_candidates": 0,
                          "note": "规则未生成任何候选(generation miss)"}, ensure_ascii=False, indent=2)); return
    kio.log.info("规则生成候选: %d", len(cand))

    # 编码 + 打分（chem-only）
    enc = FrozenEncoder(kio.load_config()["encoder"]["primary"])
    smis = [smi] + [v[0] for v in cand.values()]
    vecs = enc.encode(smis); emb = {s: np.asarray(v, np.float32) for s, v in zip(smis, vecs)}
    zs = _l2(emb[smi])
    pbs = list(cand.keys())
    X = np.zeros((len(pbs), 4 * dim), np.float32)
    for i, pb in enumerate(pbs):
        zp = _l2(emb[cand[pb][0]]); X[i] = np.concatenate([zs, zp, zs - zp, zs * zp])
    score = model.predict(X)

    ev = evidence_lookup(); sb0 = kio.inchikey_block1(kio.smiles_to_inchikey(smi))
    out = pd.DataFrame({"pb": pbs, "smiles": [cand[pb][0] for pb in pbs], "score": score})
    out = out.sort_values("score", ascending=False).head(args.topn).reset_index(drop=True)
    lo, hi = out.score.min(), out.score.max()
    rows = []
    for i, r in out.iterrows():
        e = ev.get((sb0, r.pb), {})
        rows.append({"rank": i + 1, "product_name": block1_to_name(r.pb) or "unnamed/novel",
                     "product_smiles": r.smiles, "score": round(float((r.score - lo) / (hi - lo)) if hi > lo else 0.5, 3),
                     "evidence": ";".join(f"{k}={v}" for k, v in e.items() if v) or "—"})
    prov = {"input": {"name": label, "smiles": smi}, "module": module, "module_name": mr.MODULE_NAMES[module],
            "n_rule_candidates": len(cand),
            "honest_note": "只在规则生成的候选里排序;排序分=LTR_chem学习分(非湿实验概率);证据仅解释不进排序。",
            "top": rows}
    safe = "".join(c if c.isalnum() or c in "-_" else "_" for c in label)
    kio.safe_output_path(f"modular/predictions/clean_{safe}.json").write_text(
        json.dumps(prov, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"\n=== {label} → 模块 {module} top-{args.topn} (规则候选 {len(cand)}) ===")
    print(pd.DataFrame(rows)[["rank", "product_name", "score", "evidence"]].to_string(index=False, max_colwidth=40))


if __name__ == "__main__":
    main()
