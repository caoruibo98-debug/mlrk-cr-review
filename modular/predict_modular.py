"""
predict_modular.py — 部署推理：输入物质 → 模块路由 → LTR 检索排序 + 证据/解释 overlay

链路：
  名称/SMILES → canonical → module_router 路由 A/B/C/D
   → 载入该模块 LTR 终模型 + 产物词表
   → ChemBERTa 编码底物（即时）→ 对词表每个产物打分 → top-N
   → 每个候选附：产物名 + ml 分 + 置信 tier + 适用域(AD) + 证据 overlay(已知菌/酶/文献)

诚实边界：这是**检索**（在该模块已知产物词表里排），不是生成全新结构；ml 分非湿实验概率；
证据只有"已知反应"有（answer-key），所以只作解释/优先级，不参与排序。

用法：
  python scripts/ssrf/ml_ranking_kernel/modular/predict_modular.py --name "2'-FL" --topn 12
  python scripts/ssrf/ml_ranking_kernel/modular/predict_modular.py --smiles "<SMILES>" --topn 12
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import kio  # noqa: E402
import module_router as mr  # noqa: E402
from ltr_build import LTR_OUT  # noqa: E402

MODELS = LTR_OUT / "models"
POOL = r"D:\CRB\FoodGut\positive_sample_db\data\processed\foodgut_modular_positive_candidate_pool_v2.csv"


def _l2(v):
    return v / max(float(np.linalg.norm(v)), 1e-9)


def encode_one(smiles, emb_cache):
    if smiles in emb_cache:
        return emb_cache[smiles]
    from encode_molecules import FrozenEncoder
    cfg = kio.load_config()
    fe = FrozenEncoder(cfg["encoder"]["primary"], device=cfg["train"].get("device", "cuda"))
    v = fe.encode([smiles])[0]
    emb_cache[smiles] = np.asarray(v, np.float32)
    return emb_cache[smiles]


def load_evidence_lookup():
    """(sb,pb) → 已知证据（菌/酶/文献），给检索结果做 overlay。"""
    d = pd.read_csv(POOL, low_memory=False,
                    usecols=["substrate_inchikey", "product_inchikey", "enzyme_ec", "enzyme_name",
                             "microbe_or_strain", "representative_microbes", "pmid"])
    d["sb"] = d["substrate_inchikey"].astype(str).str[:14]
    d["pb"] = d["product_inchikey"].astype(str).str[:14]
    ev = {}
    for _, r in d.iterrows():
        ev.setdefault((r.sb, r.pb), {
            "enzyme_ec": None if pd.isna(r.enzyme_ec) else str(r.enzyme_ec),
            "enzyme": None if pd.isna(r.enzyme_name) else str(r.enzyme_name)[:60],
            "microbe": None if pd.isna(r.representative_microbes) else str(r.representative_microbes)[:60],
            "pmid": None if pd.isna(r.pmid) else str(r.pmid)})
    return ev


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--name", default=None)
    ap.add_argument("--smiles", default=None)
    ap.add_argument("--topn", type=int, default=12)
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()
    kio.setup_logging()
    import xgboost as xgb
    from resolve import name_to_smiles, block1_to_name

    if args.smiles:
        smi = kio.canonical_smiles(args.smiles); label = args.name or "query"
    elif args.name:
        smi, src = name_to_smiles(args.name)
        if smi is None:
            raise SystemExit(f"无法解析名称 '{args.name}'")
        label = args.name; kio.log.info("名称解析 [%s]: %s", src, smi)
    else:
        raise SystemExit("给 --name 或 --smiles")
    if smi is None:
        raise SystemExit("SMILES 解析失败")

    module, reason = mr.classify(smi)
    kio.log.info("== %s → 模块 %s (%s) | %s", label, module, mr.MODULE_NAMES[module], reason)
    mp = MODELS / f"ltr_{module}.json"
    if not mp.exists():
        raise SystemExit(f"模块 {module} 模型未训练: {mp}（先跑 ltr_production.py）")
    meta = json.loads((MODELS / f"ltr_{module}.meta.json").read_text())
    dim = int(meta["dim"])
    model = xgb.XGBRanker(); model.load_model(str(mp))

    vocab = pd.read_parquet(LTR_OUT / "vocab.parquet")
    V = vocab[vocab["module"] == module].reset_index(drop=True)

    # 编码底物 + 词表产物（复用缓存）
    c = pd.read_parquet(kio.FEATURES_DIR / "mol_embeddings.parquet")
    emb = {s: np.asarray(v, np.float32) for s, v in zip(c["smiles"], c["embedding"])}
    zs = _l2(encode_one(smi, emb))
    # 特征矩阵（chem-only：zs,zp,diff,prod）
    X = np.zeros((len(V), 4 * dim), np.float32); ok = np.zeros(len(V), bool)
    for i, ps in enumerate(V["product_smiles"]):
        zp = emb.get(ps)
        if zp is None:
            continue
        zp = _l2(zp); X[i] = np.concatenate([zs, zp, zs - zp, zs * zp]); ok[i] = True
    score = np.full(len(V), -1e9); score[ok] = model.predict(X[ok])
    V = V.assign(ml_score=score).sort_values("ml_score", ascending=False).head(args.topn).reset_index(drop=True)

    # 适用域：底物到该模块训练底物的最大余弦
    rxn = pd.read_parquet(LTR_OUT / "reactions.parquet")
    train_subs = rxn[rxn["module"] == module]["substrate_smiles"].dropna().unique()
    sims = [float(np.dot(zs, _l2(emb[s]))) for s in train_subs if s in emb]
    ad = max(sims) if sims else float("nan")
    ad_flag = "in-domain" if ad >= 0.8 else "borderline" if ad >= 0.6 else "OUT-OF-DOMAIN(低置信)"

    # 证据 overlay
    ev = load_evidence_lookup()
    sb = kio.inchikey_block1(kio.smiles_to_inchikey(smi))

    def conf(s):
        # min-max 在本次候选内做相对 tier
        return s
    smin, smax = V["ml_score"].min(), V["ml_score"].max()
    rows = []
    for i, r in V.iterrows():
        pb = r["product_block1"]
        e = ev.get((sb, pb), {})
        nm = block1_to_name(pb) or "unnamed/novel"
        norm = (r.ml_score - smin) / (smax - smin) if smax > smin else 0.5
        rows.append({"rank": i + 1, "product_name": nm, "product_smiles": r["product_smiles"],
                     "ml_score": round(float(norm), 3),
                     "evidence": "known: " + ";".join(f"{k}={v}" for k, v in e.items() if v) if e else "—"})
    out = pd.DataFrame(rows)

    prov = {"input": {"name": label, "smiles": smi, "inchikey": kio.smiles_to_inchikey(smi)},
            "module": module, "module_name": mr.MODULE_NAMES[module], "route_reason": reason,
            "applicability_domain": {"max_cosine": None if np.isnan(ad) else round(ad, 3), "flag": ad_flag},
            "task": "retrieval over module known-product vocab (not novel-structure generation)",
            "top": out.to_dict(orient="records")}
    safe = "".join(ch if ch.isalnum() or ch in "-_" else "_" for ch in label)
    kio.safe_output_path(f"modular/predictions/predict_{safe}.json").write_text(
        json.dumps(prov, ensure_ascii=False, indent=2), encoding="utf-8")

    print(f"\n=== {label} → 模块 {module}({mr.MODULE_NAMES[module]})  AD={ad_flag}({ad:.2f}) ===")
    print(out[["rank", "product_name", "ml_score", "evidence"]].to_string(index=False, max_colwidth=42))
    print(f"\nJSON: outputs/modular/predictions/predict_{safe}.json")


if __name__ == "__main__":
    main()
