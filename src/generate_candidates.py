"""
generate_candidates.py — Phase 1a

对每个**正样本底物**跑 production 规则引擎（SMARTS 匹配 + RCLSS 打分 + RunReactants），
得到该底物的候选 (substrate→product) 边 + rclss_production_score。

产出双重用途：
  1) in-substrate hard decoys（规则生成、非真实产物的边 = 负样本）
  2) 规则先验 RCLSS 基线排序（评估对照）+ 生成召回率（规则引擎是否生成了真产物）

逐底物 checkpoint（outputs/datasets/candidates_cache/<block1>.parquet），断点可续。
最终汇总 → outputs/datasets/substrate_candidates.parquet。

⚠️ 运行较慢（~70s/底物，含 RCLSS 全打分）。建议后台运行：
  python scripts/ssrf/ml_ranking_kernel/src/generate_candidates.py [--limit N] [--force]
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import RDLogger

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1].parent))  # scripts/ssrf for engine imports
import kio

RDLogger.DisableLog("rdApp.*")

_RXN_CACHE: dict[str, object] = {}


def get_reaction(smarts: str):
    if smarts in _RXN_CACHE:
        return _RXN_CACHE[smarts]
    rxn = None
    try:
        rxn = AllChem.ReactionFromSmarts(smarts)
        if rxn is not None:
            rxn.Initialize()
    except Exception:
        rxn = None
    _RXN_CACHE[smarts] = rxn
    return rxn


def standardize_product(mol) -> str | None:
    try:
        mol = Chem.RemoveHs(mol)
        Chem.SanitizeMol(mol)
        smi = Chem.MolToSmiles(mol)
        return smi if smi else None
    except Exception:
        return None


def run_reactants(smarts: str, c_mol, c_mol_h) -> list[str]:
    """对 substrate 应用规则，返回去重的 canonical 产物 SMILES 列表。"""
    rxn = get_reaction(smarts)
    if rxn is None:
        return []
    out: set[str] = set()
    for mol in (c_mol, c_mol_h):
        if mol is None:
            continue
        try:
            product_sets = rxn.RunReactants((mol,))
        except Exception:
            continue
        for pset in product_sets:
            for p in pset:
                smi = standardize_product(p)
                if smi:
                    out.add(smi)
        if out:
            break  # 优先无 H 版本；若空再试 AddHs
    return list(out)


def generate_for_substrate(eng, filter_compound, smiles: str, name: str) -> pd.DataFrame:
    df_rules = filter_compound(eng, smiles, name, top_k=10**9, min_score=0.0,
                               verbose=False, collapse_by_rule_id=False)
    if df_rules is None or df_rules.empty:
        return pd.DataFrame()

    c_mol = Chem.MolFromSmiles(smiles)
    c_mol_h = Chem.AddHs(c_mol) if c_mol is not None else None
    sub_ik = kio.smiles_to_inchikey(smiles)
    sub_b1 = kio.inchikey_block1(sub_ik)

    rows = []
    for _, r in df_rules.iterrows():
        rule_uid = r["rule_uid"]
        score = float(r["final_score"])
        try:
            smarts = eng.master.loc[rule_uid, "rule_smarts"]
            if isinstance(smarts, pd.Series):
                smarts = smarts.iloc[0]
        except Exception:
            continue
        if not isinstance(smarts, str) or ">>" not in smarts:
            continue
        for psmi in run_reactants(smarts, c_mol, c_mol_h):
            pik = kio.smiles_to_inchikey(psmi)
            pb1 = kio.inchikey_block1(pik)
            if not pb1 or pb1 == sub_b1:
                continue
            rows.append({
                "substrate_smiles": smiles,
                "substrate_inchikey": sub_ik,
                "substrate_block1": sub_b1,
                "product_smiles": psmi,
                "product_inchikey": pik,
                "product_block1": pb1,
                "rclss_production_score": score,
                "rule_uid": rule_uid,
            })
    if not rows:
        return pd.DataFrame()
    out = pd.DataFrame(rows)
    # 每 product_block1 保留 rclss 最高的来源规则
    out = out.sort_values("rclss_production_score", ascending=False)
    out = out.drop_duplicates("product_block1", keep="first").reset_index(drop=True)
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=None, help="只处理前 N 个底物（验证用）")
    ap.add_argument("--force", action="store_true", help="重算已 checkpoint 的底物")
    ap.add_argument("--edges", default=None, help="kernel_edges.parquet 路径（取正样本底物）")
    ap.add_argument("--shard", type=int, default=0, help="本 shard 编号 [0, nshards)")
    ap.add_argument("--nshards", type=int, default=1, help="总 shard 数（多进程并行）")
    ap.add_argument("--concat-only", action="store_true", help="只把 checkpoint 汇总成 substrate_candidates.parquet")
    args = ap.parse_args()

    kio.setup_logging()
    cache_dir = kio.safe_output_path("datasets/candidates_cache/_marker").parent

    if args.concat_only:
        parts = [pd.read_parquet(p) for p in cache_dir.glob("*.parquet") if not p.name.startswith("_")]
        allc = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()
        kio.write_table(allc, kio.DATASETS_DIR / "substrate_candidates.parquet", force=True)
        kio.log.info("== CONCAT: %d substrates, %d edges ==",
                     allc["substrate_block1"].nunique() if len(allc) else 0, len(allc))
        return

    edges_path = Path(args.edges) if args.edges else (kio.DATASETS_DIR / "kernel_edges.parquet")
    edges = pd.read_parquet(edges_path)
    pos = edges[edges["y"] == 1].drop_duplicates("substrate_block1")
    subs = pos[["substrate_smiles", "substrate_block1"]].dropna().reset_index(drop=True)
    if args.limit:
        subs = subs.head(args.limit)
    if args.nshards > 1:
        subs = subs.iloc[args.shard::args.nshards].reset_index(drop=True)
        kio.log.info("shard %d/%d -> %d substrates", args.shard, args.nshards, len(subs))
    kio.log.info("Generating candidates for %d unique positive substrates", len(subs))

    # 引擎只在需要时加载
    sys.path.insert(0, str(kio.PROJECT_ROOT / "scripts" / "ssrf"))
    from rclss_similarity import RCLSSEngine
    from ssrf_filter import filter_compound
    t0 = time.time()
    eng = RCLSSEngine.load_default()
    kio.log.info("engine loaded in %.1fs", time.time() - t0)

    done = 0
    for i, row in subs.iterrows():
        b1 = row["substrate_block1"]
        ckpt = cache_dir / f"{b1}.parquet"
        if ckpt.exists() and not args.force:
            kio.log.info("[%d/%d] %s cached, skip", i + 1, len(subs), b1)
            done += 1
            continue
        t1 = time.time()
        df = generate_for_substrate(eng, filter_compound, row["substrate_smiles"], b1)
        df.to_parquet(ckpt, index=False)
        kio.log.info("[%d/%d] %s -> %d candidate edges in %.1fs",
                     i + 1, len(subs), b1, len(df), time.time() - t1)
        done += 1

    kio.log.info("shard %d done (%d substrates processed)", args.shard, done)
    # 汇总（仅单进程；多 shard 用 --concat-only 单独跑，避免写竞争）
    if args.nshards == 1:
        parts = [pd.read_parquet(p) for p in cache_dir.glob("*.parquet") if not p.name.startswith("_")]
        if parts:
            allc = pd.concat(parts, ignore_index=True)
            kio.write_table(allc, kio.DATASETS_DIR / "substrate_candidates.parquet", force=True)
            kio.log.info("== DONE: %d substrates, %d total candidate edges ==",
                         allc["substrate_block1"].nunique(), len(allc))


if __name__ == "__main__":
    main()
