"""
full_rule_diagnose.py — 用全量 reaction_rule_master(239K) 重测 generation 覆盖

修正上一版只用 MicrobeRX 10K 的缺陷。流程：
  对每条 gold+silver 反应底物 → 全量规则 SMARTS 预筛(LHS HasSubstructMatch)→ matched 规则 RunReactants
  → 生成产物 block1 集 → 真产物是否命中(hit/miss)。
分 shard 并行；按反应类型/EC/分子式变化汇总 miss。

用法(分 shard)：
  python ... full_rule_diagnose.py --shard 0 --nshards 5
  python ... full_rule_diagnose.py --concat-only      # 汇总+出表
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit import RDLogger

MOD = Path(__file__).resolve().parents[1] / "modular"
SRC = Path(__file__).resolve().parents[1] / "src"
sys.path.insert(0, str(SRC)); sys.path.insert(0, str(MOD))
import kio  # noqa: E402
from generate_candidates import run_reactants  # noqa: E402
from diagnose_generation_miss import load_clean_reactions  # noqa: E402

RDLogger.DisableLog("rdApp.*")
OUT = kio.OUTPUTS / "gen_gap"
RULE_MASTER = kio.PROJECT_ROOT / "Dataset" / "reaction_rule_master.csv"


def load_rules(max_rules=None):
    d = pd.read_csv(RULE_MASTER, usecols=["rule_smarts"], low_memory=False).dropna()
    d = d[d.rule_smarts.str.contains(">>")].drop_duplicates("rule_smarts")
    if max_rules:
        d = d.head(max_rules)
    return d.rule_smarts.tolist()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--nshards", type=int, default=1)
    ap.add_argument("--max-rules", type=int, default=None)
    ap.add_argument("--silver-cap", type=int, default=250)
    ap.add_argument("--concat-only", action="store_true")
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()
    kio.setup_logging()
    cache = kio.safe_output_path("gen_gap/fullrule_cache/_m").parent

    if args.concat_only:
        parts = [pd.read_parquet(p) for p in cache.glob("*.parquet")]
        rx = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()
        kio.write_table(rx, OUT / "fullrule_hit_miss.parquet", force=True)
        print("=== 全量规则 generation recall by module×tier ===")
        print(rx.groupby(["module", "tier"])["hit"].agg(n="size", hit="mean").round(3).to_string())
        for by in ["reaction_category", "ec_class", "formula_delta_product_minus_substrate"]:
            t = rx.groupby(["module", by]).agg(n=("hit", "size"), miss=("hit", lambda x: round(1 - x.mean(), 2))).reset_index()
            t = t[t.n >= 4].sort_values("miss", ascending=False)
            kio.write_table(t, OUT / f"fullrule_miss_by_{by}.csv", force=True)
        print("\n=== miss 最重 by reaction_category(全量规则后仍缺的)===")
        tc = rx.groupby("reaction_category").agg(n=("hit", "size"), miss=("hit", lambda x: round(1 - x.mean(), 2))).reset_index()
        print(tc[tc.n >= 5].sort_values("miss", ascending=False).head(25).to_string(index=False))
        return

    rx = load_clean_reactions()
    # silver 子采样控算力(gold 全留)
    parts = []
    for m, g in rx.groupby("module"):
        gold = g[g.tier == "gold"]; sil = g[g.tier == "silver"]
        if len(sil) > args.silver_cap:
            sil = sil.sample(n=args.silver_cap, random_state=0)
        parts.append(pd.concat([gold, sil]))
    rx = pd.concat(parts, ignore_index=True)
    subs = rx.drop_duplicates(["module", "sb"])[["module", "sb", "substrate_smiles"]].reset_index(drop=True)
    subs = subs.iloc[args.shard::args.nshards].reset_index(drop=True)
    kio.log.info("shard %d/%d: %d 底物", args.shard, args.nshards, len(subs))

    rules = load_rules(args.max_rules)
    kio.log.info("全量规则: %d 条，预解析 LHS...", len(rules))
    lhs_cache = {}
    full_of_lhs = {}
    for s in rules:
        lhs = s.split(">>")[0].strip().lstrip("(").rstrip(")")
        if lhs not in lhs_cache:
            try:
                lhs_cache[lhs] = Chem.MolFromSmarts(lhs)
            except Exception:
                lhs_cache[lhs] = None
            full_of_lhs.setdefault(lhs, s)
    pats = [(lhs, p, full_of_lhs[lhs]) for lhs, p in lhs_cache.items() if p is not None]
    kio.log.info("有效 LHS 模式: %d", len(pats))

    sub_gen = {}
    for i, r in enumerate(subs.itertuples(index=False)):
        t0 = time.time()
        cmol = Chem.MolFromSmiles(r.substrate_smiles)
        if cmol is None:
            sub_gen[(r.module, r.sb)] = set(); continue
        cmolh = Chem.AddHs(cmol); gen = set()
        for lhs, patt, full in pats:
            try:
                if cmolh.HasSubstructMatch(patt) or cmol.HasSubstructMatch(patt):
                    for psmi in run_reactants(full, cmol, cmolh):
                        pb = kio.inchikey_block1(kio.smiles_to_inchikey(psmi))
                        if pb:
                            gen.add(pb)
            except Exception:
                continue
        sub_gen[(r.module, r.sb)] = gen
        if (i + 1) % 20 == 0:
            kio.log.info("  shard%d %d/%d (last %.1fs, gen=%d)", args.shard, i + 1, len(subs), time.time() - t0, len(gen))

    sub_rx = rx[rx.set_index(["module", "sb"]).index.isin(sub_gen.keys())].copy()
    sub_rx["hit"] = [int(row.pb in sub_gen.get((row.module, row.sb), set())) for row in sub_rx.itertuples(index=False)]
    sub_rx.to_parquet(cache / f"shard_{args.shard}.parquet", index=False)
    kio.log.info("shard %d done: %d 反应, hit率 %.3f", args.shard, len(sub_rx), sub_rx.hit.mean())


if __name__ == "__main__":
    main()
