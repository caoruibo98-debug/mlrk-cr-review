"""
diagnose_generation_miss.py — #1 Step1：规则生成天花板诊断（缺什么规则）

对每条 gold+silver 反应(我们真想预测的文献/策展反应),用**全部模块规则** RunReactants,
判断真产物能否被生成(hit/miss)。然后把 miss 按反应类型/EC class/分子式变化分桶
→ 精确告诉你"哪类反应规则生不出 = 该补哪类规则 / 哪类得靠生成模型"。

不碰现有模型;只读规则 + Codex 反应,输出诊断表。
用法：python scripts/ssrf/ml_ranking_kernel/gen_gap/diagnose_generation_miss.py [--n-rules 6000] [--force]
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit import RDLogger

# 复用 modular 的工具
MOD = Path(__file__).resolve().parents[1] / "modular"
SRC = Path(__file__).resolve().parents[1] / "src"
sys.path.insert(0, str(SRC)); sys.path.insert(0, str(MOD))
import kio  # noqa: E402
from generate_candidates import run_reactants  # noqa: E402
from ltr_candidates import load_module_rules  # noqa: E402

RDLogger.DisableLog("rdApp.*")
OUT = kio.OUTPUTS / "gen_gap"
POOL = r"D:\CRB\FoodGut\positive_sample_db\data\processed\foodgut_modular_positive_candidate_pool_v2.csv"
MACRO2MOD = {"small_molecule_polyphenol": "A", "carbohydrate_glycan": "B",
             "protein_amino_acid": "C", "lipid_fat": "D"}
SILVER = {"database_positive_training_with_leakage_guard", "review_step_scope_before_strict_training"}


def load_clean_reactions() -> pd.DataFrame:
    d = pd.read_csv(POOL, low_memory=False)
    d["module"] = d["macro_module"].map(MACRO2MOD)
    d = d.dropna(subset=["module", "substrate_smiles", "product_smiles"])
    ss = d["is_single_step"].astype(str).str.lower()
    d = d[ss.isin(["yes", "derived"])]
    # tier：gold(文献) / silver(DB策展)
    tu = d["training_use_recommendation"].astype(str)
    is_gold = d["source_origin_type"].astype(str).eq("manual_literature_curated")
    d["tier"] = np.where(is_gold, "gold", np.where(tu.isin(SILVER), "silver", "weak"))
    d = d[d["tier"].isin(["gold", "silver"])].copy()
    d["sb"] = d["substrate_smiles"].map(lambda s: kio.inchikey_block1(kio.smiles_to_inchikey(kio.canonical_smiles(s) or "")))
    d["pb"] = d["product_smiles"].map(lambda s: kio.inchikey_block1(kio.smiles_to_inchikey(kio.canonical_smiles(s) or "")))
    d = d[d.sb.ne("") & d.pb.ne("") & d.sb.ne(d.pb)]
    d["ec_class"] = pd.to_numeric(d["enzyme_ec"].astype(str).str.extract(r"^(\d)")[0], errors="coerce").fillna(0).astype(int)
    keep = ["module", "tier", "sb", "pb", "substrate_smiles", "reaction_category", "reaction_type",
            "ec_class", "formula_delta_product_minus_substrate"]
    return d[keep].drop_duplicates(["sb", "pb"]).reset_index(drop=True)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-rules", type=int, default=6000, help="每模块用多少规则(默认基本=全部)")
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()
    kio.setup_logging()
    rx = load_clean_reactions()
    kio.log.info("gold+silver 反应: %d (按模块 %s)", len(rx), rx.module.value_counts().to_dict())
    rules = load_module_rules(np.random.default_rng(0), args.n_rules)
    kio.log.info("规则数/模块: %s", {m: len(v) for m, v in rules.items()})

    # 每底物 RunReactants 一次(缓存生成产物 block1 集)
    sub_gen = {}   # (module, sb) -> set(pb)
    subs = rx.drop_duplicates(["module", "sb"])[["module", "sb", "substrate_smiles"]]
    for i, r in enumerate(subs.itertuples(index=False)):
        cmol = Chem.MolFromSmiles(r.substrate_smiles)
        if cmol is None:
            sub_gen[(r.module, r.sb)] = set(); continue
        cmolh = Chem.AddHs(cmol); gen = set()
        for smarts, _ in rules[r.module]:
            for psmi in run_reactants(smarts, cmol, cmolh):
                pb = kio.inchikey_block1(kio.smiles_to_inchikey(psmi))
                if pb:
                    gen.add(pb)
        sub_gen[(r.module, r.sb)] = gen
        if (i + 1) % 200 == 0:
            kio.log.info("  RunReactants %d/%d 底物", i + 1, len(subs))

    rx["hit"] = [int(row.pb in sub_gen.get((row.module, row.sb), set())) for row in rx.itertuples(index=False)]
    kio.write_table(rx, OUT / "reaction_hit_miss.parquet", force=args.force)

    # 汇总
    print("\n=== 每模块 generation recall(hit 率)===")
    g = rx.groupby("module")["hit"].agg(["mean", "size"]).round(3)
    print(g.to_string())

    def miss_table(by):
        t = rx.groupby(["module", by]).agg(n=("hit", "size"), miss=("hit", lambda x: 1 - x.mean())).reset_index()
        t = t[t["n"] >= 3].sort_values(["module", "miss"], ascending=[True, False])
        return t

    for by, name in [("reaction_category", "反应类别"), ("ec_class", "EC class"),
                     ("formula_delta_product_minus_substrate", "分子式变化")]:
        t = miss_table(by)
        kio.write_table(t, OUT / f"miss_by_{by}.csv", force=args.force)
        print(f"\n=== miss 率 by {name}(n≥3,miss 率高=规则最缺这类)===")
        print(t.head(30).to_string(index=False))


if __name__ == "__main__":
    main()
