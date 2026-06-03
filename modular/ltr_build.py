"""
ltr_build.py — 多证据 learning-to-rank 数据装配（用 Codex foodgut positive pool）

产物检索式 LTR：给定底物 S，在该模块的"产物词表"里把真产物排前面。
负样本 = 同模块里别的底物的真实代谢物（真实但错 → 硬负，抗 decoy bias）。

特征四组：结构(ChemBERTa)、反应(差分)、菌/酶证据、文献证据。
证据只对真实反应非空；负样本(错配 S-P)证据为空 → 评估时做 chem-only vs +evidence 消融。

切分：scaffold-CV(干净集) + 时间留出(train<2023, test≥2023 = 新文献召回测试)。

输出 modular/outputs/ltr/{pairs.parquet, vocab.parquet, reactions.parquet} + 复用 emb 缓存。
用法：python scripts/ssrf/ml_ranking_kernel/modular/ltr_build.py [--neg-k 15] [--force]
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
import kio  # noqa: E402

POOL = r"D:\CRB\FoodGut\positive_sample_db\data\processed\foodgut_modular_positive_candidate_pool_v2.csv"
LTR_OUT = kio.OUTPUTS / "modular" / "ltr"
LIT_SOURCES = {"benchmark_ssrf_rclss_2026-05-22",
               "curated_literature_gap7_import_ready_v1",
               "curated_literature_non7_extension_candidates_v1"}
MACRO2MOD = {"small_molecule_polyphenol": "A", "carbohydrate_glycan": "B",
             "protein_amino_acid": "C", "lipid_fat": "D"}


def load_pool() -> pd.DataFrame:
    p = kio.resolve_input(POOL)
    d = pd.read_csv(p, low_memory=False)
    # 单步 + 结构 ok
    ss = d["is_single_step"].astype(str).str.lower()
    d = d[ss.isin(["yes", "derived"])].copy()
    d = d[d["structure_feasibility_flag"].astype(str).str.startswith("structure_parseable")].copy()
    # canonical（保证与 ChemBERTa key 一致）
    d["substrate_smiles"] = d["substrate_smiles"].map(kio.canonical_smiles)
    d["product_smiles"] = d["product_smiles"].map(kio.canonical_smiles)
    d = d.dropna(subset=["substrate_smiles", "product_smiles"]).copy()
    d["sb"] = d["substrate_smiles"].map(lambda s: kio.inchikey_block1(kio.smiles_to_inchikey(s)))
    d["pb"] = d["product_smiles"].map(lambda s: kio.inchikey_block1(kio.smiles_to_inchikey(s)))
    d = d[d.sb.ne("") & d.pb.ne("") & d.sb.ne(d.pb)].copy()
    d["module"] = d["macro_module"].map(MACRO2MOD)
    d = d.dropna(subset=["module"]).copy()
    d["year"] = pd.to_numeric(d["pubmed_year"].fillna(d["year_curated"]), errors="coerce")
    # 干净 = 文献人工策展（非规则/DB派生）→ 金标测试集
    d["is_clean"] = d["source_origin_type"].astype(str).eq("manual_literature_curated") | \
        d["source_dataset"].isin(LIT_SOURCES)
    # tier：gold(文献金标) > silver(DB策展/待审,有泄漏) > weak(规则派生,仅预训练)
    tu = d["training_use_recommendation"].astype(str)
    silver = tu.isin(["database_positive_training_with_leakage_guard",
                      "review_step_scope_before_strict_training"])
    d["tier_rank"] = np.where(d["is_clean"], 2, np.where(silver, 1, 0)).astype(int)
    # 证据特征（数值化）
    d["ev_has_ec"] = d["enzyme_ec"].notna().astype(int)
    d["ev_ec_class"] = pd.to_numeric(d["enzyme_ec"].astype(str).str.extract(r"^(\d)")[0], errors="coerce").fillna(0).astype(int)
    d["ev_microbe_n"] = np.log1p(pd.to_numeric(d["microbe_record_count"], errors="coerce").fillna(0))
    d["ev_gene_n"] = np.log1p(pd.to_numeric(d["gene_record_count"], errors="coerce").fillna(0))
    d["ev_has_pmid"] = d["pmid"].notna().astype(int)
    d["ev_has_rule"] = d["reaction_rule_id"].notna().astype(int)
    d["ev_level"] = pd.to_numeric(d["evidence_level"], errors="coerce").fillna(2)
    return d


EV_COLS = ["ev_has_ec", "ev_ec_class", "ev_microbe_n", "ev_gene_n", "ev_has_pmid", "ev_has_rule", "ev_level"]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--neg-k", type=int, default=15, help="每正样本采样多少负产物")
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()
    kio.setup_logging()
    rng = np.random.default_rng(42)

    d = load_pool()
    kio.log.info("usable rows=%d  unique reactions=%d  modules=%s",
                 len(d), d.drop_duplicates(["sb", "pb"]).shape[0], d["module"].value_counts().to_dict())

    # 反应级 is_clean = 任一重复行是文献策展（否则会被高证据的 microberx 重复行盖掉）。
    # 这同时把"clean 反应的 weak 重复行"去掉 → 训练里不会见到测试反应 = 防泄漏。
    d["is_clean"] = d.groupby(["sb", "pb"])["is_clean"].transform("max").astype(bool)
    d["tier_rank"] = d.groupby(["sb", "pb"])["tier_rank"].transform("max").astype(int)
    # 唯一反应（每 (sb,pb) 取证据最全的一行；is_clean/tier 已组内统一）
    d["_score"] = d[EV_COLS].sum(axis=1)
    rxn = d.sort_values("_score", ascending=False).drop_duplicates(["sb", "pb"]).copy()
    rxn["tier"] = rxn["tier_rank"].map({2: "gold", 1: "silver", 0: "weak"})
    rxn["substrate_scaffold"] = rxn["substrate_smiles"].map(kio.murcko_scaffold)
    kio.log.info("unique reactions: %d (clean=%d)", len(rxn), int(rxn["is_clean"].sum()))

    # 每模块产物词表（block1 → smiles）
    vocab_rows = []
    for m, g in rxn.groupby("module"):
        for b1, sm in g.drop_duplicates("pb")[["pb", "product_smiles"]].values:
            vocab_rows.append({"module": m, "product_block1": b1, "product_smiles": sm})
    vocab = pd.DataFrame(vocab_rows)
    kio.log.info("product vocab per module: %s",
                 vocab.groupby("module").size().to_dict())

    # 真产物集合（每底物，用于负采样排除 + 评估）
    true_prod = rxn.groupby(["module", "sb"])["pb"].agg(set).to_dict()

    # 训练对：正 + 负（同模块词表采样，排除该底物真产物）
    pos = rxn[["module", "sb", "substrate_smiles", "substrate_scaffold", "pb", "product_smiles",
               "is_clean", "year"] + EV_COLS].copy()
    pos["y"] = 1
    vocab_by_mod = {m: g[["product_block1", "product_smiles"]].values for m, g in vocab.groupby("module")}
    neg_rows = []
    for _, r in pos.iterrows():
        m, s = r["module"], r["sb"]
        forbid = true_prod.get((m, s), set())
        cand = vocab_by_mod[m]
        idx = rng.choice(len(cand), size=min(args.neg_k * 3, len(cand)), replace=False)
        taken = 0
        for j in idx:
            pb, psm = cand[j]
            if pb in forbid:
                continue
            neg_rows.append({"module": m, "sb": s, "substrate_smiles": r["substrate_smiles"],
                             "substrate_scaffold": r["substrate_scaffold"], "pb": pb, "product_smiles": psm,
                             "is_clean": r["is_clean"], "year": r["year"], "y": 0,
                             **{c: 0 for c in EV_COLS}})  # 负样本=错配反应→无证据
            taken += 1
            if taken >= args.neg_k:
                break
    neg = pd.DataFrame(neg_rows)
    pairs = pd.concat([pos, neg], ignore_index=True)

    # scaffold fold（仅干净底物分 fold；权重/温度留训练用）
    from sklearn.model_selection import GroupKFold
    clean_subs = rxn[rxn["is_clean"]].drop_duplicates("sb")[["sb", "substrate_scaffold"]].copy()
    clean_subs["scaf"] = clean_subs["substrate_scaffold"].fillna("NS_" + clean_subs["sb"])
    k = min(5, clean_subs["scaf"].nunique())
    clean_subs["scaffold_fold"] = -1
    if k >= 2:
        gkf = GroupKFold(n_splits=k)
        for f, (_, te) in enumerate(gkf.split(clean_subs, groups=clean_subs["scaf"])):
            clean_subs.iloc[te, clean_subs.columns.get_loc("scaffold_fold")] = f
    pairs = pairs.merge(clean_subs[["sb", "scaffold_fold"]], on="sb", how="left")
    pairs["scaffold_fold"] = pairs["scaffold_fold"].fillna(-1).astype(int)  # 非干净=训练池

    kio.write_table(pairs, LTR_OUT / "pairs.parquet", force=args.force)
    kio.write_table(vocab, LTR_OUT / "vocab.parquet", force=args.force)
    kio.write_table(rxn[["module", "sb", "substrate_smiles", "substrate_scaffold", "pb", "product_smiles",
                         "is_clean", "tier", "year"] + EV_COLS], LTR_OUT / "reactions.parquet", force=args.force)
    kio.log.info("tier counts (unique reactions): %s", rxn["tier"].value_counts().to_dict())
    kio.log.info("tier × module: \n%s", rxn.groupby(["module", "tier"]).size().to_string())
    kio.log.info("pairs=%d (pos=%d neg=%d) | clean reactions per module: %s",
                 len(pairs), int((pairs.y == 1).sum()), int((pairs.y == 0).sum()),
                 rxn[rxn.is_clean].groupby("module").size().to_dict())


if __name__ == "__main__":
    main()
