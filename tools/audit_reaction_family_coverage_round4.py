from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import RDLogger


REPO = Path(__file__).resolve().parents[1]
RAW_POOL = Path(r"D:\CRB\FoodGut\positive_sample_db\data\processed\foodgut_modular_positive_candidate_pool_v2.csv")
LTR_OUT = REPO / "outputs/modular/ltr"
OUT = REPO / "data/curation/reaction_family_coverage_round4.csv"

sys.path.insert(0, str(REPO / "src"))
import kio  # noqa: E402


RDLogger.DisableLog("rdApp.*")

MACRO2MOD = {
    "small_molecule_polyphenol": "A",
    "carbohydrate_glycan": "B",
    "protein_amino_acid": "C",
    "lipid_fat": "D",
}


def _block1_from_smiles(smiles: object) -> str:
    smi = kio.canonical_smiles(smiles)
    if not smi:
        return ""
    return kio.inchikey_block1(kio.smiles_to_inchikey(smi))


def _load_ltr_ready_pool(path: Path) -> pd.DataFrame:
    raw = pd.read_csv(path, low_memory=False)
    single_step = raw["is_single_step"].astype(str).str.lower().isin(["yes", "derived"])
    raw = raw[single_step].copy()
    raw = raw[raw["structure_feasibility_flag"].astype(str).str.startswith("structure_parseable")].copy()
    raw["sb"] = raw["substrate_smiles"].map(_block1_from_smiles)
    raw["pb"] = raw["product_smiles"].map(_block1_from_smiles)
    raw = raw[raw["sb"].ne("") & raw["pb"].ne("") & raw["sb"].ne(raw["pb"])].copy()
    raw["module"] = raw["macro_module"].map(MACRO2MOD)
    raw = raw.dropna(subset=["module"]).copy()
    raw["evidence_level_num"] = pd.to_numeric(raw.get("evidence_level"), errors="coerce").fillna(0)
    raw["family_key"] = (
        raw["macro_module"].fillna("").astype(str)
        + " :: "
        + raw["reaction_category"].fillna("").astype(str)
        + " :: "
        + raw["reaction_type"].fillna("").astype(str)
    )
    return raw


def _pair_set(df: pd.DataFrame, y: int | None = None) -> set[tuple[str, str, str]]:
    if y is not None:
        df = df[df["y"].eq(y)].copy()
    return set(zip(df["module"].astype(str), df["sb"].astype(str), df["pb"].astype(str)))


def _suggest(row: pd.Series) -> str:
    if row["reactions_pair_count"] == 0 and row["high_evidence_pair_count"] > 0:
        return "fix_positive_pool_to_reactions_filters"
    if row["reactions_pair_count"] > 0 and row["clean_full_positive_count"] == 0:
        return "add_or_import_family_rules_before_retraining"
    if row["reactions_pair_count"] > 0 and row["reactions_to_clean_full_retention"] < 0.5:
        return "trace_rule_generation_sampling_and_drop_filters"
    if row["high_evidence_pair_count"] < 3:
        return "curate_more_literature_or_database_supported_pairs"
    return "monitor_as_currently_covered_family"


def build_coverage(raw_pool: Path, repo: Path) -> pd.DataFrame:
    pool = _load_ltr_ready_pool(raw_pool)
    reactions = pd.read_parquet(repo / "outputs/modular/ltr/reactions.parquet")
    clean_full = pd.read_parquet(repo / "outputs/modular/ltr/clean_candidates_full.parquet")
    clean = pd.read_parquet(repo / "outputs/modular/ltr/clean_candidates.parquet")
    fullrule_path = repo / "outputs/gen_gap/fullrule_hit_miss.parquet"
    fullrule = pd.read_parquet(fullrule_path) if fullrule_path.exists() else pd.DataFrame()

    reactions_pairs = set(zip(reactions["module"], reactions["sb"], reactions["pb"]))
    clean_full_pos = _pair_set(clean_full, y=1)
    clean_pos = _pair_set(clean, y=1)
    if fullrule.empty:
        fullrule_seen: set[tuple[str, str, str]] = set()
        fullrule_hit: set[tuple[str, str, str]] = set()
    else:
        fullrule_seen = set(zip(fullrule["module"], fullrule["sb"], fullrule["pb"]))
        fullrule_hit = set(zip(fullrule.loc[fullrule["hit"].eq(1), "module"], fullrule.loc[fullrule["hit"].eq(1), "sb"], fullrule.loc[fullrule["hit"].eq(1), "pb"]))

    pair_meta = pool.sort_values("evidence_level_num", ascending=False).drop_duplicates(
        ["module", "sb", "pb", "family_key"]
    )
    pair_tuples = list(zip(pair_meta["module"], pair_meta["sb"], pair_meta["pb"]))
    pair_meta["in_reactions"] = [x in reactions_pairs for x in pair_tuples]
    pair_meta["in_clean_full_positive"] = [x in clean_full_pos for x in pair_tuples]
    pair_meta["in_clean_positive"] = [x in clean_pos for x in pair_tuples]
    pair_meta["fullrule_seen"] = [x in fullrule_seen for x in pair_tuples]
    pair_meta["fullrule_hit"] = [x in fullrule_hit for x in pair_tuples]
    pair_meta["manual_literature"] = pair_meta["source_origin_type"].astype(str).eq("manual_literature_curated")
    pair_meta["high_evidence"] = pair_meta["manual_literature"] | pair_meta["evidence_level_num"].ge(3)

    grouped = (
        pair_meta.groupby(["module", "macro_module", "reaction_category", "reaction_type", "family_key"], dropna=False)
        .agg(
            raw_pair_count=("pb", "size"),
            reactions_pair_count=("in_reactions", "sum"),
            clean_full_positive_count=("in_clean_full_positive", "sum"),
            clean_positive_count=("in_clean_positive", "sum"),
            fullrule_seen_count=("fullrule_seen", "sum"),
            fullrule_hit_count=("fullrule_hit", "sum"),
            high_evidence_pair_count=("high_evidence", "sum"),
            manual_literature_pair_count=("manual_literature", "sum"),
        )
        .reset_index()
    )
    for num, den, out_col in [
        ("reactions_pair_count", "raw_pair_count", "raw_to_reactions_retention"),
        ("clean_full_positive_count", "reactions_pair_count", "reactions_to_clean_full_retention"),
        ("clean_positive_count", "reactions_pair_count", "reactions_to_clean_retention"),
        ("fullrule_hit_count", "fullrule_seen_count", "fullrule_generation_recall_seen"),
    ]:
        grouped[out_col] = np.where(grouped[den].astype(float) > 0, grouped[num] / grouped[den], np.nan)

    grouped["missing_from_clean_full"] = grouped["reactions_pair_count"] - grouped["clean_full_positive_count"]
    grouped["round4_recommended_action"] = grouped.apply(_suggest, axis=1)
    grouped = grouped.sort_values(
        ["missing_from_clean_full", "high_evidence_pair_count", "raw_pair_count"],
        ascending=False,
    )
    return grouped


def main() -> None:
    parser = argparse.ArgumentParser(description="Build reaction-family coverage audit for Round 4.")
    parser.add_argument("--raw-pool", type=Path, default=RAW_POOL)
    parser.add_argument("--repo", type=Path, default=REPO)
    parser.add_argument("--out", type=Path, default=OUT)
    args = parser.parse_args()

    coverage = build_coverage(args.raw_pool, args.repo)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    coverage.to_csv(args.out, index=False)
    summary = {
        "rows": len(coverage),
        "families_with_reactions": int(coverage["reactions_pair_count"].gt(0).sum()),
        "families_missing_clean_full": int(coverage["missing_from_clean_full"].gt(0).sum()),
        "high_evidence_families": int(coverage["high_evidence_pair_count"].gt(0).sum()),
        "top_actions": coverage["round4_recommended_action"].value_counts().to_dict(),
    }
    print(summary)
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
