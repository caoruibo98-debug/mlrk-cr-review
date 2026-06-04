from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

import pandas as pd


REPO = Path(__file__).resolve().parents[1]
DEFAULT_POSITIVE_POOL = Path(
    r"D:\CRB\FoodGut\positive_sample_db\data\processed\foodgut_modular_positive_candidate_pool_v2.csv"
)
DEFAULT_ORIGINAL_KERNEL = Path(r"D:\CRB\Food models\scripts\ssrf\ml_ranking_kernel")


PUBMED_SPOT_CHECKED = {
    "30113130",
    "24909569",
    "22113864",
    "23542626",
    "11086885",
    "17196483",
    "21457417",
    "11722907",
    "37789701",
    "35118817",
}

PUBCHEM_SPOT_CHECKED_INCHIKEYS = {
    "ZQSIJRDFPHDXIC",  # daidzein
    "ADFCQWZHKCXPAJ",  # equol
    "CWVRJTMFETXNAD",  # chlorogenic acid
    "AFSDNFLWKVMVRB",  # ellagic acid
    "RIUPLDUFZCXCHM",  # urolithin A
    "REFJWTPEDVJJIY",  # quercetin
    "FTVWIRXFELQLPI",  # naringenin
    "HVDGDHBAMCBBLR",  # enterolactone
}


@dataclass(frozen=True)
class FamilyQuery:
    family: str
    terms: tuple[str, ...]
    priority_terms: tuple[str, ...]


FAMILY_QUERIES = [
    FamilyQuery(
        "phenolic_ester_hydrolysis",
        ("chlorogenic", "caffeic", "quinic", "hydroxycinnamate_ester", "ester hydrolysis"),
        ("chlorogenic", "caffeic"),
    ),
    FamilyQuery(
        "ellagitannin_urolithin_multistep",
        ("urolithin", "ellagic", "gordonibacter", "ellagibacter", "enterocloster"),
        ("urolithin", "ellagic"),
    ),
    FamilyQuery(
        "isoflavone_reductive_metabolism",
        ("equol", "daidzein", "dihydrodaidzein", "tetrahydrodaidzein", "slackia"),
        ("equol", "tetrahydrodaidzein"),
    ),
    FamilyQuery(
        "isoflavone_reduction",
        ("dihydrodaidzein", "dihydrogenistein", "isoflavone-specific reductase", "daidzein reductase"),
        ("dihydrodaidzein", "dihydrogenistein"),
    ),
    FamilyQuery(
        "isoflavone_c_glycoside_conversion",
        ("puerarin", "c-glycoside", "c-glucoside", "dorea"),
        ("puerarin", "c-glycoside"),
    ),
    FamilyQuery(
        "lignan_deglucosylation",
        ("secoisolariciresinol diglucoside", "sdg", "lignan", "deglycosylation"),
        ("secoisolariciresinol diglucoside", "sdg"),
    ),
    FamilyQuery(
        "lignan_multistep_demethylation_dehydroxylation",
        ("secoisolariciresinol", "enterodiol", "demethyl", "dehydroxyl", "lignan"),
        ("secoisolariciresinol", "enterodiol"),
    ),
    FamilyQuery(
        "lignan_oxidation",
        ("enterodiol", "enterolactone", "lactonifactor", "lactonization"),
        ("enterodiol", "enterolactone"),
    ),
    FamilyQuery(
        "flavanol_ring_fission",
        ("catechin", "epicatechin", "valerolactone", "flavan-3-ol", "flavan_3_ol"),
        ("catechin", "valerolactone"),
    ),
    FamilyQuery(
        "flavanone_ring_fission",
        ("naringenin", "naringin", "hydroxyphenylpropionic", "phenylpropionic"),
        ("naringenin", "phenylpropionic"),
    ),
    FamilyQuery(
        "flavonol_ring_fission",
        ("quercetin", "taxifolin", "alphitonin", "dihydroxyphenylacetic", "eubacterium ramulus"),
        ("quercetin", "dihydroxyphenylacetic"),
    ),
    FamilyQuery(
        "ellagitannin_hydrolysis",
        ("punicalagin", "punicalin", "ellagitannin", "ellagic acid"),
        ("punicalagin", "ellagic acid"),
    ),
]


NEGATIVE_STRATEGIES = [
    {
        "negative_type": "unlabeled_generated_nonmatch",
        "use_as_y0": False,
        "logic": "Generated candidates that are not in the positive source should stay unlabeled unless a reviewed source rules them out.",
    },
    {
        "negative_type": "class_contrast_negative",
        "use_as_y0": "class_only",
        "logic": "Use reactions from a neighboring class as negatives for reaction-family classification, not as universal molecule-level negatives.",
    },
    {
        "negative_type": "context_negative_nonproducer",
        "use_as_y0": "context_model_only",
        "logic": "A non-producer strain/community can be negative only when the model includes strain/community context.",
    },
    {
        "negative_type": "hard_negative_after_source_review",
        "use_as_y0": "reviewed_only",
        "logic": "A plausible generated product may become a hard negative only after source review shows the supported product is different under the same context.",
    },
]


def _block1(value: object) -> str:
    if pd.isna(value):
        return ""
    return str(value).split("-")[0].strip()


def _clean_text(value: object) -> str:
    if pd.isna(value):
        return ""
    text = str(value).strip()
    if text.endswith(".0") and text[:-2].isdigit():
        text = text[:-2]
    return text


def _source_id(row: pd.Series) -> str:
    pmid = _clean_text(row.get("pmid"))
    doi = _clean_text(row.get("doi_curated"))
    ids = []
    if pmid:
        ids.append(f"PMID:{pmid}")
    if doi:
        ids.append(f"DOI:{doi}")
    return ";".join(ids)


def _evidence_score(df: pd.DataFrame) -> pd.Series:
    level = pd.to_numeric(df.get("evidence_level", 0), errors="coerce").fillna(0)
    confidence = df.get("confidence_level", "").fillna("").astype(str).str.lower()
    bonus = confidence.str.contains("direct|strain|high").astype(int)
    return level * 10 + bonus


def _contains_terms(blob: pd.Series, terms: tuple[str, ...]) -> pd.Series:
    mask = pd.Series(False, index=blob.index)
    for term in terms:
        mask |= blob.str.contains(term.lower(), regex=False, na=False)
    return mask


def build_positive_manifest(args: argparse.Namespace) -> pd.DataFrame:
    pos = pd.read_csv(args.positive_pool, low_memory=False)
    reactions = pd.read_parquet(args.original_kernel / "outputs/modular/ltr/reactions.parquet")
    clean = pd.read_parquet(args.original_kernel / "outputs/modular/ltr/clean_candidates_full.parquet")
    seeds = pd.read_csv(REPO / "data/reaction_family_gap_seed_manifest.tsv", sep="\t")

    pos = pos.copy()
    pos["sb"] = pos["substrate_inchikey"].map(_block1)
    pos["pb"] = pos["product_inchikey"].map(_block1)
    pos["source_id"] = pos.apply(_source_id, axis=1)
    pos["evidence_score"] = _evidence_score(pos)

    reaction_pairs = set(zip(reactions["sb"], reactions["pb"]))
    clean_positive_pairs = set(zip(clean.loc[clean["y"].eq(1), "sb"], clean.loc[clean["y"].eq(1), "pb"]))
    seed_pairs = set(zip(seeds["substrate_inchikey"].map(_block1), seeds["product_inchikey"].map(_block1)))

    text_cols = [
        "substrate_name",
        "product_name",
        "reaction_category",
        "reaction_type",
        "enzyme_name",
        "microbe_or_strain",
        "source_dataset",
        "paper_title_curated",
        "pmid",
        "doi_curated",
    ]
    available_text_cols = [c for c in text_cols if c in pos.columns]
    blob = pos[available_text_cols].fillna("").astype(str).agg(" | ".join, axis=1).str.lower()

    rows: list[dict[str, object]] = []
    for fq in FAMILY_QUERIES:
        sub = pos[_contains_terms(blob, fq.terms)].copy()
        if not sub.empty:
            evidence_level = pd.to_numeric(sub.get("evidence_level", 0), errors="coerce").fillna(0)
            has_source_id = sub["source_id"].fillna("").astype(str).str.len().gt(0)
            sub = sub[evidence_level.ge(args.min_evidence_level) & has_source_id].copy()
        if sub.empty:
            rows.append(
                {
                    "row_id": f"round2_{fq.family}_no_local_match",
                    "gap_family": fq.family,
                    "candidate_status": "absent_from_local_positive_pool",
                    "substrate_name": "",
                    "product_name": "",
                    "source_dataset": "",
                    "source_id": "",
                    "validation_status": "external_curation_required",
                    "training_allowed_round2": False,
                    "recommended_next_action": (
                        "extract exact substrate/product pairs from cited external source before training; "
                        "local matches were absent or below the Round 2 evidence/source-id threshold"
                    ),
                }
            )
            continue

        priority_mask = _contains_terms(blob.loc[sub.index], fq.priority_terms)
        sub["priority_hit"] = priority_mask.astype(int)
        sub = (
            sub.drop_duplicates(["sb", "pb", "source_dataset", "source_id"])
            .sort_values(["priority_hit", "evidence_score"], ascending=[False, False])
            .head(args.max_per_family)
        )

        for _, row in sub.iterrows():
            pair = (row["sb"], row["pb"])
            is_holdout = pair in seed_pairs
            in_reactions = pair in reaction_pairs
            in_clean_positive = pair in clean_positive_pairs
            pmid = _clean_text(row.get("pmid"))
            pubmed_checked = bool(pmid and pmid in PUBMED_SPOT_CHECKED)
            pubchem_checked = row["sb"] in PUBCHEM_SPOT_CHECKED_INCHIKEYS or row["pb"] in PUBCHEM_SPOT_CHECKED_INCHIKEYS

            if is_holdout:
                status = "do_not_train_exact_holdout"
                next_action = "replace with independent non-holdout examples"
            elif in_clean_positive:
                status = "already_in_clean_candidate_pool"
                next_action = "use as coverage control; do not count as new gap fix"
            elif in_reactions:
                status = "pipeline_recovery_candidate"
                next_action = "trace why positive reaction is absent from clean_candidates_full before retraining"
            else:
                status = "positive_pool_to_reactions_candidate"
                next_action = "verify source and decide whether ltr_build filters should admit this reaction"

            rows.append(
                {
                    "row_id": f"round2_{fq.family}_{len(rows)+1:04d}",
                    "gap_family": fq.family,
                    "candidate_status": status,
                    "substrate_name": _clean_text(row.get("substrate_name")),
                    "substrate_smiles": _clean_text(row.get("substrate_smiles")),
                    "substrate_inchikey": _clean_text(row.get("substrate_inchikey")),
                    "product_name": _clean_text(row.get("product_name")),
                    "product_smiles": _clean_text(row.get("product_smiles")),
                    "product_inchikey": _clean_text(row.get("product_inchikey")),
                    "reaction_category": _clean_text(row.get("reaction_category")),
                    "reaction_type": _clean_text(row.get("reaction_type")),
                    "macro_module": _clean_text(row.get("macro_module")),
                    "source_dataset": _clean_text(row.get("source_dataset")),
                    "source_id": _clean_text(row.get("source_id")),
                    "evidence_type": _clean_text(row.get("evidence_type")),
                    "evidence_level": _clean_text(row.get("evidence_level")),
                    "confidence_level": _clean_text(row.get("confidence_level")),
                    "enzyme_name": _clean_text(row.get("enzyme_name")),
                    "microbe_or_strain": _clean_text(row.get("microbe_or_strain")),
                    "in_reactions_parquet": in_reactions,
                    "in_clean_candidates_full_positive": in_clean_positive,
                    "exact_holdout_pair": is_holdout,
                    "pubmed_spot_checked": pubmed_checked,
                    "pubchem_spot_checked_block1": pubchem_checked,
                    "validation_status": "spot_checked" if pubmed_checked or pubchem_checked else "local_source_present_needs_external_check",
                    "training_allowed_round2": False,
                    "recommended_next_action": next_action,
                }
            )

    return pd.DataFrame(rows)


def build_negative_strategy() -> pd.DataFrame:
    rows = []
    for fq in FAMILY_QUERIES:
        for strategy in NEGATIVE_STRATEGIES:
            rows.append(
                {
                    "gap_family": fq.family,
                    "negative_type": strategy["negative_type"],
                    "use_as_y0": strategy["use_as_y0"],
                    "logic": strategy["logic"],
                    "validation_required": "source review; do not infer true negatives from absence",
                }
            )
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description="Build Round 2 sample-gap manifests.")
    parser.add_argument("--positive-pool", type=Path, default=DEFAULT_POSITIVE_POOL)
    parser.add_argument("--original-kernel", type=Path, default=DEFAULT_ORIGINAL_KERNEL)
    parser.add_argument("--max-per-family", type=int, default=12)
    parser.add_argument("--min-evidence-level", type=int, default=3)
    parser.add_argument("--out-dir", type=Path, default=REPO / "data/curation")
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    positive = build_positive_manifest(args)
    negative = build_negative_strategy()

    positive_path = args.out_dir / "reaction_gap_positive_manifest_round2.csv"
    negative_path = args.out_dir / "reaction_gap_negative_strategy_round2.csv"
    positive.to_csv(positive_path, index=False)
    negative.to_csv(negative_path, index=False)

    summary = {
        "positive_rows": len(positive),
        "negative_strategy_rows": len(negative),
        "families": positive["gap_family"].nunique(),
        "candidate_status_counts": positive["candidate_status"].value_counts(dropna=False).to_dict(),
        "validation_status_counts": positive["validation_status"].value_counts(dropna=False).to_dict(),
    }
    print(summary)
    print(f"wrote {positive_path}")
    print(f"wrote {negative_path}")


if __name__ == "__main__":
    main()
