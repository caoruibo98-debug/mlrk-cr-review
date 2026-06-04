from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


REPO = Path(__file__).resolve().parents[1]
RECOVERY = REPO / "data/curation/master_rule_recovery_candidates_round5.csv"
OUT = REPO / "data/curation/rule_recovery_decisions_round5.csv"


def _uniq(values: pd.Series) -> str:
    cleaned = [str(v) for v in values.dropna().unique() if str(v).strip() and str(v).lower() != "nan"]
    return ";".join(sorted(cleaned))


def _decision(row: pd.Series) -> tuple[str, str, str]:
    sources = set(str(row["rule_sources"]).split(";")) if row["rule_sources"] else set()
    has_retrorules = "retrorules" in sources
    has_microberx = "microberx" in sources
    has_ec = bool(str(row["rule_ec_numbers"]).strip())
    if has_retrorules and has_ec:
        return (
            "priority_rule_import_candidate",
            "RetroRules rule has EC/source-reaction provenance and regenerated the target product.",
            "Check RetroRules/MetaNetX source reaction, license, direction, and exact literature pair before training.",
        )
    if has_microberx and not has_ec:
        return (
            "rule_candidate_needs_source_recovery",
            "MicrobeRX rule regenerated the target product, but the master row lacks EC and source-reaction metadata.",
            "Recover MicrobeRX source reaction or replace with RetroRules/Rhea/BRENDA/MetaNetX rule before production import.",
        )
    return (
        "manual_review_required",
        "A master rule regenerated the target product but provenance is incomplete or mixed.",
        "Review rule source and exact biological evidence before importing.",
    )


def build_decisions(recovery_path: Path) -> pd.DataFrame:
    rec = pd.read_csv(recovery_path)
    hits = rec[rec["generated_target_product"].eq(True)].copy()
    grouped = (
        hits.groupby(["module", "sb", "pb", "substrate_name", "product_name"], dropna=False)
        .agg(
            hit_rule_count=("rule_id", "size"),
            rule_sources=("rule_source", _uniq),
            rule_ec_numbers=("rule_ec_number", _uniq),
            rule_ec_classes=("rule_ec_class", _uniq),
            rule_legacy_ids=("rule_legacy_id", _uniq),
            rule_ids=("rule_id", _uniq),
            max_rule_score=("rule_score", "max"),
            min_rule_diameter=("rule_diameter", "min"),
            gap_families=("gap_families", "first"),
            round4_row_ids=("round4_row_ids", "first"),
        )
        .reset_index()
    )
    decisions = grouped.apply(_decision, axis=1, result_type="expand")
    grouped["round5_decision"] = decisions[0]
    grouped["decision_reason"] = decisions[1]
    grouped["remaining_verification"] = decisions[2]
    grouped["positive_sample_status"] = "not_added_by_rule_alone"
    grouped["negative_sample_status"] = "not_applicable"
    grouped["next_pipeline_step"] = grouped["round5_decision"].map(
        {
            "priority_rule_import_candidate": "build_rule_import_manifest_with_source_reaction_check",
            "rule_candidate_needs_source_recovery": "recover_missing_rule_provenance_or_find_database_alternative",
            "manual_review_required": "manual_rule_and_literature_review",
        }
    )
    return grouped.sort_values(["round5_decision", "module", "substrate_name", "product_name"]).reset_index(drop=True)


def main() -> None:
    parser = argparse.ArgumentParser(description="Summarize recovered master rules into import decisions.")
    parser.add_argument("--recovery", type=Path, default=RECOVERY)
    parser.add_argument("--out", type=Path, default=OUT)
    args = parser.parse_args()
    out = build_decisions(args.recovery)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out, index=False)
    print(
        {
            "rows": len(out),
            "decisions": out["round5_decision"].value_counts().to_dict(),
        }
    )
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
