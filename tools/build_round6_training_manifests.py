from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


REPO = Path(__file__).resolve().parents[1]
CURATION = REPO / "data/curation"


def _block1(value: object) -> str:
    if pd.isna(value):
        return ""
    return str(value).split("-")[0].strip()


def _uniq(values: pd.Series) -> str:
    out = []
    for value in values.dropna().astype(str):
        value = value.strip()
        if value and value.lower() != "nan" and value not in out:
            out.append(value)
    return ";".join(sorted(out))


def _best_rule_rows(recovery: pd.DataFrame) -> pd.DataFrame:
    hits = recovery[recovery["generated_target_product"].eq(True)].copy()
    hits["has_ec"] = hits["rule_ec_number"].notna() & hits["rule_ec_number"].astype(str).ne("nan")
    hits["source_rank"] = hits["rule_source"].map({"retrorules": 0, "microberx": 1}).fillna(9)
    hits["ec_rank"] = hits["has_ec"].map({True: 0, False: 1})
    hits["diameter_rank"] = pd.to_numeric(hits["rule_diameter"], errors="coerce").fillna(999)
    hits["score_rank"] = -pd.to_numeric(hits["rule_score"], errors="coerce").fillna(0)
    hits = hits.sort_values(["module", "sb", "pb", "source_rank", "ec_rank", "diameter_rank", "score_rank"])
    cols = [
        "module",
        "sb",
        "pb",
        "rule_id",
        "rule_smarts",
        "rule_source",
        "rule_legacy_id",
        "rule_ec_number",
        "rule_ec_class",
        "rule_diameter",
        "rule_score",
        "rule_applicable_structure_type",
    ]
    return hits.drop_duplicates(["module", "sb", "pb"])[cols]


def _verification_summary(lit: pd.DataFrame) -> pd.DataFrame:
    lit = lit.copy()
    lit["sb"] = lit["pair_key"].astype(str).str.split("__").str[0]
    lit["pb"] = lit["pair_key"].astype(str).str.split("__").str[1]
    return (
        lit.groupby(["sb", "pb"], dropna=False)
        .agg(
            evidence_source_ids=("source_id", _uniq),
            evidence_titles=("source_title", _uniq),
            evidence_types=("evidence_type", _uniq),
            exact_reaction_support=("exact_reaction_supported", _uniq),
            positive_sample_use=("positive_sample_use", _uniq),
            rule_import_use=("rule_import_use", _uniq),
            verification_statuses=("verification_status", _uniq),
            caveats=("caveat", _uniq),
            next_actions=("next_action", _uniq),
        )
        .reset_index()
    )


def _classify_rule_import(row: pd.Series) -> tuple[str, bool, str, str]:
    decision = str(row["round5_decision"])
    verification = str(row.get("verification_statuses", ""))
    exact_support = str(row.get("exact_reaction_support", ""))
    if decision == "priority_rule_import_candidate":
        return (
            "candidate_after_source_reaction_check",
            False,
            "RetroRules recovered the pair, but source reaction, direction, license, and exact biological pair must be checked first.",
            "verify_retrorules_source_reaction_then_import_to_evidence_tiered_rule_pool",
        )
    if decision == "rule_candidate_needs_source_recovery" and "verified_pubmed_summary" in verification and "yes" in exact_support:
        return (
            "source_recovery_required_before_import",
            False,
            "Pair has promising PubMed support, but recovered rule is MicrobeRX-only or lacks EC/source-reaction provenance.",
            "recover_microberx_source_or_replace_with_rhea_brenda_metanetx_retrorules_rule",
        )
    return (
        "blocked_until_external_evidence_or_rule_provenance",
        False,
        "Rule or literature provenance remains too weak for production import.",
        "continue_literature_database_verification",
    )


def build_rule_import_manifest(decisions: pd.DataFrame, recovery: pd.DataFrame, lit: pd.DataFrame) -> pd.DataFrame:
    best = _best_rule_rows(recovery)
    ver = _verification_summary(lit)
    out = decisions.merge(best, on=["module", "sb", "pb"], how="left").merge(ver, on=["sb", "pb"], how="left")
    status = out.apply(_classify_rule_import, axis=1, result_type="expand")
    out["round6_rule_import_status"] = status[0]
    out["import_allowed_round6"] = status[1]
    out["import_blocker"] = status[2]
    out["next_rule_action"] = status[3]
    out["training_allowed_round6"] = False
    out["label_if_generated"] = "candidate_only_not_positive_by_rule"
    out["holdout_guard"] = "keep_round2_holdouts_out_of_training"
    return out


def _positive_status(row: pd.Series) -> tuple[str, bool, str]:
    support = str(row["exact_reaction_supported"]).lower()
    use = str(row["positive_sample_use"]).lower()
    decision = str(row["rule_recovery_decision"])
    if "yes" in support and "allowed" in use and decision == "priority_rule_import_candidate":
        return (
            "eligible_after_exact_extraction_and_rule_source_check",
            False,
            "Exact evidence appears promising and recovered rule has RetroRules provenance, but source extraction and rule source check are still required.",
        )
    if "yes" in support and "allowed" in use:
        return (
            "eligible_after_exact_extraction_and_rule_provenance_recovery",
            False,
            "Exact evidence appears promising, but the recovered rule needs source/EC provenance or replacement.",
        )
    if "mechanism" in support or "supportive" in support:
        return (
            "mechanism_support_only_not_training_positive",
            False,
            "This supports chemistry or biology, but not an exact positive label.",
        )
    return (
        "not_training_positive_in_round6",
        False,
        "Evidence is missing, not exact, or not rechecked enough for training.",
    )


def build_positive_manifest(lit: pd.DataFrame, decisions: pd.DataFrame) -> pd.DataFrame:
    out = lit.copy()
    out["sb"] = out["pair_key"].astype(str).str.split("__").str[0]
    out["pb"] = out["pair_key"].astype(str).str.split("__").str[1]
    dec_cols = ["module", "sb", "pb", "round5_decision", "rule_sources", "rule_ec_numbers", "rule_ids"]
    out = out.merge(decisions[dec_cols], on=["sb", "pb"], how="left")
    status = out.apply(_positive_status, axis=1, result_type="expand")
    out["round6_positive_status"] = status[0]
    out["training_allowed_round6"] = status[1]
    out["training_blocker"] = status[2]
    out["required_missing_fields"] = "exact_table_or_text_extraction;compound_mapping_confirmation;source_reaction_or_rule_provenance;holdout_check;chinese_name"
    out["label_type_round6"] = "positive_candidate_not_yet_training_positive"
    return out


def build_negative_policy() -> pd.DataFrame:
    rows = [
        {
            "label_type": "hard_decoy",
            "definition": "Same-substrate candidate generated by a rule but not present as a curated positive for that substrate.",
            "allowed_for_training": True,
            "allowed_for_external_eval": False,
            "required_fields": "substrate;generated_product;generating_rule_id;candidate_generation_run;positive_pool_version",
            "prohibited_claim": "Do not claim the product is biologically impossible or experimentally negative.",
            "recommended_use": "Learning-to-rank contrast within the same substrate.",
        },
        {
            "label_type": "conditional_negative",
            "definition": "Defined organism, strain, community, or genome-scale model lacks pathway/gene/reaction support under a stated context.",
            "allowed_for_training": "only_in_context_model",
            "allowed_for_external_eval": "only_for_same_context",
            "required_fields": "organism_or_model;context;missing_gene_or_reaction;database_version;confidence_reason",
            "prohibited_claim": "Do not generalize this absence to all gut microbiomes.",
            "recommended_use": "Context-aware ranking or organism-specific product design.",
        },
        {
            "label_type": "assay_negative",
            "definition": "A paper reports no conversion under explicit assay conditions.",
            "allowed_for_training": True,
            "allowed_for_external_eval": True,
            "required_fields": "PMID_or_DOI;assay_conditions;substrate;product_tested;organism_or_enzyme;detection_limit",
            "prohibited_claim": "Do not apply outside the reported assay conditions.",
            "recommended_use": "Highest-quality negative label when available.",
        },
        {
            "label_type": "unlabeled",
            "definition": "No curated positive or negative evidence.",
            "allowed_for_training": False,
            "allowed_for_external_eval": False,
            "required_fields": "compound_pair;search_scope;database_versions",
            "prohibited_claim": "Never treat absence from search as a true negative.",
            "recommended_use": "PU-learning pool, future curation queue, or web-app low-confidence candidate.",
        },
        {
            "label_type": "prohibited_not_found_negative",
            "definition": "A pair is labeled negative only because it was not found in PubMed/Rhea/VMH/other databases.",
            "allowed_for_training": False,
            "allowed_for_external_eval": False,
            "required_fields": "not_allowed",
            "prohibited_claim": "This label type is forbidden.",
            "recommended_use": "Do not use.",
        },
    ]
    return pd.DataFrame(rows)


def build_pipeline_queue(rule_manifest: pd.DataFrame, positive_manifest: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for _, row in rule_manifest.iterrows():
        rows.append(
            {
                "queue_type": "rule_import",
                "module": row.get("module", ""),
                "substrate_name": row.get("substrate_name", ""),
                "product_name": row.get("product_name", ""),
                "status": row.get("round6_rule_import_status", ""),
                "priority": "high" if row.get("round6_rule_import_status") == "candidate_after_source_reaction_check" else "medium",
                "next_action": row.get("next_rule_action", ""),
                "source_ids": row.get("evidence_source_ids", ""),
            }
        )
    for _, row in positive_manifest.iterrows():
        rows.append(
            {
                "queue_type": "positive_curation",
                "module": row.get("module", ""),
                "substrate_name": row.get("substrate_name", ""),
                "product_name": row.get("product_name", ""),
                "status": row.get("round6_positive_status", ""),
                "priority": "high" if "eligible" in str(row.get("round6_positive_status", "")) else "medium",
                "next_action": row.get("next_action", ""),
                "source_ids": row.get("source_id", ""),
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description="Build Round 6 production-gated training manifests.")
    parser.add_argument("--curation-dir", type=Path, default=CURATION)
    args = parser.parse_args()

    decisions = pd.read_csv(args.curation_dir / "rule_recovery_decisions_round5.csv")
    recovery = pd.read_csv(args.curation_dir / "master_rule_recovery_candidates_round5.csv")
    lit = pd.read_csv(args.curation_dir / "literature_verification_round5.csv")

    rule_manifest = build_rule_import_manifest(decisions, recovery, lit)
    positive_manifest = build_positive_manifest(lit, decisions)
    negative_policy = build_negative_policy()
    queue = build_pipeline_queue(rule_manifest, positive_manifest)

    outputs = {
        "rule_import_manifest_round6.csv": rule_manifest,
        "positive_sample_expansion_manifest_round6.csv": positive_manifest,
        "negative_label_policy_round6.csv": negative_policy,
        "round6_pipeline_queue.csv": queue,
    }
    args.curation_dir.mkdir(parents=True, exist_ok=True)
    for name, df in outputs.items():
        df.to_csv(args.curation_dir / name, index=False)

    print(
        {
            "rule_import_rows": len(rule_manifest),
            "rule_import_status": rule_manifest["round6_rule_import_status"].value_counts(dropna=False).to_dict(),
            "positive_rows": len(positive_manifest),
            "positive_status": positive_manifest["round6_positive_status"].value_counts(dropna=False).to_dict(),
            "negative_policy_rows": len(negative_policy),
            "queue_rows": len(queue),
        }
    )
    for name in outputs:
        print(f"wrote {args.curation_dir / name}")


if __name__ == "__main__":
    main()
