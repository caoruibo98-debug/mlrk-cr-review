#!/usr/bin/env python
"""Build Round7 source-reaction and sample expansion audit manifests.

Round7 does not add training labels. It checks whether recovered rules and
candidate positive pairs can be traced to existing biochemical databases or
literature before any future import.
"""

from __future__ import annotations

import csv
import re
from collections import defaultdict
from pathlib import Path

import pandas as pd


REPO = Path(__file__).resolve().parents[1]
CURATION = REPO / "data" / "curation"
DATASET = Path(r"D:\CRB\Food models\Dataset")
MNX = DATASET / "metanetx" / "2025-09-11"
RETRORULES = DATASET / "retrorules"


ROUND6_RULES = CURATION / "rule_import_manifest_round6.csv"
ROUND6_POSITIVES = CURATION / "positive_sample_expansion_manifest_round6.csv"
ROUND6_NEGATIVES = CURATION / "negative_label_policy_round6.csv"
ROUND7_SOURCE_CHECK = CURATION / "source_reaction_check_round7.csv"
ROUND7_EVIDENCE_QUEUE = CURATION / "exact_evidence_extraction_queue_round7.csv"
ROUND7_NEGATIVE_PLAN = CURATION / "negative_sample_acquisition_plan_round7.csv"
ROUND7_EXTERNAL_MODELS = CURATION / "external_model_reaction_coverage_round7.csv"


def split_semicolon(value: object) -> list[str]:
    if value is None or pd.isna(value):
        return []
    return [x.strip() for x in str(value).split(";") if x.strip()]


def first_mnxr(legacy_id: object) -> tuple[str, str]:
    """Return base MNXR id and optional template compound id."""
    if legacy_id is None or pd.isna(legacy_id):
        return "", ""
    text = str(legacy_id).strip()
    match = re.match(r"^(MNXR\d+)(?:_(MNXM\d+))?$", text)
    if not match:
        return "", ""
    return match.group(1), match.group(2) or ""


def strip_rule_prefix(rule_id: object) -> str:
    if rule_id is None or pd.isna(rule_id):
        return ""
    text = str(rule_id).strip()
    return text.replace("retrorules_", "", 1)


def read_mnx_tables() -> tuple[pd.DataFrame, pd.DataFrame, dict[str, list[str]]]:
    prop = pd.read_csv(
        MNX / "reac_prop_2025-09-11.tsv",
        sep="\t",
        comment="#",
        header=None,
        names=["mnxr", "equation", "description", "balance", "ec", "source"],
        dtype=str,
        keep_default_na=False,
    )
    xref = pd.read_csv(
        MNX / "reac_xref_2025-09-11.tsv",
        sep="\t",
        comment="#",
        header=None,
        names=["xref", "mnxr", "xref_description"],
        dtype=str,
        keep_default_na=False,
    )
    depr = pd.read_csv(
        MNX / "reac_depr_2025-09-11.tsv",
        sep="\t",
        comment="#",
        header=None,
        names=["deprecated_mnxr", "current_mnxr", "version"],
        dtype=str,
        keep_default_na=False,
    )
    depr_map: dict[str, list[str]] = defaultdict(list)
    for row in depr.itertuples(index=False):
        if row.deprecated_mnxr and row.current_mnxr:
            depr_map[row.deprecated_mnxr].append(row.current_mnxr)
    return prop, xref, depr_map


def aggregate_mnx(
    base_mnxr: str,
    prop: pd.DataFrame,
    xref: pd.DataFrame,
    depr_map: dict[str, list[str]],
) -> dict[str, str | bool | int]:
    if not base_mnxr:
        return {
            "resolved_mnxr_ids": "",
            "deprecation_status": "no_mnxr_legacy_id",
            "metanetx_found": False,
            "mnx_equations": "",
            "mnx_descriptions": "",
            "mnx_ecs": "",
            "xref_count": 0,
            "rhea_xrefs": "",
            "kegg_xrefs": "",
            "metacyc_xrefs": "",
            "seed_xrefs": "",
        }

    current_ids = [base_mnxr]
    deprecation_status = "current_or_unmapped"
    if base_mnxr not in set(prop["mnxr"]) and base_mnxr in depr_map:
        current_ids = sorted(set(depr_map[base_mnxr]))
        deprecation_status = "deprecated_mapped_to_current"

    p = prop[prop["mnxr"].isin(current_ids)].copy()
    x = xref[xref["mnxr"].isin(current_ids)].copy()

    def join_unique(series: pd.Series, limit: int = 12) -> str:
        values = [str(v) for v in series.tolist() if str(v).strip()]
        return ";".join(sorted(set(values))[:limit])

    def join_xref(prefix: str) -> str:
        rows = x[x["xref"].str.lower().str.startswith(prefix.lower())]
        return join_unique(rows["xref"], limit=20)

    return {
        "resolved_mnxr_ids": ";".join(current_ids),
        "deprecation_status": deprecation_status,
        "metanetx_found": not p.empty,
        "mnx_equations": join_unique(p["equation"], limit=6),
        "mnx_descriptions": join_unique(p["description"], limit=6),
        "mnx_ecs": join_unique(p["ec"], limit=12),
        "xref_count": int(len(x)),
        "rhea_xrefs": join_xref("rhea:") or join_xref("rhea.comp:"),
        "kegg_xrefs": join_xref("kegg.reaction:"),
        "metacyc_xrefs": join_xref("metacyc.reaction:"),
        "seed_xrefs": join_xref("seed.reaction:"),
    }


def source_status(row: pd.Series, mnx_info: dict[str, str | bool | int], rr_info: pd.Series | None) -> tuple[str, str]:
    rule_source = str(row.get("rule_source", "")).lower()
    exact_support = str(row.get("exact_reaction_support", "")).lower()
    has_mnx = bool(mnx_info["metanetx_found"])
    has_curated_xref = any(
        str(mnx_info.get(col, "")).strip()
        for col in ["rhea_xrefs", "kegg_xrefs", "metacyc_xrefs"]
    )

    if "retrorules" in rule_source:
        if has_mnx and has_curated_xref:
            if "mechanism_support_only" in exact_support:
                return (
                    "source_reaction_verified_but_pair_not_exact_positive",
                    "do_not_train_positive; keep as mechanism-support or candidate-generation evidence",
                )
            return (
                "source_reaction_candidate_verified_not_imported",
                "manual_direction_license_and_exact_pair_check_before_import",
            )
        if has_mnx:
            return (
                "metanetx_reaction_found_without_high_value_xref",
                "manual_database_review_before_import",
            )
        if str(mnx_info["deprecation_status"]) == "deprecated_mapped_to_current":
            return (
                "deprecated_mapping_unresolved_in_current_reac_prop",
                "manual_metanetx_deprecation_review",
            )
        return ("source_reaction_missing", "do_not_import_until_source_reaction_found")

    if "microberx" in rule_source:
        return (
            "microberx_rule_needs_source_recovery",
            "replace_or_annotate_with_ec_database_or_exact_literature_rule_before_training",
        )

    return ("unknown_rule_source", "do_not_import_until_rule_source_is_normalized")


def build_source_check() -> pd.DataFrame:
    rules = pd.read_csv(ROUND6_RULES)
    prop, xref, depr_map = read_mnx_tables()
    rr = pd.read_csv(RETRORULES / "retrorules_index.csv", dtype=str, keep_default_na=False)
    rr_by_id = rr.set_index("rule_id", drop=False)

    rows: list[dict[str, object]] = []
    for row in rules.itertuples(index=False):
        s = pd.Series(row._asdict())
        base_mnxr, template_mnxm = first_mnxr(s.get("rule_legacy_id"))
        mnx_info = aggregate_mnx(base_mnxr, prop, xref, depr_map)
        clean_rule_id = strip_rule_prefix(s.get("rule_id"))
        rr_info = rr_by_id.loc[clean_rule_id] if clean_rule_id in rr_by_id.index else None
        status, recommendation = source_status(s, mnx_info, rr_info)

        rows.append(
            {
                "round": 7,
                "module": s.get("module", ""),
                "sb": s.get("sb", ""),
                "pb": s.get("pb", ""),
                "substrate_name": s.get("substrate_name", ""),
                "product_name": s.get("product_name", ""),
                "round6_rule_import_status": s.get("round6_rule_import_status", ""),
                "round6_training_allowed": s.get("training_allowed_round6", ""),
                "round6_import_blocker": s.get("import_blocker", ""),
                "rule_id": s.get("rule_id", ""),
                "rule_source": s.get("rule_source", ""),
                "rule_ec_number": s.get("rule_ec_number", ""),
                "rule_ec_class": s.get("rule_ec_class", ""),
                "rule_legacy_id": s.get("rule_legacy_id", ""),
                "base_mnxr_id": base_mnxr,
                "template_mnxm_id": template_mnxm,
                "resolved_mnxr_ids": mnx_info["resolved_mnxr_ids"],
                "deprecation_status": mnx_info["deprecation_status"],
                "metanetx_found": mnx_info["metanetx_found"],
                "mnx_equations": mnx_info["mnx_equations"],
                "mnx_descriptions": mnx_info["mnx_descriptions"],
                "mnx_ecs": mnx_info["mnx_ecs"],
                "xref_count": mnx_info["xref_count"],
                "rhea_xrefs": mnx_info["rhea_xrefs"],
                "kegg_xrefs": mnx_info["kegg_xrefs"],
                "metacyc_xrefs": mnx_info["metacyc_xrefs"],
                "seed_xrefs": mnx_info["seed_xrefs"],
                "retrorules_index_found": rr_info is not None,
                "retrorules_direction": "" if rr_info is None else rr_info.get("reaction_direction", ""),
                "retrorules_relative_direction": "" if rr_info is None else rr_info.get("rule_relative_direction", ""),
                "retrorules_usage": "" if rr_info is None else rr_info.get("rule_usage", ""),
                "rule_score": s.get("rule_score", ""),
                "rule_diameter": s.get("rule_diameter", ""),
                "exact_reaction_support": s.get("exact_reaction_support", ""),
                "evidence_source_ids": s.get("evidence_source_ids", ""),
                "source_reaction_check_status": status,
                "round7_import_recommendation": recommendation,
                "import_allowed_round7": False,
                "training_allowed_round7": False,
                "reason_training_still_blocked": (
                    "Round7 verifies source reaction traceability only; exact biological pair, "
                    "direction, stereochemistry, license, holdout, and split checks remain required."
                ),
            }
        )
    return pd.DataFrame(rows)


def build_evidence_queue(source_check: pd.DataFrame) -> pd.DataFrame:
    positives = pd.read_csv(ROUND6_POSITIVES)
    status_by_pair = {
        (r.sb, r.pb): r.source_reaction_check_status for r in source_check.itertuples(index=False)
    }
    import_by_pair = {
        (r.sb, r.pb): r.round7_import_recommendation for r in source_check.itertuples(index=False)
    }
    rows = []
    for row in positives.itertuples(index=False):
        pair = (row.sb, row.pb)
        source_ids = split_semicolon(row.source_id)
        rows.append(
            {
                "round": 7,
                "pair_key": row.pair_key,
                "module": row.module,
                "gap_family": row.gap_family,
                "substrate_name": row.substrate_name,
                "product_name": row.product_name,
                "source_db": row.source_db,
                "source_id": row.source_id,
                "source_title": row.source_title,
                "evidence_type": row.evidence_type,
                "round6_positive_status": row.round6_positive_status,
                "source_reaction_check_status": status_by_pair.get(pair, "no_rule_source_check_row"),
                "round7_import_recommendation": import_by_pair.get(pair, "not_applicable"),
                "exact_text_or_table_needed": True,
                "compound_mapping_needed": True,
                "organism_or_enzyme_needed": True,
                "direction_needed": True,
                "stereochemistry_needed": True,
                "assay_context_needed": True,
                "negative_evidence_needed": (
                    "only_if_paper_reports_tested_no_conversion; absence from databases is not a negative"
                ),
                "pmid_or_doi_tokens": ";".join(source_ids),
                "round7_training_decision": "not_training_positive_yet",
                "minimum_fields_before_training": (
                    "substrate_inchikey;product_inchikey;substrate_smiles;product_smiles;"
                    "exact_source_sentence_or_table;PMID_or_DOI;organism_or_enzyme;"
                    "direction;assay_context;evidence_tier;holdout_guard;license_note"
                ),
            }
        )
    return pd.DataFrame(rows)


def build_negative_plan() -> pd.DataFrame:
    policy = pd.read_csv(ROUND6_NEGATIVES)
    plan_rows = [
        {
            "round": 7,
            "negative_source_type": "same_substrate_hard_decoy",
            "evidence_basis": "rule-generated candidate product for a substrate that is not in curated positives",
            "can_train": True,
            "claim_boundary": "ranking contrast only; not a biological impossibility claim",
            "required_fields": "sb;candidate_pb;candidate_smiles;generating_rule_id;generation_run_id;positive_pool_version",
            "best_use": "LTR within-substrate contrast",
            "risk": "unknown positives become mislabeled negatives",
            "control": "PU-learning framing; down-weight decoys; keep high-evidence positives out of decoy pool",
        },
        {
            "round": 7,
            "negative_source_type": "assay_negative",
            "evidence_basis": "paper explicitly reports no conversion under stated assay conditions",
            "can_train": True,
            "claim_boundary": "negative only under reported organism/enzyme/assay context",
            "required_fields": "PMID_or_DOI;substrate;product_tested;organism_or_enzyme;conditions;detection_limit",
            "best_use": "highest-quality negative labels and external evaluation",
            "risk": "rare in literature; context may not transfer",
            "control": "store context fields and never generalize across gut microbiome",
        },
        {
            "round": 7,
            "negative_source_type": "genome_or_model_absence",
            "evidence_basis": "defined strain/GEM lacks enzyme/reaction/pathway support",
            "can_train": "context_only",
            "claim_boundary": "conditional negative for a named organism/model and database version",
            "required_fields": "organism_or_model;database_version;missing_gene_or_reaction;confidence_reason",
            "best_use": "organism-aware ranking or product design constraints",
            "risk": "annotation incompleteness is not true absence",
            "control": "separate from global product negatives; use only with organism context",
        },
        {
            "round": 7,
            "negative_source_type": "database_absence",
            "evidence_basis": "pair not found in Rhea/MetaNetX/PubMed/local positives",
            "can_train": False,
            "claim_boundary": "not evidence of non-reaction",
            "required_fields": "search_query;database_versions;date;not_found_scope",
            "best_use": "triage for manual review only",
            "risk": "creates false negatives and hides unknown biology",
            "control": "do not train as y=0; label as unlabeled",
        },
    ]
    out = pd.DataFrame(plan_rows)
    if not policy.empty:
        out["round6_policy_reference"] = ";".join(policy["label_type"].astype(str).tolist())
    return out


def build_external_model_coverage() -> pd.DataFrame:
    rows = [
        {
            "round": 7,
            "source_or_model": "ECREACT / RXN biocatalysis model",
            "github_or_public_url": "https://github.com/rxn4chemistry/biocatalysis-model",
            "reaction_scope": "enzyme-catalysed reactions across all 7 EC first-level classes",
            "data_sources": "Rhea;BRENDA;PathBank;MetaNetX",
            "reported_scale_or_signal": "published model reports forward top-1 accuracy 49.6%; dataset merged public biochemical reactions",
            "what_it_teaches_this_project": "Our four modules do not yet show balanced EC-class/family coverage or a credible external split; use EC/source-token fields as evidence, not as automatic positives.",
            "can_directly_import_samples": "only if local license and molecule mapping are checked",
            "candidate_fields_to_align": "rxn_smiles;ec;source;substrate_smiles;product_smiles;source_database;license",
        },
        {
            "round": 7,
            "source_or_model": "RetroRules / RetroPathRL",
            "github_or_public_url": "https://github.com/brsynth/RetroPathRL",
            "reaction_scope": "mono-component reaction SMARTS templates from biochemical reactions; rule radius/diameter controls specificity",
            "data_sources": "RetroRules;MetaNetX;Rhea-derived reaction identifiers",
            "reported_scale_or_signal": "RetroRules 2026 public article reports >1M templates and thousands of EC numbers",
            "what_it_teaches_this_project": "Rules are candidate generators. A rule firing is not a positive label unless its specific pair has database/literature evidence.",
            "can_directly_import_samples": "rules maybe; labels no",
            "candidate_fields_to_align": "Rule_ID;Reaction_ID;Rule_SMARTS;EC;Rule_usage;source_reaction_id;direction",
        },
        {
            "round": 7,
            "source_or_model": "gapseq",
            "github_or_public_url": "https://github.com/jotech/gapseq",
            "reaction_scope": "bacterial metabolic pathways, transporters, genome-scale networks, gap filling",
            "data_sources": "MNXref;MetaCyc;ModelSEED;KEGG;BRENDA;UniProt;TCDB;BiGG",
            "reported_scale_or_signal": "paper validates enzyme activity and metabolic phenotype predictions against experimental/literature data",
            "what_it_teaches_this_project": "Production needs organism/genome-aware evidence and context-specific negatives, not only molecule-pair ranking.",
            "can_directly_import_samples": "not as generic positives; useful for organism-context evidence and conditional negatives",
            "candidate_fields_to_align": "organism;gene_or_EC;reaction_id;database_version;pathway;confidence",
        },
        {
            "round": 7,
            "source_or_model": "GutBugDB / GutBug",
            "github_or_public_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC11810598/",
            "reaction_scope": "gut microbiome-mediated biotransformation potential for biotic and xenobiotic molecules",
            "data_sources": "gut bacterial genomes; EC assignments; PubChem molecules",
            "reported_scale_or_signal": "database built from hundreds of gut genomes and >300k enzyme assignments",
            "what_it_teaches_this_project": "Good source for organism/enzyme context, but predictions are not wet-lab positives.",
            "can_directly_import_samples": "only as weak context/evidence, not gold positives",
            "candidate_fields_to_align": "compound_pubchem;predicted_EC;organism;protein;prediction_score",
        },
        {
            "round": 7,
            "source_or_model": "gutSMASH",
            "github_or_public_url": "https://github.com/victoriapascal/gutsmash",
            "reaction_scope": "anaerobic gut microbial metabolic gene clusters and specialized primary metabolism",
            "data_sources": "curated MGC detection rules; known clusters; genome annotations",
            "reported_scale_or_signal": "rules validated with curated dataset in publication",
            "what_it_teaches_this_project": "This is gene-cluster potential, not direct substrate-product labels; useful for biological plausibility and future product design.",
            "can_directly_import_samples": "no direct pair labels",
            "candidate_fields_to_align": "genome;MGC_type;gene_function;known_cluster_match;predicted_metabolic_function",
        },
    ]
    return pd.DataFrame(rows)


def main() -> None:
    CURATION.mkdir(parents=True, exist_ok=True)
    source_check = build_source_check()
    evidence_queue = build_evidence_queue(source_check)
    negative_plan = build_negative_plan()
    external_models = build_external_model_coverage()

    source_check.to_csv(ROUND7_SOURCE_CHECK, index=False, quoting=csv.QUOTE_MINIMAL)
    evidence_queue.to_csv(ROUND7_EVIDENCE_QUEUE, index=False, quoting=csv.QUOTE_MINIMAL)
    negative_plan.to_csv(ROUND7_NEGATIVE_PLAN, index=False, quoting=csv.QUOTE_MINIMAL)
    external_models.to_csv(ROUND7_EXTERNAL_MODELS, index=False, quoting=csv.QUOTE_MINIMAL)

    print(f"wrote {ROUND7_SOURCE_CHECK} rows={len(source_check)}")
    print(source_check["source_reaction_check_status"].value_counts(dropna=False).to_string())
    print(f"wrote {ROUND7_EVIDENCE_QUEUE} rows={len(evidence_queue)}")
    print(f"wrote {ROUND7_NEGATIVE_PLAN} rows={len(negative_plan)}")
    print(f"wrote {ROUND7_EXTERNAL_MODELS} rows={len(external_models)}")


if __name__ == "__main__":
    main()
