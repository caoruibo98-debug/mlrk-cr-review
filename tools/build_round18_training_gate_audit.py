#!/usr/bin/env python
"""Round18 training-readiness gate audit.

Round17 identified priority reaction families and external leads. Round18 turns
those leads into production gates: what is allowed, what is blocked, and which
verification step is missing before a positive or negative sample can be used.

This script is audit-only. It does not modify training data, model code, rules,
or candidate-generation outputs.
"""

from __future__ import annotations

import html
import json
import re
from pathlib import Path

import pandas as pd


REPO = Path(__file__).resolve().parents[1]
CURATION = REPO / "data" / "curation"
DOCS = REPO / "docs" / "reviews"

ROUND16_ROUTE = CURATION / "round16_clean_path_trace.csv"
ROUND17_COVERAGE = CURATION / "round17_reaction_family_coverage_matrix.csv"
ROUND17_POSITIVE = CURATION / "round17_positive_expansion_queue.csv"
ROUND17_NEGATIVE = CURATION / "round17_negative_sampling_logic.csv"
ROUND17_RHEA = CURATION / "round17_rhea_query_triage.csv"
ROUND17_PUBMED = CURATION / "round17_pubmed_query_triage.csv"
ROUND17_GITHUB = CURATION / "round17_github_external_coverage_contract.csv"
CHEBI_RAW = CURATION / "chebi_round18_raw"

REACTIONS = REPO / "outputs" / "modular" / "ltr" / "reactions.parquet"
CLEAN_FULL = REPO / "outputs" / "modular" / "ltr" / "clean_candidates_full.parquet"
CLEAN = REPO / "outputs" / "modular" / "ltr" / "clean_candidates.parquet"

OUT_READINESS = CURATION / "round18_training_readiness_gate.csv"
OUT_P0 = CURATION / "round18_p0_route_repair_worklist.csv"
OUT_POS = CURATION / "round18_positive_candidate_gate.csv"
OUT_NEG = CURATION / "round18_negative_evidence_gate.csv"
OUT_SCHEMA = CURATION / "round18_training_sample_schema_contract.csv"
OUT_NEXT = CURATION / "round18_next_iteration_workplan.csv"
OUT_CHEBI = CURATION / "round18_chebi_structure_status.csv"
OUT_DOC = DOCS / "training_gate_round18.md"


def read_json(path: Path) -> dict:
    if not path.exists():
        return {}
    return json.loads(path.read_text(encoding="utf-8"))


def clean_text(value: object) -> str:
    text = "" if value is None else str(value)
    text = re.sub(r"<[^>]+>", "", text)
    return html.unescape(re.sub(r"\s+", " ", text).strip())


def block1(value: object) -> str:
    text = clean_text(value)
    return text.split("-", 1)[0][:14] if text else ""


def add_pair_key(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out["pair_key"] = out["sb"].astype(str) + "__" + out["pb"].astype(str)
    return out


def parse_chebi_status() -> pd.DataFrame:
    rows = []
    seen_ids: set[str] = set()
    for path in sorted(CHEBI_RAW.glob("CHEBI_*.json")):
        chebi_id = path.stem.replace("CHEBI_", "")
        seen_ids.add(chebi_id)
        raw = read_json(path)
        structure = raw.get("default_structure") or {}
        chem = raw.get("chemical_data") or {}
        smiles = structure.get("smiles", "")
        inchikey = structure.get("standard_inchi_key", "")
        rows.append(
            {
                "round": 18,
                "chebi_id": f"CHEBI:{chebi_id}",
                "name": clean_text(raw.get("name", "")),
                "ascii_name": clean_text(raw.get("ascii_name", "")),
                "stars": raw.get("stars", ""),
                "formula": chem.get("formula", ""),
                "charge": chem.get("charge", ""),
                "smiles": smiles,
                "inchikey": inchikey,
                "inchikey_block1": block1(inchikey),
                "structure_status": "mapped" if smiles and inchikey else "missing_structure",
                "raw_file": str(path.relative_to(REPO)),
            }
        )
    if "58899" not in seen_ids:
        rows.append(
            {
                "round": 18,
                "chebi_id": "CHEBI:58899",
                "name": "",
                "ascii_name": "",
                "stars": "",
                "formula": "",
                "charge": "",
                "smiles": "",
                "inchikey": "",
                "inchikey_block1": "",
                "structure_status": "lookup_failed_or_missing_raw",
                "raw_file": "",
            }
        )
    return pd.DataFrame(rows)


def build_p0_worklist() -> pd.DataFrame:
    route = pd.read_csv(ROUND16_ROUTE)
    p0_families = set(
        pd.read_csv(ROUND17_COVERAGE)
        .query("priority_level == 'P0'")["family"]
        .astype(str)
    )
    rows = []
    for rec in route[route["priority_family"].isin(p0_families)].itertuples(index=False):
        if rec.loss_stage_round16 == "present_in_final_clean":
            gate = "pass_existing_final_clean"
            action = "provenance enrichment only"
        elif rec.loss_stage_round16 == "target_lost_between_clean_full_and_clean":
            gate = "blocked_final_clean_filter_or_variant"
            action = "rebuild or explain final clean filtering so every validated true product for the substrate is retained or intentionally held out"
        elif rec.loss_stage_round16 == "absent_from_clean_full_but_fullrule_hit":
            gate = "blocked_generator_rule_pool_route"
            action = "wire source-traceable rule into deployment-like clean generator and dry-run target substrate"
        else:
            gate = "blocked_unknown_route"
            action = "trace local positive path before retraining"
        rows.append(
            {
                "round": 18,
                "source_id": rec.source_id,
                "priority_family": rec.priority_family,
                "pair_key": rec.pair_key,
                "substrate_name": rec.substrate_name,
                "product_name": rec.product_name,
                "loss_stage_round16": rec.loss_stage_round16,
                "reactions_pair_rows": rec.reactions_pair_rows,
                "clean_full_pair_rows": rec.clean_full_pair_rows,
                "clean_pair_rows": rec.clean_pair_rows,
                "fullrule_pair_hit": rec.fullrule_pair_hit,
                "route_gate_status": gate,
                "training_allowed_round18": False,
                "required_repair": action,
            }
        )
    return pd.DataFrame(rows)


def candidate_exactness(row: pd.Series) -> tuple[str, str, str]:
    source = str(row.get("source_type", ""))
    status = str(row.get("candidate_status", ""))
    text = str(row.get("candidate_text", ""))
    if source == "Rhea" and status == "approved_exact_or_specific_reaction_lead":
        return (
            "database_equation_present",
            "partial_pass",
            "Rhea gives an approved balanced equation; still needs participant structure mapping and local de-dup.",
        )
    if source == "Rhea" and status == "generic_rule_template_context":
        return (
            "generic_reaction_template",
            "fail_for_positive_label",
            "Generic participants such as 'an alcohol' or 'a glucuronoside' are rule context, not a specific positive pair.",
        )
    if source == "PubMed" and status == "enzyme_or_mechanism_lead":
        return (
            "paper_mechanism_lead",
            "needs_manual_extraction",
            "The title suggests an enzyme/mechanism; full paper or abstract must expose exact substrate-product pairs.",
        )
    if source == "PubMed":
        if any(k in text.lower() for k in ["review", "axis", "levels", "biomarker", "supplementation"]):
            return (
                "context_or_association",
                "fail_for_positive_label",
                "This is biological context unless exact conversion evidence is extracted.",
            )
        return (
            "paper_search_lead",
            "needs_manual_extraction",
            "Search hit needs manual reaction extraction before it can become a candidate label.",
        )
    return "unknown", "fail_for_positive_label", "No trusted candidate type was identified."


def build_positive_gate(chebi: pd.DataFrame) -> pd.DataFrame:
    positive = pd.read_csv(ROUND17_POSITIVE)
    rhea = pd.read_csv(ROUND17_RHEA)
    chebi_mapped = chebi["structure_status"].eq("mapped").sum()
    chebi_total = len(chebi)
    rows = []
    for rec in positive.itertuples(index=False):
        s = pd.Series(rec._asdict())
        exactness, gate_status, note = candidate_exactness(s)
        source_id = str(s["source_id"])
        related_rhea = rhea[rhea["rhea_id"].astype(str).eq(source_id)] if source_id.startswith("RHEA:") else pd.DataFrame()
        chebi_ready = ""
        if not related_rhea.empty:
            chebi_ready = "round18_chebi_cache_available_partial"
        elif source_id.startswith("PMID:"):
            chebi_ready = "not_applicable_until_exact_pair_extracted"
        else:
            chebi_ready = "not_checked"
        rows.append(
            {
                "round": 18,
                "priority_family": s["priority_family"],
                "source_type": s["source_type"],
                "source_id": source_id,
                "candidate_status_round17": s["candidate_status"],
                "candidate_text": s["candidate_text"],
                "exactness_level": exactness,
                "exactness_gate_status": gate_status,
                "structure_gate_status": chebi_ready,
                "local_dedup_gate_status": "not_started",
                "generator_route_gate_status": "not_started",
                "split_leakage_gate_status": "not_started",
                "training_allowed_round18": False,
                "why_blocked": note,
                "chebi_round18_cache_summary": f"{chebi_mapped}/{chebi_total} cached ChEBI records have structures; failed/missing records must be rechecked.",
                "raw_file": s["raw_file"],
            }
        )
    return pd.DataFrame(rows)


def build_negative_gate() -> pd.DataFrame:
    neg = pd.read_csv(ROUND17_NEGATIVE)
    rows = []
    for rec in neg.itertuples(index=False):
        if rec.negative_type == "hard_decoy_rule_generated_nontruth":
            gate = "allowed_as_ranking_decoy_only_after_leakage_check"
            status = "conditional"
        else:
            gate = "blocked_until_explicit_no_conversion_source"
            status = "blocked"
        rows.append(
            {
                "round": 18,
                "priority_family": rec.priority_family,
                "negative_type": rec.negative_type,
                "negative_gate_status": status,
                "allowed_training_scope": gate,
                "required_evidence": rec.source_requirements,
                "forbidden_shortcut": rec.forbidden_shortcut,
                "anti_cheat_check": "report decoy source separately from assay negatives; never pool them as the same label type",
            }
        )
    rows.extend(
        [
            {
                "round": 18,
                "priority_family": "all",
                "negative_type": "generator_miss",
                "negative_gate_status": "forbidden",
                "allowed_training_scope": "none",
                "required_evidence": "A generator miss only says the current rule path failed.",
                "forbidden_shortcut": "Do not label a true-world absence from model non-generation.",
                "anti_cheat_check": "Separate generator recall errors from biological negative evidence.",
            },
            {
                "round": 18,
                "priority_family": "all",
                "negative_type": "database_absence",
                "negative_gate_status": "forbidden",
                "allowed_training_scope": "none",
                "required_evidence": "Absence from Rhea/PubMed/VMH/MicrobeRX is not a no-conversion assay.",
                "forbidden_shortcut": "Do not use missing database records as negative labels.",
                "anti_cheat_check": "Require explicit no-conversion text before assay-negative status.",
            },
        ]
    )
    return pd.DataFrame(rows)


def build_readiness_gate(p0: pd.DataFrame, positive: pd.DataFrame, negative: pd.DataFrame) -> pd.DataFrame:
    coverage = pd.read_csv(ROUND17_COVERAGE)
    github = pd.read_csv(ROUND17_GITHUB)
    rows = []
    for rec in coverage.itertuples(index=False):
        fam = rec.family
        route_blockers = p0[p0["priority_family"].eq(fam)]
        pos = positive[positive["priority_family"].eq(fam)]
        neg = negative[negative["priority_family"].eq(fam)]
        if not route_blockers.empty:
            route_gate = "blocked"
        elif rec.dominant_gap_type == "gold_test_visibility_gap":
            route_gate = "pass_for_route_but_eval_gap"
        else:
            route_gate = "not_primary"
        source_gate = (
            "needs_manual_exact_pair_extraction"
            if not pos.empty
            else "no_round18_external_candidate"
        )
        negative_gate = (
            "ranking_decoys_conditional; assay_negatives_missing"
            if not neg.empty
            else "not_defined"
        )
        if route_gate == "blocked":
            decision = "do_not_retrain_fix_route_first"
        elif rec.priority_level == "P1":
            decision = "do_not_retrain_build_family_balanced_gold_eval"
        else:
            decision = "do_not_retrain_monitor_or_split_subfamily"
        rows.append(
            {
                "round": 18,
                "priority_family": fam,
                "priority_level": rec.priority_level,
                "dominant_gap_type": rec.dominant_gap_type,
                "clean_pairs": rec.clean_pairs,
                "clean_gold_pairs": rec.clean_gold_pairs,
                "route_gate": route_gate,
                "source_exactness_gate": source_gate,
                "structure_mapping_gate": "partial_chebi_cache_for_round17_rhea_leads" if not pos.empty else "not_started",
                "generator_route_gate": "blocked_or_not_checked",
                "negative_gate": negative_gate,
                "external_model_contract": "; ".join(github["repo_or_dataset"].astype(str).tolist()),
                "training_allowed_round18": False,
                "decision": decision,
                "next_action": rec.next_action,
            }
        )
    return pd.DataFrame(rows)


def build_schema_contract() -> pd.DataFrame:
    fields = [
        ("record_id", "machine", "stable unique row id", "required"),
        ("label_type", "machine", "positive / hard_decoy / assay_negative", "required"),
        ("negative_type", "machine", "hard_decoy_rule_generated_nontruth / explicit_assay_no_conversion / blank for positive", "required_for_negative"),
        ("reaction_family", "machine", "controlled family name from coverage matrix", "required"),
        ("module", "machine", "A/B/C/D routing module", "required"),
        ("substrate_name_en", "human", "English compound name", "required"),
        ("substrate_name_zh", "human", "Chinese compound name or blank with translation status", "recommended_for_web_app"),
        ("substrate_food_alias", "human", "food/common-source term for user-facing lookup", "recommended_for_web_app"),
        ("substrate_smiles", "machine", "canonical SMILES", "required"),
        ("substrate_inchikey", "machine", "full InChIKey", "required"),
        ("substrate_inchikey_block1", "machine", "first InChIKey block used by current model", "required"),
        ("product_name_en", "human", "English product name", "required"),
        ("product_name_zh", "human", "Chinese product name or blank with translation status", "recommended_for_web_app"),
        ("product_smiles", "machine", "canonical SMILES", "required"),
        ("product_inchikey", "machine", "full InChIKey", "required"),
        ("product_inchikey_block1", "machine", "first InChIKey block used by current model", "required"),
        ("reaction_smiles", "machine", "substrate>>product or source-native reaction SMILES", "required_before_model_export"),
        ("ec_number", "evidence", "full EC or class if available", "required_when_claimed"),
        ("source_database", "evidence", "Rhea/ChEBI/PubChem/BRENDA/MetaNetX/VMH/MicrobeRX/PubMed", "required"),
        ("source_id", "evidence", "RHEA id, PMID, DOI, database accession", "required"),
        ("evidence_sentence", "evidence", "short extracted text supporting conversion/no-conversion", "required_for_gold_or_assay_negative"),
        ("microbe_or_strain", "evidence", "organism/community/strain", "required_when_microbiome_specific"),
        ("assay_condition", "evidence", "condition/time/analytical method", "required_for_assay_negative"),
        ("rule_id", "generator", "RetroRules/Rhea/overlay source-traceable rule id", "required_for_rule_generated_candidate"),
        ("rule_smarts", "generator", "source-backed SMARTS, if promoted", "required_for_rule_promotion"),
        ("generator_hit", "generator", "whether deployment-like generator can produce product", "required_before_training"),
        ("split_group", "evaluation", "scaffold/product-disjoint group", "required_before_eval"),
        ("training_decision", "governance", "allowed/blocked and reason", "required"),
    ]
    return pd.DataFrame(
        [
            {
                "round": 18,
                "field_name": name,
                "field_audience": audience,
                "definition": definition,
                "requirement_level": level,
            }
            for name, audience, definition, level in fields
        ]
    )


def build_next_workplan(readiness: pd.DataFrame) -> pd.DataFrame:
    rows = []
    priority_rows = readiness[readiness["priority_level"].isin(["P0", "P1"])]
    for rec in priority_rows.itertuples(index=False):
        if rec.route_gate == "blocked":
            work = "route_repair"
            success = "target positives reach final clean_candidates or documented holdout table"
        elif rec.source_exactness_gate == "needs_manual_exact_pair_extraction":
            work = "manual_exact_pair_extraction"
            success = "candidate has exact substrate/product, structure mapping, local de-dup, and generator route status"
        else:
            work = "family_balanced_eval_design"
            success = "family has sufficient gold/silver evaluation rows with product/scaffold leakage controls"
        rows.append(
            {
                "round": 18,
                "priority_family": rec.priority_family,
                "next_iteration_work": work,
                "why_this_next": rec.decision,
                "success_criterion": success,
                "recommended_branch": f"round19-{rec.priority_family[:36].replace('_', '-')}",
                "do_not_do": "do not import labels or retrain until the gate success criterion is met",
            }
        )
    return pd.DataFrame(rows)


def write_doc(
    readiness: pd.DataFrame,
    p0: pd.DataFrame,
    pos: pd.DataFrame,
    neg: pd.DataFrame,
    chebi: pd.DataFrame,
) -> None:
    ready_counts = readiness.groupby(["decision", "training_allowed_round18"]).size().to_dict()
    pos_counts = pos.groupby(["exactness_gate_status", "source_type"]).size().to_dict() if not pos.empty else {}
    neg_counts = neg.groupby(["negative_gate_status", "negative_type"]).size().to_dict() if not neg.empty else {}
    chebi_counts = chebi["structure_status"].value_counts().to_dict()
    top = readiness[readiness["priority_level"].isin(["P0", "P1"])][
        ["priority_family", "priority_level", "route_gate", "source_exactness_gate", "negative_gate", "decision"]
    ]

    doc = f"""# Round18 Training Gate Audit

## Working conclusion

No new sample is allowed into training yet. Round18 turns Round17 leads into gates and keeps `training_allowed_round18=False` across the board.

The main reason is not lack of raw rows. It is missing production-grade evidence routing:

1. P0 families still have route blockers.
2. P1 families need family-balanced gold evaluation before retraining.
3. External Rhea/PubMed leads need exact-pair extraction, ChEBI/PubChem mapping, local de-duplication, generator-route checks, and leakage-safe split assignment.
4. Negative labels remain the biggest risk area: hard decoys and assay negatives must not be mixed.

## P0/P1 gate snapshot

```text
{top.to_string(index=False)}
```

## Gate counts

Training decisions:

```text
{json.dumps({str(k): v for k, v in ready_counts.items()}, ensure_ascii=False, indent=2)}
```

Positive candidate exactness:

```text
{json.dumps({str(k): v for k, v in pos_counts.items()}, ensure_ascii=False, indent=2)}
```

Negative evidence:

```text
{json.dumps({str(k): v for k, v in neg_counts.items()}, ensure_ascii=False, indent=2)}
```

ChEBI Round18 structure cache:

```text
{json.dumps(chebi_counts, ensure_ascii=False, indent=2)}
```

## Next production step

Round19 should not retrain. It should pick one P0 route family and make a dry-run fix:

- bile-acid deconjugation: explain or repair final clean filtering so validated multiple true products are not silently lost.
- hydroxycinnamate reduction/hydrolysis: wire a source-traceable rule path for chlorogenate -> trans-caffeate and rerun a target-only clean generation dry run.

Only after those route gates pass should we promote P1 external leads into exact-pair curation.

## Written artifacts

- `data/curation/round18_training_readiness_gate.csv`
- `data/curation/round18_p0_route_repair_worklist.csv`
- `data/curation/round18_positive_candidate_gate.csv`
- `data/curation/round18_negative_evidence_gate.csv`
- `data/curation/round18_training_sample_schema_contract.csv`
- `data/curation/round18_next_iteration_workplan.csv`
- `data/curation/round18_chebi_structure_status.csv`
"""
    DOCS.mkdir(parents=True, exist_ok=True)
    OUT_DOC.write_text(doc, encoding="utf-8")


def main() -> None:
    CURATION.mkdir(parents=True, exist_ok=True)
    DOCS.mkdir(parents=True, exist_ok=True)

    chebi = parse_chebi_status()
    p0 = build_p0_worklist()
    pos = build_positive_gate(chebi)
    neg = build_negative_gate()
    readiness = build_readiness_gate(p0, pos, neg)
    schema = build_schema_contract()
    next_plan = build_next_workplan(readiness)

    chebi.to_csv(OUT_CHEBI, index=False)
    p0.to_csv(OUT_P0, index=False)
    pos.to_csv(OUT_POS, index=False)
    neg.to_csv(OUT_NEG, index=False)
    readiness.to_csv(OUT_READINESS, index=False)
    schema.to_csv(OUT_SCHEMA, index=False)
    next_plan.to_csv(OUT_NEXT, index=False)
    write_doc(readiness, p0, pos, neg, chebi)

    print(f"wrote {OUT_READINESS.relative_to(REPO)} rows={len(readiness)}")
    print(f"wrote {OUT_P0.relative_to(REPO)} rows={len(p0)}")
    print(f"wrote {OUT_POS.relative_to(REPO)} rows={len(pos)}")
    print(f"wrote {OUT_NEG.relative_to(REPO)} rows={len(neg)}")
    print(f"wrote {OUT_SCHEMA.relative_to(REPO)} rows={len(schema)}")
    print(f"wrote {OUT_DOC.relative_to(REPO)}")


if __name__ == "__main__":
    main()
