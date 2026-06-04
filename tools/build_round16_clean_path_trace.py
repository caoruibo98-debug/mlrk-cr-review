#!/usr/bin/env python
"""Round16 clean-candidate path and source-backed expansion audit.

This round answers a narrow production question:

Are the current blockers mainly missing reaction families, or are known
source-backed positives being lost between the positive pool, rule generation,
and the final clean candidate table?

The script is audit-only. It writes review CSVs and a markdown note, but does
not change training labels, generator code, or model artifacts.
"""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger


REPO = Path(__file__).resolve().parents[1]
SRC = REPO / "src"
MODULAR = REPO / "modular"
sys.path.insert(0, str(SRC))
sys.path.insert(0, str(MODULAR))

import kio  # noqa: E402
from generate_candidates import run_reactants  # noqa: E402
from ltr_candidates import load_module_rules  # noqa: E402


RDLogger.DisableLog("rdApp.*")

CURATION = REPO / "data" / "curation"
DOCS = REPO / "docs" / "reviews"
RHEA_RAW = CURATION / "rhea_round16_raw"
NCBI_RAW = CURATION / "ncbi_round16_raw"

ROUND15_PAIRS = CURATION / "round15_candidate_pair_gate.csv"
REACTIONS = REPO / "outputs" / "modular" / "ltr" / "reactions.parquet"
CLEAN_FULL = REPO / "outputs" / "modular" / "ltr" / "clean_candidates_full.parquet"
CLEAN = REPO / "outputs" / "modular" / "ltr" / "clean_candidates.parquet"
FULLRULE = REPO / "outputs" / "gen_gap" / "fullrule_hit_miss.parquet"
METRICS = REPO / "outputs" / "modular" / "ltr" / "clean2_metrics.csv"

OUT_LOCAL = CURATION / "round16_clean_path_trace.csv"
OUT_RULES = CURATION / "round16_rule_sampling_sensitivity.csv"
OUT_RHEA = CURATION / "round16_rhea_query_triage.csv"
OUT_PUBMED = CURATION / "round16_pubmed_query_triage.csv"
OUT_MODELS = CURATION / "round16_external_model_position.csv"
OUT_DECISIONS = CURATION / "round16_training_decisions.csv"
OUT_PIPELINE = CURATION / "round16_collect_verify_modify_pipeline.csv"
OUT_DOC = DOCS / "clean_path_trace_round16.md"


MODEL_STANDARDS = [
    {
        "reference_model_or_dataset": "ECREACT / RXN for biocatalysis",
        "source_url": "https://github.com/rxn4chemistry/biocatalysis-model",
        "coverage_claim": "62,222 unique enzymatic reaction-EC combinations aggregated from Rhea, BRENDA, PathBank, and MetaNetX.",
        "fields_or_representation": "rxn_smiles, ec, source; reaction SMILES with EC tokens at multiple EC levels.",
        "lesson_for_foodgut": "Keep EC/source fields explicit and split by product/reaction to avoid product leakage.",
    },
    {
        "reference_model_or_dataset": "EnzymeMap",
        "source_url": "https://github.com/hesther/enzymemap",
        "coverage_claim": "Large curated and atom-mapped enzymatic reaction dataset with processed_reactions.csv.gz and correction/validation scripts.",
        "fields_or_representation": "Balanced/atom-mapped enzymatic reactions; EC-oriented scripts; processed reaction CSV.",
        "lesson_for_foodgut": "Require correction, validation, and atom/structure checks before adding reactions as strict positives.",
    },
    {
        "reference_model_or_dataset": "gapseq",
        "source_url": "https://github.com/jotech/gapseq",
        "coverage_claim": "Curated bacterial metabolism database with about 15,150 reactions and 8,446 metabolites; benchmarked against phenotype data.",
        "fields_or_representation": "Reaction database, transporter database, sequence-homology evidence, gap-filling trace flags.",
        "lesson_for_foodgut": "Production needs trace flags explaining why every reaction was admitted, not only a generated product score.",
    },
    {
        "reference_model_or_dataset": "FoodGut current branch",
        "source_url": "local:outputs/modular/ltr",
        "coverage_claim": "10,550 unique positive reactions; 80 gold positives; final clean candidates contain 2,139 positive rows; clean2 metrics use n_pts=10 per module/method.",
        "fields_or_representation": "InChIKey block1 sb/pb pairs, canonical SMILES, module, tier, EC/source evidence features.",
        "lesson_for_foodgut": "The current total sample count is not the main bottleneck; strict gold/test breadth and generator-route recall are.",
    },
]


QUERY_FAMILIES = {
    "bile_salt_hydrolase": "bile_acid_deconjugation_and_lipid_context",
    "taurocholate_hydrolase": "bile_acid_deconjugation_and_lipid_context",
    "glycocholate_hydrolase": "bile_acid_deconjugation_and_lipid_context",
    "chlorogenate_hydrolase": "hydroxycinnamate_reduction_and_hydrolysis",
    "caffeate_reduction": "hydroxycinnamate_reduction_and_hydrolysis",
    "quercetin_degradation": "polyphenol_ring_fission",
    "catechin_gut_microbiota": "polyphenol_ring_fission",
    "daidzein_reductase": "isoflavone_reductive_and_glycoside_metabolism",
    "equol": "isoflavone_reductive_and_glycoside_metabolism",
    "urolithin_dehydroxylation": "urolithin_dehydroxylation",
    "fucosyllactose_hydrolase": "glycoside_and_hmo_hydrolysis",
    "pinoresinol_reductase": "lignan_redox_and_deglycosylation",
}

PUBMED_QUERY_FAMILIES = {
    "gut_microbiota_bile_salt_hydrolase_glycocholate_cholate": "bile_acid_deconjugation_and_lipid_context",
    "gut_microbiota_chlorogenic_acid_caffeic_acid_quinic_acid_hydrolysis": "hydroxycinnamate_reduction_and_hydrolysis",
    "gut_microbiota_quercetin_ring_fission_dihydroxyphenylacetic_acid": "polyphenol_ring_fission",
    "catechin_valerolactone_gut_microbiota_reaction": "polyphenol_ring_fission",
    "urolithin_c_urolithin_a_dehydroxylation_gut_microbiota": "urolithin_dehydroxylation",
    "2_fucosyllactose_fucose_bifidobacterium_fucosidase": "glycoside_and_hmo_hydrolysis",
    "daidzein_equol_non_producer_no_conversion_bacteria": "isoflavone_reductive_and_glycoside_metabolism",
    "pinoresinol_lariciresinol_human_intestinal_bacteria": "lignan_redox_and_deglycosylation",
}


def read_json(path: Path) -> dict:
    if not path.exists():
        return {}
    return json.loads(path.read_text(encoding="utf-8"))


def csv_join(values: list[object] | set[object]) -> str:
    return ";".join(sorted({str(v) for v in values if str(v) and str(v) != "nan"}))


def table_with_pair_key(df: pd.DataFrame) -> pd.DataFrame:
    if df.empty:
        return df.copy()
    out = df.copy()
    out["pair_key"] = out["sb"].astype(str) + "__" + out["pb"].astype(str)
    return out


def source_counts(df: pd.DataFrame, target: pd.Series) -> dict[str, object]:
    if df.empty:
        return {
            "substrate_rows": 0,
            "pair_rows": 0,
            "positive_rows": 0,
            "positive_products_for_substrate": "",
            "candidate_products_for_substrate": 0,
        }
    sub = df[(df["module"].astype(str) == str(target["module"])) & (df["sb"].astype(str) == str(target["sb"]))]
    pair = sub[sub["pb"].astype(str) == str(target["pb"])]
    pos = sub[sub["y"].eq(1)] if "y" in sub.columns else sub
    return {
        "substrate_rows": len(sub),
        "pair_rows": len(pair),
        "positive_rows": len(pair[pair["y"].eq(1)]) if "y" in pair.columns else len(pair),
        "positive_products_for_substrate": csv_join(pos["pb"].tolist()) if "pb" in pos.columns else "",
        "candidate_products_for_substrate": sub["pb"].nunique() if "pb" in sub.columns else 0,
    }


def infer_loss_stage(row: dict[str, object]) -> tuple[str, str]:
    if row["clean_pair_rows"] > 0:
        return "present_in_final_clean", "The pair is already available to the final ranker."
    if row["clean_full_pair_rows"] > 0 and row["clean_pair_rows"] == 0:
        return (
            "target_lost_between_clean_full_and_clean",
            "The positive reaches clean_candidates_full but is absent from final clean_candidates; trace final filtering or generation variant.",
        )
    if row["reactions_pair_rows"] > 0 and row["clean_full_pair_rows"] == 0 and row["fullrule_pair_hit"]:
        return (
            "absent_from_clean_full_but_fullrule_hit",
            "The positive exists locally and full-rule diagnosis can hit it, but deployment-like clean generation cannot.",
        )
    if row["reactions_pair_rows"] > 0 and row["clean_full_pair_rows"] == 0:
        return (
            "absent_from_clean_full",
            "The positive exists locally but does not reach rule-generated clean candidates.",
        )
    return "not_in_local_positive_pool", "Candidate would need source-backed import before generator/ranker evaluation."


def build_local_trace() -> pd.DataFrame:
    pairs = pd.read_csv(ROUND15_PAIRS)
    rxn = table_with_pair_key(pd.read_parquet(REACTIONS))
    clean_full = table_with_pair_key(pd.read_parquet(CLEAN_FULL))
    clean = table_with_pair_key(pd.read_parquet(CLEAN))
    fullrule = table_with_pair_key(pd.read_parquet(FULLRULE)) if FULLRULE.exists() else pd.DataFrame()

    rows = []
    for target in pairs.itertuples(index=False):
        t = pd.Series(target._asdict())
        r = source_counts(rxn, t)
        cf = source_counts(clean_full, t)
        c = source_counts(clean, t)
        fr_pair = (
            fullrule[fullrule["pair_key"].astype(str).eq(str(t["pair_key"]))]
            if not fullrule.empty and "pair_key" in fullrule.columns
            else pd.DataFrame()
        )
        row = {
            "round": 16,
            "source_id": t["source_id"],
            "priority_family": t["priority_family"],
            "module": t["module"],
            "substrate_name": t["substrate_name"],
            "product_name": t["product_name"],
            "pair_key": t["pair_key"],
            "reactions_pair_rows": r["pair_rows"],
            "reactions_substrate_rows": r["substrate_rows"],
            "reactions_positive_products_for_substrate": r["positive_products_for_substrate"],
            "clean_full_pair_rows": cf["pair_rows"],
            "clean_full_substrate_rows": cf["substrate_rows"],
            "clean_full_positive_rows_for_pair": cf["positive_rows"],
            "clean_full_positive_products_for_substrate": cf["positive_products_for_substrate"],
            "clean_pair_rows": c["pair_rows"],
            "clean_substrate_rows": c["substrate_rows"],
            "clean_positive_rows_for_pair": c["positive_rows"],
            "clean_positive_products_for_substrate": c["positive_products_for_substrate"],
            "fullrule_pair_rows": len(fr_pair),
            "fullrule_pair_hit": bool(fr_pair["hit"].astype(bool).any()) if "hit" in fr_pair.columns else False,
        }
        row["clean_full_lost_positive_products"] = csv_join(
            set(str(row["clean_full_positive_products_for_substrate"]).split(";"))
            - set(str(row["clean_positive_products_for_substrate"]).split(";"))
        )
        stage, reason = infer_loss_stage(row)
        row["loss_stage_round16"] = stage
        row["production_blocker_reason"] = reason
        rows.append(row)
    return pd.DataFrame(rows)


def build_rule_sensitivity(targets: pd.DataFrame) -> pd.DataFrame:
    rxn = table_with_pair_key(pd.read_parquet(REACTIONS))
    test_sizes = [400, 600, 800, 1200, 2000, 10000]
    rows = []
    for n_rules in test_sizes:
        for require_ec in [False, True]:
            rng = np.random.default_rng(0)
            rules = load_module_rules(rng, n_rules)
            if require_ec:
                rules = {m: [(s, e) for s, e in vals if int(e) != 0] for m, vals in rules.items()}
            for t in targets.itertuples(index=False):
                m = str(t.module)
                local = rxn[(rxn["module"].astype(str).eq(m)) & (rxn["sb"].astype(str).eq(str(t.sb)))]
                smiles = t.generator_input_substrate_smiles or ""
                if not smiles and not local.empty:
                    smiles = str(local["substrate_smiles"].iloc[0])
                mol = Chem.MolFromSmiles(smiles) if smiles else None
                generated: dict[str, list[int]] = {}
                if mol is not None:
                    molh = Chem.AddHs(mol)
                    for smarts, ecc in rules.get(m, []):
                        for psmi in run_reactants(smarts, mol, molh):
                            pb = kio.inchikey_block1(kio.smiles_to_inchikey(psmi))
                            if not pb or pb == str(t.sb):
                                continue
                            generated.setdefault(pb, []).append(int(ecc))
                hit_ec = sorted(set(generated.get(str(t.pb), [])))
                rows.append(
                    {
                        "round": 16,
                        "pair_key": t.pair_key,
                        "source_id": t.source_id,
                        "module": m,
                        "priority_family": t.priority_family,
                        "n_rules_requested": n_rules,
                        "require_ec": bool(require_ec),
                        "module_rules_loaded": len(rules.get(m, [])),
                        "generated_product_count": len(generated),
                        "target_generated": str(t.pb) in generated,
                        "target_rule_count": len(generated.get(str(t.pb), [])),
                        "target_hit_ec_classes": csv_join(hit_ec),
                        "rule_path_interpretation": rule_path_interpretation(str(t.pair_key), bool(require_ec), hit_ec),
                    }
                )
    return pd.DataFrame(rows)


def rule_path_interpretation(pair_key: str, require_ec: bool, hit_ec: list[int]) -> str:
    if not hit_ec:
        return "target not generated under this rule sample/filter"
    if require_ec:
        return "target generated with EC-bearing rules"
    if set(hit_ec) == {0}:
        return "target generated only by EC=0 rules; production must recover provenance or explicitly allow unannotated rules"
    return "target generated by at least one rule path in this setting"


def build_rhea_triage() -> pd.DataFrame:
    rows = []
    for path in sorted(RHEA_RAW.glob("*.json")):
        query_slug = path.stem
        raw = read_json(path)
        results = raw.get("results", [])
        if not results:
            rows.append(
                {
                    "round": 16,
                    "query_slug": query_slug,
                    "priority_family": QUERY_FAMILIES.get(query_slug, "unknown"),
                    "rhea_id": "",
                    "equation": "",
                    "status": "no_rhea_result",
                    "balanced": "",
                    "triage_decision": "needs_literature_or_other_database",
                    "training_use": "not_a_sample",
                    "raw_file": str(path.relative_to(REPO)),
                }
            )
            continue
        for rec in results:
            rid = f"RHEA:{rec.get('id')}"
            equation = str(rec.get("equation", ""))
            if rid in {"RHEA:16309", "RHEA:19353", "RHEA:20689", "RHEA:69683"}:
                decision = "source_backed_existing_positive"
                training_use = "do_not_duplicate; use for route repair or provenance"
            elif rec.get("status") == "approved" and rec.get("balanced") is True:
                decision = "approved_biochemical_context"
                training_use = "candidate_only_after exact substrate/product mapping and local de-dup"
            else:
                decision = "weak_context"
                training_use = "not_a_sample"
            rows.append(
                {
                    "round": 16,
                    "query_slug": query_slug,
                    "priority_family": QUERY_FAMILIES.get(query_slug, "unknown"),
                    "rhea_id": rid,
                    "equation": equation,
                    "status": rec.get("status", ""),
                    "balanced": rec.get("balanced", ""),
                    "triage_decision": decision,
                    "training_use": training_use,
                    "raw_file": str(path.relative_to(REPO)),
                }
            )
    return pd.DataFrame(rows)


def build_pubmed_triage() -> pd.DataFrame:
    summary = read_json(NCBI_RAW / "round16_pubmed_esummary.json")
    result = summary.get("result", {})
    title_by_pmid = {
        pmid: result.get(pmid, {}).get("title", "")
        for pmid in result.get("uids", [])
    }
    journal_by_pmid = {
        pmid: result.get(pmid, {}).get("fulljournalname", result.get(pmid, {}).get("source", ""))
        for pmid in result.get("uids", [])
    }
    rows = []
    for path in sorted(NCBI_RAW.glob("*_esearch.json")):
        query_slug = re.sub(r"_esearch$", "", path.stem)
        raw = read_json(path)
        ids = raw.get("esearchresult", {}).get("idlist", [])
        if not ids:
            rows.append(
                {
                    "round": 16,
                    "query_slug": query_slug,
                    "priority_family": PUBMED_QUERY_FAMILIES.get(query_slug, "unknown"),
                    "pmid": "",
                    "title": "",
                    "journal": "",
                    "triage_decision": "search_miss",
                    "training_use": "not_a_sample",
                    "raw_file": str(path.relative_to(REPO)),
                }
            )
            continue
        for pmid in ids:
            title = title_by_pmid.get(pmid, "")
            text = title.lower()
            if any(term in text for term in ["conversion", "metabolism", "hydrolase", "fucosidase", "transform", "catechin-converting"]):
                decision = "literature_lead_for_manual_exact_pair_review"
                training_use = "candidate_only_after exact substrate/product/evidence extraction"
            elif any(term in text for term in ["non producer", "no conversion", "antibiotics", "interindividual"]):
                decision = "possible_assay_negative_or_cohort_context"
                training_use = "negative_candidate_only_if no-conversion assay conditions are explicit"
            else:
                decision = "context_only"
                training_use = "not_a_sample"
            rows.append(
                {
                    "round": 16,
                    "query_slug": query_slug,
                    "priority_family": PUBMED_QUERY_FAMILIES.get(query_slug, "unknown"),
                    "pmid": pmid,
                    "title": title,
                    "journal": journal_by_pmid.get(pmid, ""),
                    "triage_decision": decision,
                    "training_use": training_use,
                    "raw_file": str(path.relative_to(REPO)),
                }
            )
    return pd.DataFrame(rows)


def build_decisions(local: pd.DataFrame, rules: pd.DataFrame, rhea: pd.DataFrame, pubmed: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for rec in local.itertuples(index=False):
        rs = rules[rules["pair_key"].eq(rec.pair_key)]
        generated_any = bool(rs["target_generated"].any()) if not rs.empty else False
        generated_with_ec = bool(rs[rs["require_ec"]]["target_generated"].any()) if not rs.empty else False
        if rec.loss_stage_round16 == "present_in_final_clean":
            decision = "do_not_import_existing_positive_already_reaches_ranker"
            action = "provenance enrichment only"
        elif rec.loss_stage_round16 == "target_lost_between_clean_full_and_clean":
            decision = "do_not_import_existing_positive_fix_final_clean_filtering"
            action = "trace why final clean keeps only part of the true-product set for this substrate"
        elif rec.loss_stage_round16 == "absent_from_clean_full_but_fullrule_hit":
            decision = "do_not_import_existing_positive_wire_source_traceable_rule"
            action = "promote or overlay the source-traceable rule path and rerun clean candidate generation"
        else:
            decision = "not_enough_evidence_for_training"
            action = "verify source, structure, and generator route before import"
        rows.append(
            {
                "round": 16,
                "source_id": rec.source_id,
                "pair_key": rec.pair_key,
                "priority_family": rec.priority_family,
                "training_allowed_round16": False,
                "decision": decision,
                "local_loss_stage": rec.loss_stage_round16,
                "target_generated_any_round16": generated_any,
                "target_generated_with_ec_round16": generated_with_ec,
                "recommended_next_action": action,
            }
        )

    for family in sorted(set(rhea["priority_family"]) | set(pubmed["priority_family"])):
        if family in {"unknown"}:
            continue
        rhea_hits = rhea[(rhea["priority_family"].eq(family)) & rhea["rhea_id"].astype(str).ne("")]
        pubmed_hits = pubmed[(pubmed["priority_family"].eq(family)) & pubmed["pmid"].astype(str).ne("")]
        if not rhea_hits.empty or not pubmed_hits.empty:
            rows.append(
                {
                    "round": 16,
                    "source_id": "external_screen",
                    "pair_key": "",
                    "priority_family": family,
                    "training_allowed_round16": False,
                    "decision": "external_leads_need_collect_verify_modify_cycle",
                    "local_loss_stage": "",
                    "target_generated_any_round16": "",
                    "target_generated_with_ec_round16": "",
                    "recommended_next_action": "extract exact substrate/product pairs, map ChEBI/PubChem/InChIKey, de-dup against reactions, then run generator-route gate",
                }
            )
    return pd.DataFrame(rows)


def build_pipeline() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "step_order": 1,
                "pipeline_step": "collect",
                "required_fields": "source_id; query; database; candidate_equation; substrate_name; product_name; EC; organism_or_microbe; PMID/DOI",
                "acceptance_gate": "Candidate must come from Rhea/ChEBI/PubChem/BRENDA/MetaNetX/VMH/MicrobeRX or a specific paper, never from model imagination.",
                "output_table": "roundXX_external_query_triage.csv",
            },
            {
                "step_order": 2,
                "pipeline_step": "verify_truth",
                "required_fields": "exact substrate; exact product; reaction direction; single-step flag; structure parse status; InChIKey block1; evidence sentence",
                "acceptance_gate": "Strict positive requires exact substrate-product conversion, not broad pathway context.",
                "output_table": "roundXX_positive_candidate_gate.csv",
            },
            {
                "step_order": 3,
                "pipeline_step": "deduplicate_and_route",
                "required_fields": "raw_pool_pair_rows; reactions_pair_rows; clean_full_pair_rows; clean_pair_rows; fullrule_hit; rule provenance",
                "acceptance_gate": "If already in reactions, do not import; repair generator/clean route first.",
                "output_table": "roundXX_clean_path_trace.csv",
            },
            {
                "step_order": 4,
                "pipeline_step": "modify_branch_only",
                "required_fields": "source-backed rule SMARTS; rule id; EC; Rhea/MetaNetX/RetroRules link; dry-run generated products",
                "acceptance_gate": "Rule promotion is allowed only after source-traceable dry run; no hallucinated reaction rules.",
                "output_table": "roundXX_rule_overlay_dryrun.csv",
            },
            {
                "step_order": 5,
                "pipeline_step": "negative_sampling",
                "required_fields": "negative_type; source; substrate; non-product; assay_condition_or_decoy_generation_rule; leakage check",
                "acceptance_gate": "No-hit is not a biological negative; use hard decoys for ranking, and assay negatives only with explicit no-conversion evidence.",
                "output_table": "roundXX_negative_candidate_gate.csv",
            },
            {
                "step_order": 6,
                "pipeline_step": "summarize_and_freeze",
                "required_fields": "commit; diff; CSV list; metrics caveat; unresolved gaps",
                "acceptance_gate": "Freeze every iteration before retraining; high metrics must be tied to adequate n and family coverage.",
                "output_table": "docs/reviews/roundXX.md",
            },
        ]
    )


def write_doc(local: pd.DataFrame, rules: pd.DataFrame, rhea: pd.DataFrame, pubmed: pd.DataFrame, decisions: pd.DataFrame) -> None:
    metrics = pd.read_csv(METRICS)
    npts = csv_join(metrics["n_pts"].unique().tolist()) if "n_pts" in metrics.columns else ""
    local_counts = local["loss_stage_round16"].value_counts().to_dict()
    rhea_counts = rhea["triage_decision"].value_counts().to_dict() if not rhea.empty else {}
    pubmed_counts = pubmed["triage_decision"].value_counts().to_dict() if not pubmed.empty else {}
    ec_only_hits = rules.groupby(["pair_key", "require_ec"])["target_generated"].any().reset_index()

    doc = f"""# Round16 Clean Candidate Path Trace

## Working conclusion

Ray's suspicion is directionally right, but it needs a sharper split:

1. The model does have reaction-family coverage gaps.
2. The immediate production blocker is even more concrete: source-backed positives already present in `reactions.parquet` can be lost before final `clean_candidates.parquet`.
3. The current high metrics are not production evidence because `clean2_metrics.csv` evaluates tiny panels (`n_pts={npts}`).

## Local path result

Loss-stage counts:

```text
{json.dumps(local_counts, ensure_ascii=False, indent=2)}
```

The four Round15 Rhea candidates are not new labels. Three of them are repair targets:

- bile-acid deconjugation positives reach `clean_candidates_full` but are absent from final clean candidates.
- chlorogenate hydrolysis is in `reactions.parquet` and full-rule diagnosis can hit it, but deployment-like clean generation misses it.
- daidzein glycoside hydrolysis already reaches final clean candidates and should only receive provenance enrichment.

## Rule sensitivity result

The rule-sampling table shows whether a target appears only under broad/unannotated rules or under EC-bearing rules. A target that only appears under EC=0 rules is not safe for production unless rule provenance is recovered.

```text
{ec_only_hits.to_string(index=False)}
```

## External retrieval result

Rhea triage counts:

```text
{json.dumps(rhea_counts, ensure_ascii=False, indent=2)}
```

PubMed triage counts:

```text
{json.dumps(pubmed_counts, ensure_ascii=False, indent=2)}
```

Rhea is strong for exact biochemical reactions such as bile-acid deconjugation and chlorogenate hydrolysis. It is weak for several gut-microbiome transformation families where the evidence lives in primary papers or curated microbiome databases rather than Rhea exact reactions.

## Production implication

Do not retrain yet. First fix route recall and evidence gates:

1. Preserve all true products for substrates that have multiple validated products, or document the final clean filtering rule.
2. Promote only source-traceable rules for known full-rule hits that are absent from deployment-like clean generation.
3. Expand gold/test sets by reaction family, not by random database duplicates.
4. Treat no-hit generator outputs as hard decoys, not biological negatives.
5. Add assay negatives only when a paper/database explicitly reports no conversion under defined microbe/condition/time.

## Written artifacts

- `data/curation/round16_clean_path_trace.csv`
- `data/curation/round16_rule_sampling_sensitivity.csv`
- `data/curation/round16_rhea_query_triage.csv`
- `data/curation/round16_pubmed_query_triage.csv`
- `data/curation/round16_external_model_position.csv`
- `data/curation/round16_training_decisions.csv`
- `data/curation/round16_collect_verify_modify_pipeline.csv`
"""
    DOCS.mkdir(parents=True, exist_ok=True)
    OUT_DOC.write_text(doc, encoding="utf-8")


def main() -> None:
    CURATION.mkdir(parents=True, exist_ok=True)
    DOCS.mkdir(parents=True, exist_ok=True)

    local = build_local_trace()
    targets = pd.read_csv(ROUND15_PAIRS)
    rules = build_rule_sensitivity(targets)
    rhea = build_rhea_triage()
    pubmed = build_pubmed_triage()
    models = pd.DataFrame(MODEL_STANDARDS)
    decisions = build_decisions(local, rules, rhea, pubmed)
    pipeline = build_pipeline()

    local.to_csv(OUT_LOCAL, index=False)
    rules.to_csv(OUT_RULES, index=False)
    rhea.to_csv(OUT_RHEA, index=False)
    pubmed.to_csv(OUT_PUBMED, index=False)
    models.to_csv(OUT_MODELS, index=False)
    decisions.to_csv(OUT_DECISIONS, index=False)
    pipeline.to_csv(OUT_PIPELINE, index=False)
    write_doc(local, rules, rhea, pubmed, decisions)

    print(f"wrote {OUT_LOCAL.relative_to(REPO)} rows={len(local)}")
    print(f"wrote {OUT_RULES.relative_to(REPO)} rows={len(rules)}")
    print(f"wrote {OUT_RHEA.relative_to(REPO)} rows={len(rhea)}")
    print(f"wrote {OUT_PUBMED.relative_to(REPO)} rows={len(pubmed)}")
    print(f"wrote {OUT_DOC.relative_to(REPO)}")


if __name__ == "__main__":
    main()
