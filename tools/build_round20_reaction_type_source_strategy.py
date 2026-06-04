#!/usr/bin/env python
"""Round20 reaction-type gap and source-acquisition strategy.

This is an audit-only script. It writes CSV manifests that separate:

1. local reaction-family retention gaps;
2. external source roles for source-backed positive expansion;
3. negative-sample logic that does not pretend unknown biology is false;
4. the collect -> verify -> modify -> summarize -> CSV -> next-round pipeline.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
from rdkit import RDLogger


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "curation"
DOC = ROOT / "docs" / "reviews" / "reaction_type_gap_and_sample_strategy_round20.md"
LTR_OUT = ROOT / "outputs" / "modular" / "ltr"

sys.path.insert(0, str(ROOT / "modular"))
sys.path.insert(0, str(ROOT / "src"))
import ltr_build  # noqa: E402


RDLogger.DisableLog("rdApp.*")


def pair_key(df: pd.DataFrame) -> pd.Series:
    return df["module"].astype(str) + "__" + df["sb"].astype(str) + "__" + df["pb"].astype(str)


def joined(values: pd.Series) -> str:
    return ";".join(sorted({str(v) for v in values.dropna() if str(v) and str(v) != "nan"}))


def risk_label(row: pd.Series) -> str:
    labels = []
    if row["gold_pairs"] == 0 and row["silver_pairs"] == 0:
        labels.append("no_gold_or_silver_eval_anchor")
    if row["clean_final_pos"] == 0 and row["reactions_pairs"] >= 20:
        labels.append("complete_final_clean_loss")
    elif row["clean_final_retention_vs_reactions"] < 0.05 and row["reactions_pairs"] >= 20:
        labels.append("severe_final_clean_loss")
    if row["clean_full_pos"] > row["clean_final_pos"]:
        labels.append("clean_full_to_final_route_loss")
    if row["ec_supported_pairs"] == 0 and row["reactions_pairs"] >= 20:
        labels.append("missing_ec_mapping")
    if row["rule_supported_pairs"] > 0 and row["ec_supported_pairs"] == 0:
        labels.append("rule_only_not_training_ready")
    return ";".join(labels) if labels else "monitor"


def recommended_source(row: pd.Series) -> str:
    category = str(row["reaction_category"]).lower()
    module = str(row["module"])
    if "flavonoid" in category or "polyphenol" in category or "ring" in category:
        return "Rhea/ChEBI exact search + EnzymeMap/ECReact for EC-class coverage + PubMed manual curation"
    if module == "B" or "glyco" in category or "carbohydrate" in category:
        return "CAZy/CAZac + Rhea/ChEBI + literature-curated glycosidase reactions"
    if "bile" in category or "steroid" in category or "lipid" in category:
        return "Rhea exact reactions + VMH/gapseq microbial route evidence + PubMed BSH/HSDH/bai pathway papers"
    if "amino" in category or module == "C":
        return "Rhea/ECReact/EnzymeMap + gapseq/KEGG/MetaCyc pathway membership"
    return "Rhea first; EnzymeMap/ECReact for mapped enzyme reaction support; RetroRules only as rule provenance"


def build_gap_matrix() -> pd.DataFrame:
    pool = ltr_build.load_pool()
    reactions = pd.read_parquet(LTR_OUT / "reactions.parquet")
    clean_full = pd.read_parquet(LTR_OUT / "clean_candidates_full.parquet")
    clean = pd.read_parquet(LTR_OUT / "clean_candidates.parquet")

    pool = pool.copy()
    pool["_pair_key"] = pair_key(pool)
    reactions = reactions.copy()
    clean_full = clean_full.copy()
    clean = clean.copy()
    reactions["_pair_key"] = pair_key(reactions)
    clean_full["_pair_key"] = pair_key(clean_full)
    clean["_pair_key"] = pair_key(clean)

    pool["tier"] = pool["tier_rank"].map({2: "gold", 1: "silver", 0: "weak"}).fillna("weak")
    pool["_evidence_score"] = pool[
        ["tier_rank", "ev_has_ec", "ev_has_rule", "ev_has_pmid", "ev_microbe_n", "ev_gene_n"]
    ].sum(axis=1)
    meta = (
        pool.sort_values("_evidence_score", ascending=False)
        .drop_duplicates("_pair_key")
        .copy()
    )
    meta["in_reactions"] = meta["_pair_key"].isin(set(reactions["_pair_key"]))
    meta["in_clean_full_pos"] = meta["_pair_key"].isin(set(clean_full.loc[clean_full["y"].eq(1), "_pair_key"]))
    meta["in_clean_pos"] = meta["_pair_key"].isin(set(clean.loc[clean["y"].eq(1), "_pair_key"]))

    group = meta.groupby(["module", "reaction_category"], dropna=False)
    out = group.agg(
        unique_pool_pairs=("_pair_key", "nunique"),
        reactions_pairs=("in_reactions", "sum"),
        clean_full_pos=("in_clean_full_pos", "sum"),
        clean_final_pos=("in_clean_pos", "sum"),
        gold_pairs=("tier", lambda s: int((s == "gold").sum())),
        silver_pairs=("tier", lambda s: int((s == "silver").sum())),
        weak_pairs=("tier", lambda s: int((s == "weak").sum())),
        ec_supported_pairs=("ev_has_ec", "sum"),
        rule_supported_pairs=("ev_has_rule", "sum"),
        pmid_supported_pairs=("ev_has_pmid", "sum"),
        representative_reaction_types=("reaction_type", joined),
        source_origins=("source_origin_type", joined),
        training_use_recommendations=("training_use_recommendation", joined),
    ).reset_index()
    out["route_loss_reactions_to_final"] = out["reactions_pairs"] - out["clean_final_pos"]
    out["clean_final_retention_vs_reactions"] = (
        out["clean_final_pos"] / out["reactions_pairs"].replace(0, pd.NA)
    ).fillna(0).round(4)
    out["reaction_category"] = out["reaction_category"].fillna("missing_reaction_category")
    out["gap_labels_round20"] = out.apply(risk_label, axis=1)
    out["recommended_external_source_path"] = out.apply(recommended_source, axis=1)
    out["training_allowed_round20"] = False
    return out.sort_values(
        ["module", "route_loss_reactions_to_final", "reactions_pairs"],
        ascending=[True, False, False],
    )


def build_external_sources() -> pd.DataFrame:
    rows = [
        {
            "source_name": "Rhea",
            "url": "https://www.rhea-db.org/",
            "source_role": "primary positive verification",
            "use_for_positive": "exact substrate/product/co-product reaction, ChEBI participants, EC, references, UniProt links",
            "use_for_rules": "validate EC provenance and direction before promoting a rule",
            "use_for_negative": "exclude candidates already known as curated reactions",
            "guardrail": "Rhea hit is a reaction truth source, not proof that a specific gut species performs it unless enzyme/microbe evidence is attached",
            "priority": "P0",
        },
        {
            "source_name": "ECReact / rxn4chemistry biocatalysis-model",
            "url": "https://github.com/rxn4chemistry/biocatalysis-model",
            "source_role": "broad enzymatic reaction coverage benchmark",
            "use_for_positive": "EC-labelled reaction SMILES across all seven EC classes after exact structure matching",
            "use_for_rules": "compare our missing EC classes and reaction templates against a published model dataset",
            "use_for_negative": "remove generated decoys that appear as known enzyme reactions",
            "guardrail": "not microbiome-specific; source database and structure mapping must be retained",
            "priority": "P1",
        },
        {
            "source_name": "EnzymeMap",
            "url": "https://github.com/hesther/enzymemap",
            "source_role": "mapped and corrected enzymatic reaction expansion",
            "use_for_positive": "atom-mapped/balanced enzyme reactions when exact participants and EC can be aligned",
            "use_for_rules": "check atom-mapping plausibility for transformations before rule promotion",
            "use_for_negative": "screen decoys against known mapped enzymatic reactions",
            "guardrail": "BRENDA-derived entries still require exact participant and direction checks before import",
            "priority": "P1",
        },
        {
            "source_name": "RetroRules",
            "url": "https://retrorules.org/docs",
            "source_role": "candidate generation rule source",
            "use_for_positive": "not direct labels; only source reactions behind templates can become candidate positives",
            "use_for_rules": "reaction SMARTS, radius/specificity, EC, source reaction, sequence-aware score",
            "use_for_negative": "do not call non-hit templates negative; use as candidate-generation universe only",
            "guardrail": "template existence is not biological occurrence and must not be imported as a positive label by itself",
            "priority": "P1",
        },
        {
            "source_name": "gapseq / ModelSEED / KEGG / MetaCyc",
            "url": "https://github.com/jotech/gapseq",
            "source_role": "microbial pathway and genome-route evidence",
            "use_for_positive": "attach organism/pathway/EC/gene support to reactions already verified elsewhere",
            "use_for_rules": "prioritize rules whose EC/gene evidence exists in gut-relevant microbes",
            "use_for_negative": "conditional negatives from absence of enzyme/pathway under a named organism/genome context",
            "guardrail": "genome-model absence is condition-specific and not a universal biochemical negative",
            "priority": "P1",
        },
        {
            "source_name": "VMH",
            "url": "https://www.vmh.life/",
            "source_role": "gut microbiome reaction and metabolite relevance",
            "use_for_positive": "gut microbe reaction/metabolite/gene context when exact reaction IDs can be mapped",
            "use_for_rules": "connect food/gut metabolites to shared human-gut microbiome nomenclature",
            "use_for_negative": "condition-specific absence or unsupported route flag only",
            "guardrail": "VMH context supports relevance, but exact substrate/product identity still needs reaction-source verification",
            "priority": "P1",
        },
        {
            "source_name": "CAZy / CAZac",
            "url": "https://www.cazy.org/",
            "source_role": "carbohydrate-active enzyme coverage",
            "use_for_positive": "glycosidase/glycosyltransferase/polysaccharide lyase/carbohydrate esterase family support",
            "use_for_rules": "repair module B reaction-type coverage and enzyme-family labels",
            "use_for_negative": "not a negative source; use to avoid falsely rejecting plausible glycan reactions",
            "guardrail": "enzyme family annotation needs substrate specificity before creating exact reaction labels",
            "priority": "P1",
        },
        {
            "source_name": "gutMGene",
            "url": "https://pubmed.ncbi.nlm.nih.gov/39475181",
            "source_role": "microbe-metabolite-gene association context",
            "use_for_positive": "human/mouse gut microbe-metabolite-gene support, not exact reaction labels",
            "use_for_rules": "rank gut relevance and web-app explanatory metadata",
            "use_for_negative": "not suitable as a biochemical negative source",
            "guardrail": "association evidence is not substrate-product conversion evidence",
            "priority": "P2",
        },
    ]
    return pd.DataFrame(rows)


def build_pipeline() -> pd.DataFrame:
    rows = [
        {
            "step_order": 1,
            "stage": "collect",
            "question": "Which local reaction families have low clean-route retention or no gold/silver anchor?",
            "input_artifacts": "round20_reaction_type_gap_matrix.csv; reactions.parquet; clean_candidates.parquet",
            "required_fields": "module; reaction_category; reaction_type; sb; pb; substrate_name; product_name; source_dataset",
            "pass_gate": "gap family selected with explicit route-loss/sample-evidence reason",
            "output_artifact": "roundXX_collection_worklist.csv",
        },
        {
            "step_order": 2,
            "stage": "verify_truth",
            "question": "Is this exact substrate-product reaction present in a curated source?",
            "input_artifacts": "Rhea; ChEBI; EnzymeMap; ECReact; PubMed; UniProt; VMH/gapseq when microbe context is needed",
            "required_fields": "external_reaction_id; equation; participant_chebi_ids; inchikeys; ec; direction; references",
            "pass_gate": "exact participant identity plus curated reaction/source reference; no analogy-only import",
            "output_artifact": "roundXX_verified_positive_candidates.csv",
        },
        {
            "step_order": 3,
            "stage": "modify_manifest_only",
            "question": "Should the sample become training, route repair, evaluation-only, or blocked?",
            "input_artifacts": "verified positive candidates; current reactions/clean route status",
            "required_fields": "in_reactions; in_clean_full; in_clean_final; route_loss_cause; leakage_risk",
            "pass_gate": "training_allowed is false until route/generator/split gates pass",
            "output_artifact": "roundXX_positive_candidate_gate.csv",
        },
        {
            "step_order": 4,
            "stage": "negative_logic",
            "question": "What kind of negative is defensible?",
            "input_artifacts": "generated candidates; Rhea/ECReact/EnzymeMap screens; organism/pathway evidence; assay papers",
            "required_fields": "negative_type; condition_scope; screened_known_reaction_sources; reason_not_positive",
            "pass_gate": "unknown is not treated as false; negative class is named and condition-scoped",
            "output_artifact": "roundXX_negative_candidate_gate.csv",
        },
        {
            "step_order": 5,
            "stage": "route_dryrun",
            "question": "Can the generator produce the verified target with at least one decoy under production-like filters?",
            "input_artifacts": "candidate rules; clean generator; exact target panel",
            "required_fields": "target_generated; target_ec_class; decoy_count; final_clean_eligible; rule_hash",
            "pass_gate": "target_generated=True and decoy_count>0 without disabling anti-cheat filters globally",
            "output_artifact": "roundXX_route_dryrun.csv",
        },
        {
            "step_order": 6,
            "stage": "summarize_and_next_round",
            "question": "What changed, what remains blocked, and what should be collected next?",
            "input_artifacts": "all round CSVs; review markdown; git diff",
            "required_fields": "decision; blocker; next_gate; commit; training_allowed",
            "pass_gate": "CSV plus review doc committed before next round starts",
            "output_artifact": "docs/reviews/roundXX_summary.md",
        },
    ]
    return pd.DataFrame(rows)


def build_negative_logic() -> pd.DataFrame:
    rows = [
        {
            "negative_type": "ranker_decoy_generated_same_substrate",
            "definition": "candidate product generated for the same substrate but not among known positives",
            "acceptable_use": "within-substrate ranking negative / hard decoy",
            "not_allowed_claim": "biochemically impossible reaction",
            "verification_required": "screen against Rhea, ECReact, EnzymeMap, and current positive pool",
            "condition_scope": "unlabeled unless screened; model-training weight should reflect uncertainty",
        },
        {
            "negative_type": "organism_context_absent_pathway",
            "definition": "enzyme/pathway absent in a named organism or genome model under explicit criteria",
            "acceptable_use": "condition-specific microbe-route negative",
            "not_allowed_claim": "universal negative for all gut microbiota",
            "verification_required": "gapseq/ModelSEED/KEGG/MetaCyc route criteria plus genome/protein evidence threshold",
            "condition_scope": "organism/strain/genome-specific",
        },
        {
            "negative_type": "assay_tested_no_conversion",
            "definition": "paper reports substrate tested and no product/low conversion under assay conditions",
            "acceptable_use": "strongest negative, but assay-condition-specific",
            "not_allowed_claim": "reaction can never occur",
            "verification_required": "PubMed paper, organism/enzyme, substrate, condition, detection method, negative sentence",
            "condition_scope": "enzyme/strain/assay-specific",
        },
        {
            "negative_type": "direction_or_cofactor_invalid_decoy",
            "definition": "candidate violates curated direction/cofactor/stoichiometry for the selected context",
            "acceptable_use": "structural/process decoy after Rhea/ChEBI direction and balance check",
            "not_allowed_claim": "false across all biochemical contexts",
            "verification_required": "reaction balance, EC class, direction, cofactor availability",
            "condition_scope": "route-condition-specific",
        },
        {
            "negative_type": "invalid_structure_or_unmapped",
            "definition": "SMILES/InChIKey parse failure, unbalanced participant, R-group, polymer ambiguity",
            "acceptable_use": "exclude from training/evaluation",
            "not_allowed_claim": "negative sample",
            "verification_required": "RDKit parse status and ChEBI/PubChem/HMDB mapping status",
            "condition_scope": "data-quality exclusion",
        },
    ]
    return pd.DataFrame(rows)


def write_doc(gap: pd.DataFrame, sources: pd.DataFrame, pipeline: pd.DataFrame, negatives: pd.DataFrame) -> None:
    DOC.parent.mkdir(parents=True, exist_ok=True)
    severe_all = gap[
        gap["gap_labels_round20"].str.contains("complete_final_clean_loss|severe_final_clean_loss", na=False)
        & gap["reactions_pairs"].ge(20)
    ].copy()
    severe = (
        severe_all.sort_values(
            ["module", "clean_final_retention_vs_reactions", "route_loss_reactions_to_final"],
            ascending=[True, True, False],
        )
        .groupby("module", group_keys=False)
        .head(5)
    )
    text = "# Round20 Reaction-Type Gap And Sample Strategy\n\n"
    text += "## Working conclusion\n\n"
    text += (
        "Ray's diagnosis is mostly right, but the production blocker is more precise than 'sample type is not enough'. "
        "The current model fails because reaction-family coverage, source-backed exactness, generator route retention, "
        "negative-label semantics, and evaluation size are all weak. Reaction-type gaps are one major cause, but route "
        "loss can erase positives that already exist in `reactions.parquet`.\n\n"
    )
    text += "## Highest-risk local families\n\n```text\n"
    cols = [
        "module",
        "reaction_category",
        "reactions_pairs",
        "clean_full_pos",
        "clean_final_pos",
        "gold_pairs",
        "silver_pairs",
        "clean_final_retention_vs_reactions",
        "gap_labels_round20",
    ]
    text += severe[cols].to_string(index=False)
    text += "\n```\n\n"
    text += "## External source roles\n\n```text\n"
    text += sources[["source_name", "source_role", "priority", "guardrail"]].to_string(index=False)
    text += "\n```\n\n"
    text += "## Negative sample logic\n\n```text\n"
    text += negatives[["negative_type", "acceptable_use", "not_allowed_claim", "condition_scope"]].to_string(index=False)
    text += "\n```\n\n"
    text += "## Pipeline\n\n```text\n"
    text += pipeline[["step_order", "stage", "pass_gate", "output_artifact"]].to_string(index=False)
    text += "\n```\n\n"
    text += "## Written artifacts\n\n"
    text += "- `data/curation/round20_reaction_type_gap_matrix.csv`\n"
    text += "- `data/curation/round20_external_training_source_candidates.csv`\n"
    text += "- `data/curation/round20_negative_sample_logic.csv`\n"
    text += "- `data/curation/round20_sample_acquisition_pipeline.csv`\n"
    DOC.write_text(text, encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    gap = build_gap_matrix()
    sources = build_external_sources()
    pipeline = build_pipeline()
    negatives = build_negative_logic()
    outputs = {
        "round20_reaction_type_gap_matrix.csv": gap,
        "round20_external_training_source_candidates.csv": sources,
        "round20_sample_acquisition_pipeline.csv": pipeline,
        "round20_negative_sample_logic.csv": negatives,
    }
    for name, df in outputs.items():
        path = OUT / name
        df.to_csv(path, index=False)
        print(f"wrote {path.relative_to(ROOT)} rows={len(df)}")
    write_doc(gap, sources, pipeline, negatives)
    print(f"wrote {DOC.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
