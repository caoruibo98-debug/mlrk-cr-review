#!/usr/bin/env python
"""Build Round14 sample-expansion review tables.

Round14 answers Ray's current review question:

- Is reaction-family coverage really a blocker?
- Which external databases/repositories support adding source-backed samples?
- Which candidate positives and negatives are safe enough for the next curation
  round, and which are only search leads?

This script does not modify model code, training labels, or candidate generator
outputs. It summarizes local evidence plus raw PubMed/Rhea screens into
auditable CSVs and one markdown review.
"""

from __future__ import annotations

import json
import re
import textwrap
import xml.etree.ElementTree as ET
from collections import defaultdict
from pathlib import Path

import pandas as pd


REPO = Path(__file__).resolve().parents[1]
CURATION = REPO / "data" / "curation"
DOCS = REPO / "docs" / "reviews"
NCBI_RAW = CURATION / "ncbi_round14_raw"
RHEA_RAW = CURATION / "rhea_round14_raw"

ROUND10_GAPS = CURATION / "reaction_family_gap_priorities_round10.csv"
ROUND11_GENERATOR = CURATION / "known_gold_generator_trace_round11.csv"
ROUND11_DECISIONS = CURATION / "round11_training_gate_decisions.csv"
ROUND12_OVERLAY = CURATION / "rule_promotion_overlay_round12.csv"
ROUND13_DRYRUN = CURATION / "rule_overlay_dryrun_round13.csv"
METRICS = REPO / "outputs" / "modular" / "ltr" / "clean2_metrics.csv"

OUT_STANDARDS = CURATION / "round14_external_model_standards.csv"
OUT_PUBMED_QUERY = CURATION / "round14_pubmed_query_audit.csv"
OUT_PUBMED_TRIAGE = CURATION / "round14_pubmed_candidate_triage.csv"
OUT_RHEA = CURATION / "round14_rhea_reaction_triage.csv"
OUT_BLOCKERS = CURATION / "round14_reaction_type_blocker_diagnosis.csv"
OUT_POSITIVE = CURATION / "round14_positive_sample_candidates.csv"
OUT_NEGATIVE = CURATION / "round14_negative_sample_strategy.csv"
OUT_PIPELINE = CURATION / "round14_collect_verify_modify_pipeline.csv"
OUT_DOC = DOCS / "reaction_sample_expansion_round14.md"


PUBMED_QUERY_SPECS = {
    "polyphenol_ring_fission_catechin": {
        "priority_family": "polyphenol_ring_fission",
        "query": "gut microbiota catechin valerolactone phenylvalerolactone exact substrate product",
        "role": "strict_exact_pair_search",
    },
    "polyphenol_ring_fission_catechin_broad": {
        "priority_family": "polyphenol_ring_fission",
        "query": "catechin valerolactone gut microbiota",
        "role": "broad_literature_screen",
    },
    "polyphenol_ring_fission_quercetin": {
        "priority_family": "polyphenol_ring_fission",
        "query": "gut microbiota quercetin dihydroxyphenylacetic acid ring fission",
        "role": "broad_literature_screen",
    },
    "urolithin_c_to_a": {
        "priority_family": "urolithin_dehydroxylation",
        "query": "urolithin C urolithin A gut microbiota dehydroxylation",
        "role": "strict_exact_pair_search",
    },
    "urolithin_c_to_a_broad": {
        "priority_family": "urolithin_dehydroxylation",
        "query": "urolithin C urolithin A",
        "role": "broad_literature_screen",
    },
    "hmo_2fl_fucose": {
        "priority_family": "glycoside_and_hmo_hydrolysis",
        "query": "2-fucosyllactose L-fucose Bifidobacterium infantis hydrolysis",
        "role": "strict_exact_pair_search",
    },
    "hmo_2fl_fucose_broad": {
        "priority_family": "glycoside_and_hmo_hydrolysis",
        "query": "2 fucosyllactose fucose bifidobacterium",
        "role": "broad_literature_screen",
    },
    "hmo_fucosidase_broad": {
        "priority_family": "glycoside_and_hmo_hydrolysis",
        "query": "fucosyllactose fucosidase Bifidobacterium",
        "role": "broad_literature_screen",
    },
    "pinoresinol_lariciresinol": {
        "priority_family": "lignan_redox_and_deglycosylation",
        "query": "pinoresinol lariciresinol human intestinal microflora",
        "role": "known_exact_pair_support",
    },
    "isoxanthohumol_8pn": {
        "priority_family": "prenylflavonoid_o_demethylation",
        "query": "isoxanthohumol 8-prenylnaringenin human intestinal bacteria",
        "role": "known_exact_pair_support",
    },
    "hydroxycinnamate_reduction": {
        "priority_family": "hydroxycinnamate_reduction_and_hydrolysis",
        "query": "gut microbiota caffeic acid dihydrocaffeic acid reduction",
        "role": "broad_literature_screen",
    },
    "isoflavone_equol_negative": {
        "priority_family": "isoflavone_reductive_and_glycoside_metabolism",
        "query": "daidzein equol non-producer no conversion gut bacteria",
        "role": "negative_evidence_search",
    },
    "isoflavone_equol_no_conversion_broad": {
        "priority_family": "isoflavone_reductive_and_glycoside_metabolism",
        "query": "daidzein equol no conversion",
        "role": "negative_evidence_search",
    },
    "isoflavone_equol_nonproducer_broad": {
        "priority_family": "isoflavone_reductive_and_glycoside_metabolism",
        "query": "equol non producer daidzein conversion",
        "role": "negative_evidence_search",
    },
}


RHEA_QUERY_FAMILY = {
    "catechin": "polyphenol_ring_fission",
    "phenylvalerolactone": "polyphenol_ring_fission",
    "quercetin": "polyphenol_ring_fission",
    "urolithin": "urolithin_dehydroxylation",
    "fucosyllactose": "glycoside_and_hmo_hydrolysis",
    "pinoresinol": "lignan_redox_and_deglycosylation",
    "lariciresinol": "lignan_redox_and_deglycosylation",
    "isoxanthohumol": "prenylflavonoid_o_demethylation",
    "8_prenylnaringenin": "prenylflavonoid_o_demethylation",
    "caffeic_acid": "hydroxycinnamate_reduction_and_hydrolysis",
    "chlorogenate": "hydroxycinnamate_reduction_and_hydrolysis",
    "dihydrocaffeic_acid": "hydroxycinnamate_reduction_and_hydrolysis",
    "daidzein": "isoflavone_reductive_and_glycoside_metabolism",
    "equol": "isoflavone_reductive_and_glycoside_metabolism",
    "bile_salt_hydrolase": "bile_acid_deconjugation_and_lipid_context",
}


KNOWN_EXISTING_EXACT_OR_SUPPORT_PMIDS = {
    "31138818": "existing_gold_exact_or_support: 2-FL -> L-fucose",
    "12736449": "existing_gold_exact_or_support: pinoresinol -> lariciresinol",
    "16772450": "existing_gold_exact_or_support: isoxanthohumol -> 8-prenylnaringenin",
    "41797252": "supporting_context: urolithin dehydroxylation",
}


KNOWN_MISS_PAIR_FAMILY = {
    "HHXMEXZVPJFAIJ__RIUPLDUFZCXCHM": "urolithin_dehydroxylation",
    "SNFSYLYCDAVZGP__SHZGCJCMOBCMKK": "glycoside_and_hmo_hydrolysis",
    "HGXBRUKMWQGOIE__MHXCIKYXNYCMHY": "lignan_redox_and_deglycosylation",
    "YKGCBLWILMDSAV__LPEPZZAVFJPLNZ": "prenylflavonoid_o_demethylation",
}


def safe_read_csv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path) if path.exists() else pd.DataFrame()


def clean_text(value: object) -> str:
    return re.sub(r"\s+", " ", "" if value is None else str(value)).strip()


def parse_pubmed_xml(path: Path) -> dict[str, dict[str, str]]:
    if not path.exists():
        return {}
    root = ET.parse(path).getroot()
    records: dict[str, dict[str, str]] = {}
    for article in root.findall(".//PubmedArticle"):
        pmid = clean_text(article.findtext(".//PMID"))
        title = clean_text("".join(article.findtext(".//ArticleTitle") or ""))
        abstract_parts = []
        for node in article.findall(".//Abstract/AbstractText"):
            label = node.attrib.get("Label", "")
            text = clean_text("".join(node.itertext()))
            if label and text:
                abstract_parts.append(f"{label}: {text}")
            elif text:
                abstract_parts.append(text)
        doi_values = []
        for node in article.findall("./PubmedData/ArticleIdList/ArticleId"):
            if node.attrib.get("IdType", "").lower() == "doi" and node.text:
                doi_values.append(clean_text(node.text))
        records[pmid] = {
            "pmid": pmid,
            "title": title,
            "journal": clean_text(article.findtext(".//Journal/Title")),
            "year": clean_text(article.findtext(".//PubDate/Year")),
            "doi": ";".join(sorted(set(doi_values))),
            "abstract": " ".join(abstract_parts),
        }
    return records


def load_pubmed_query_audit() -> tuple[pd.DataFrame, dict[str, set[str]]]:
    rows = []
    pmid_to_queries: dict[str, set[str]] = defaultdict(set)
    for path in sorted(NCBI_RAW.glob("*_esearch.json")):
        query_id = path.name.removesuffix("_esearch.json")
        spec = PUBMED_QUERY_SPECS.get(query_id, {})
        raw = json.loads(path.read_text(encoding="utf-8"))
        result = raw.get("esearchresult", {})
        ids = [str(v) for v in result.get("idlist", [])]
        for pmid in ids:
            pmid_to_queries[pmid].add(query_id)
        rows.append(
            {
                "round": 14,
                "query_id": query_id,
                "priority_family": spec.get("priority_family", "unmapped"),
                "query_role": spec.get("role", "unmapped"),
                "query": spec.get("query", ""),
                "pubmed_count_available": int(result.get("count", 0) or 0),
                "pubmed_ids_returned": len(ids),
                "returned_pmids": ";".join(ids),
                "query_translation": result.get("querytranslation", ""),
                "raw_file": str(path.relative_to(REPO)),
                "interpretation": "no_hit_is_not_negative" if not ids else "screen_hits_require_exact_pair_verification",
            }
        )
    return pd.DataFrame(rows), pmid_to_queries


def classify_pubmed_record(record: dict[str, str], query_ids: set[str]) -> tuple[str, str]:
    pmid = record.get("pmid", "")
    text = f"{record.get('title', '')} {record.get('abstract', '')}".lower()
    if pmid in KNOWN_EXISTING_EXACT_OR_SUPPORT_PMIDS:
        return "known_existing_gold_or_supporting_context", KNOWN_EXISTING_EXACT_OR_SUPPORT_PMIDS[pmid]
    if any(PUBMED_QUERY_SPECS.get(q, {}).get("role") == "negative_evidence_search" for q in query_ids):
        if "no conversion" in text or "non-producer" in text or "non producer" in text:
            return (
                "negative_evidence_lead_not_label",
                "Potential negative or phenotype-stratification lead; needs exact substrate, organism/strain, assay condition, detection method, and no-product statement.",
            )
        return "context_only_not_negative_label", "Equol/non-producer literature hit without explicit no-conversion edge in title/abstract."
    if any(term in text for term in ["metabol", "conversion", "converted", "hydrolysis", "degradation", "catabol", "ring"]):
        return (
            "candidate_exact_pair_extraction_needed",
            "Biotransformation language present, but title/abstract is not enough for a training edge; extract tables/full text.",
        )
    return "literature_context_only", "Useful search context, not a direct training label."


def build_pubmed_triage(records: dict[str, dict[str, str]], pmid_to_queries: dict[str, set[str]]) -> pd.DataFrame:
    rows = []
    for pmid, query_ids in sorted(pmid_to_queries.items(), key=lambda x: int(x[0])):
        record = records.get(pmid, {"pmid": pmid})
        families = sorted({PUBMED_QUERY_SPECS.get(q, {}).get("priority_family", "unmapped") for q in query_ids})
        status, note = classify_pubmed_record(record, query_ids)
        rows.append(
            {
                "round": 14,
                "pmid": pmid,
                "priority_families": ";".join(families),
                "query_ids": ";".join(sorted(query_ids)),
                "title": record.get("title", ""),
                "journal": record.get("journal", ""),
                "year": record.get("year", ""),
                "doi": record.get("doi", ""),
                "triage_status_round14": status,
                "training_label_allowed_round14": False,
                "negative_label_allowed_round14": False,
                "reason_or_next_check": note,
                "pubmed_url": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/",
            }
        )
    return pd.DataFrame(rows)


def rhea_action_for(query_id: str, equation: str) -> tuple[str, str, bool]:
    lower = equation.lower()
    if query_id in {"urolithin", "equol", "isoxanthohumol", "phenylvalerolactone"}:
        return "rhea_no_direct_support_for_gap", "No direct Rhea reaction for this gut-specific edge in Round14 query.", False
    if query_id == "quercetin" and any(
        marker in lower
        for marker in [
            "glucoside + h2o",
            "quercitrin + h2o",
            "s-adenosyl-l-methionine",
            "sulfate",
            "udp",
        ]
    ):
        return (
            "off_family_quercetin_not_ring_fission",
            "Rhea hit is quercetin glycoside/conjugation chemistry, not gut polyphenol ring-fission; do not count it as solving the ring-fission gap.",
            False,
        )
    direct_good = [
        "chlorogenate + h2o = l-quinate + (e)-caffeate",
        "daidzein 7-o-beta-d-glucoside + h2o = daidzein + beta-d-glucose",
        "taurochenodeoxycholate + h2o = chenodeoxycholate + taurine",
        "glycocholate + h2o = cholate + glycine",
    ]
    if any(pattern in lower for pattern in direct_good):
        return (
            "source_backed_positive_candidate_after_structure_mapping",
            "Approved balanced Rhea reaction; can become a positive candidate only after structure mapping, food-gut relevance, split guard, and generator recall checks.",
            True,
        )
    if "pinoresinol" in lower and "lariciresinol" in lower:
        return (
            "source_rule_template_candidate_direction_review",
            "Rhea supports the redox pair but listed equation direction is opposite of target generator repair; direction/stereochemistry and organism context review required.",
            False,
        )
    if "fucosyllactose" in lower or "fuc-" in lower:
        return (
            "biosynthesis_not_hydrolysis_for_target",
            "Rhea hit is fucosyltransferase/biosynthesis-like, not the target 2-FL hydrolysis edge.",
            False,
        )
    return (
        "context_or_rule_template_only",
        "Reaction may inform rule coverage but is not yet a source-backed food-gut positive.",
        False,
    )


def build_rhea_triage() -> pd.DataFrame:
    rows = []
    for path in sorted(RHEA_RAW.glob("*_rhea.json")):
        query_id = path.name.removesuffix("_rhea.json")
        raw = json.loads(path.read_text(encoding="utf-8"))
        family = RHEA_QUERY_FAMILY.get(query_id, "unmapped")
        results = raw.get("results", [])
        if not results:
            status, note, allowed = rhea_action_for(query_id, "")
            rows.append(
                {
                    "round": 14,
                    "query_id": query_id,
                    "priority_family": family,
                    "rhea_id": "",
                    "equation": "",
                    "status": "",
                    "balanced": "",
                    "transport": "",
                    "rhea_action_round14": status,
                    "positive_candidate_allowed_round14": allowed,
                    "reason_or_next_check": note,
                    "rhea_url": "",
                    "raw_file": str(path.relative_to(REPO)),
                }
            )
            continue
        for item in results:
            status, note, allowed = rhea_action_for(query_id, item.get("equation", ""))
            rows.append(
                {
                    "round": 14,
                    "query_id": query_id,
                    "priority_family": family,
                    "rhea_id": item.get("id", ""),
                    "equation": item.get("equation", ""),
                    "status": item.get("status", ""),
                    "balanced": item.get("balanced", ""),
                    "transport": item.get("transport", ""),
                    "rhea_action_round14": status,
                    "positive_candidate_allowed_round14": allowed,
                    "reason_or_next_check": note,
                    "rhea_url": f"https://www.rhea-db.org/rhea/{item.get('id', '')}",
                    "raw_file": str(path.relative_to(REPO)),
                }
            )
    return pd.DataFrame(rows)


def build_external_standards() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "round": 14,
                "standard_source": "ECREACT / RXN biocatalysis-model",
                "source_url": "https://github.com/rxn4chemistry/biocatalysis-model",
                "evidence": "Data fields include rxn_smiles, ec, and source; source databases are Rhea, BRENDA, PathBank, and MetaNetX; all seven first-level EC classes are represented.",
                "standard_for_foodgut": "A training edge should carry reaction SMILES or mapped substrate/product structures, EC/database source, source table/version, and provenance.",
            },
            {
                "round": 14,
                "standard_source": "RetroPathRL / RetroRules",
                "source_url": "https://github.com/brsynth/RetroPathRL",
                "evidence": "Rule input requires Rule_ID, Reaction_ID, Rule_SMARTS; optional fields include substrates, products, EC, and rule usage.",
                "standard_for_foodgut": "A SMARTS rule is a generator object, not a biological label. Promote rules only with source-reaction provenance and direction/license review.",
            },
            {
                "round": 14,
                "standard_source": "gapseq",
                "source_url": "https://github.com/jotech/gapseq",
                "evidence": "Combines pathway prediction, reaction databases, transporter inference, model construction, and gap filling; README lists source-license constraints.",
                "standard_for_foodgut": "Genome/model/pathway absence can support conditional negatives only under versioned organism/scope, never global negatives.",
            },
            {
                "round": 14,
                "standard_source": "EnzymeMap / enzymatic reaction curation literature",
                "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC10718068/",
                "evidence": "Large enzymatic reaction datasets expand beyond ECREACT but still rely on curated/validated reaction sources.",
                "standard_for_foodgut": "Dataset scale helps coverage, but validation still depends on exact reaction identity and source traceability.",
            },
        ]
    )


def build_blocker_diagnosis(
    gaps: pd.DataFrame,
    pubmed_query: pd.DataFrame,
    rhea: pd.DataFrame,
) -> pd.DataFrame:
    rows = []
    generator = safe_read_csv(ROUND11_GENERATOR)
    if not generator.empty and "priority_family" not in generator.columns:
        generator = generator.copy()
        generator["priority_family"] = generator["pair_key"].map(KNOWN_MISS_PAIR_FAMILY).fillna("")
    for item in gaps.itertuples(index=False):
        family = item.priority_family
        pubmed_hits = int(pubmed_query.loc[pubmed_query["priority_family"].eq(family), "pubmed_ids_returned"].sum())
        rhea_hits = int((rhea["priority_family"].eq(family) & rhea["rhea_id"].astype(str).ne("")).sum())
        rhea_family_positive = rhea[
            rhea["priority_family"].eq(family)
            & rhea["positive_candidate_allowed_round14"].astype(bool)
        ].drop_duplicates(["rhea_id", "equation"])
        rhea_positive = int(len(rhea_family_positive))
        known_miss = ""
        if not generator.empty:
            known_miss = ";".join(generator.loc[generator["priority_family"].eq(family), "pair_key"].astype(str).tolist())
        if known_miss:
            blocker_type = "generator_recall_blocker_plus_family_coverage"
            judgment = "Existing gold edge is missing from clean candidates; fix generator before retraining or duplicate imports."
        elif rhea_positive:
            blocker_type = "source_backed_sample_coverage_gap"
            judgment = "Rhea provides source-backed reactions that can expand positives after mapping and generator checks."
        elif rhea_hits or pubmed_hits:
            blocker_type = "evidence_extraction_gap"
            judgment = "There is external evidence, but it is not yet exact/source-mapped enough for labels."
        else:
            blocker_type = "source_sparsity_gap"
            judgment = "Public exact-reaction evidence is sparse; do not invent reactions or treat no-hit as negative."
        rows.append(
            {
                "round": 14,
                "priority_rank": item.priority_rank,
                "module": item.module,
                "priority_family": family,
                "local_fullrule_recall": item.fullrule_recall,
                "challenge_case_count": item.challenge_case_count,
                "known_gold_missing_clean_candidate_keys": known_miss or item.round9_known_gold_missing_clean_candidate_keys,
                "pubmed_ids_returned_round14": pubmed_hits,
                "rhea_reactions_returned_round14": rhea_hits,
                "rhea_positive_candidates_round14": rhea_positive,
                "blocker_type_round14": blocker_type,
                "judgment_round14": judgment,
                "production_action": "fix_generator_then_curate_labels" if known_miss else "collect_verify_map_before_training",
            }
        )
    return pd.DataFrame(rows)


def build_positive_candidates(rhea: pd.DataFrame, pubmed: pd.DataFrame) -> pd.DataFrame:
    rows = []
    allowed = rhea[rhea["positive_candidate_allowed_round14"].astype(bool)].copy()
    allowed = allowed.drop_duplicates(["priority_family", "rhea_id", "equation"])
    for item in allowed.itertuples(index=False):
        rows.append(
            {
                "round": 14,
                "candidate_type": "rhea_source_backed_reaction",
                "priority_family": item.priority_family,
                "source_id": f"RHEA:{item.rhea_id}",
                "source_url": item.rhea_url,
                "reaction_or_pair": item.equation,
                "evidence_strength": "database_approved_balanced_reaction",
                "training_allowed_now": False,
                "why_not_yet": "Needs structure mapping to local SMILES/InChIKey, food-gut relevance annotation, split/holdout guard, generator recall check, and duplicate check.",
                "next_round_action": "map structures and test whether clean generator can enumerate the product",
            }
        )
    known = pubmed[pubmed["triage_status_round14"].eq("known_existing_gold_or_supporting_context")]
    for item in known.itertuples(index=False):
        rows.append(
            {
                "round": 14,
                "candidate_type": "pubmed_known_existing_gold_or_support",
                "priority_family": item.priority_families,
                "source_id": f"PMID:{item.pmid}",
                "source_url": item.pubmed_url,
                "reaction_or_pair": item.reason_or_next_check,
                "evidence_strength": "already_known_or_supporting_context",
                "training_allowed_now": False,
                "why_not_yet": "Do not duplicate-import existing gold/supporting rows; generator recall repair is the correct action.",
                "next_round_action": "use source only for rule/source recovery and generator tests",
            }
        )
    if not rows:
        rows.append(
            {
                "round": 14,
                "candidate_type": "none",
                "priority_family": "",
                "source_id": "",
                "source_url": "",
                "reaction_or_pair": "",
                "evidence_strength": "",
                "training_allowed_now": False,
                "why_not_yet": "No Round14 candidate passed even source-backed candidate gates.",
                "next_round_action": "continue exact extraction before training",
            }
        )
    return pd.DataFrame(rows)


def build_negative_strategy(pubmed: pd.DataFrame) -> pd.DataFrame:
    negative_leads = pubmed[
        pubmed["triage_status_round14"].isin(
            ["negative_evidence_lead_not_label", "context_only_not_negative_label"]
        )
    ]
    return pd.DataFrame(
        [
            {
                "round": 14,
                "negative_route": "assay_negative",
                "candidate_sources_round14": ";".join(negative_leads["pmid"].astype(str).head(12).tolist()),
                "allowed_for_training_now": False,
                "acceptance_gate": "Paper states exact substrate, exact organism/strain/community, assay condition/time, detection method, and explicit no-product/no-conversion outcome.",
                "reject_if": "The paper only says non-producer phenotype, no PubMed/Rhea hit, or lower abundance without exact assay edge.",
                "use_in_model": "biological negative label only after manual exact extraction",
            },
            {
                "round": 14,
                "negative_route": "conditional_genome_or_model_absence",
                "candidate_sources_round14": "gapseq/ModelSEED/MGnify/strain genome models; not fetched as labels in Round14",
                "allowed_for_training_now": False,
                "acceptance_gate": "Versioned genome/model lacks required enzyme/pathway under defined strain and substrate scope; recorded as conditional, not global.",
                "reject_if": "Missing annotation is used as universal impossibility or without database version/scope.",
                "use_in_model": "context feature or conditional negative with scope columns",
            },
            {
                "round": 14,
                "negative_route": "hard_decoy_for_ranking",
                "candidate_sources_round14": "clean generator outputs and source-backed rule overlays",
                "allowed_for_training_now": True,
                "acceptance_gate": "Generated by allowed rules, not known positive, structurally plausible, split-safe, and labeled as ranking decoy not biological false.",
                "reject_if": "Decoy overlaps known/source-backed positive or future validation queue.",
                "use_in_model": "ranking contrastive negative only",
            },
            {
                "round": 14,
                "negative_route": "pubmed_or_rhea_no_hit",
                "candidate_sources_round14": "zero-hit strict queries and sparse Rhea queries",
                "allowed_for_training_now": False,
                "acceptance_gate": "None: no-hit is only a search gap.",
                "reject_if": "Always reject as biological negative.",
                "use_in_model": "unlabeled or priority for source search",
            },
        ]
    )


def build_pipeline() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "stage": 1,
                "stage_name": "collect",
                "input": "priority family + local generator miss + database/literature query",
                "output": "raw PubMed/Rhea/GitHub evidence files and query audit rows",
                "pass_gate": "query, date, raw file, source URL saved",
            },
            {
                "stage": 2,
                "stage_name": "verify_truth",
                "input": "candidate paper/database row",
                "output": "exact substrate-product, direction, organism/enzyme/context, source id",
                "pass_gate": "exact edge is visible in source; broad review/context rejected",
            },
            {
                "stage": 3,
                "stage_name": "map_structures",
                "input": "accepted exact edge",
                "output": "SMILES, InChIKey, names, synonyms, database ids",
                "pass_gate": "substrate/product structures parse and match local chemistry",
            },
            {
                "stage": 4,
                "stage_name": "modify_candidate_source_or_rule",
                "input": "source-backed pair or source-backed rule",
                "output": "import manifest or generator overlay manifest",
                "pass_gate": "license, direction, stereochemistry, duplicate, and split guards pass",
            },
            {
                "stage": 5,
                "stage_name": "generator_recall_check",
                "input": "candidate pair and deployment-faithful generator",
                "output": "true product in clean candidates or blocked reason",
                "pass_gate": "true product is enumerable before LTR retraining",
            },
            {
                "stage": 6,
                "stage_name": "freeze_and_review",
                "input": "CSV manifests + generated candidate outputs",
                "output": "commit hash, diff, model-quality review",
                "pass_gate": "AI/human review finds no leakage, duplicate import, or unsupported label",
            },
            {
                "stage": 7,
                "stage_name": "train_and_evaluate",
                "input": "gated positives + scoped negatives + hard decoys",
                "output": "split-safe metrics and failure panel",
                "pass_gate": "metrics reported with n, CI, family coverage, and challenge-panel failures",
            },
        ]
    )


def write_doc(
    standards: pd.DataFrame,
    pubmed_query: pd.DataFrame,
    pubmed: pd.DataFrame,
    rhea: pd.DataFrame,
    blockers: pd.DataFrame,
    positive: pd.DataFrame,
    negative: pd.DataFrame,
) -> None:
    metrics = safe_read_csv(METRICS)
    npts = ""
    if not metrics.empty and "n_pts" in metrics.columns:
        npts = ", ".join(
            f"{row.module}/{row.method}:n={row.n_pts}"
            for row in metrics[["module", "method", "n_pts"]].drop_duplicates().itertuples(index=False)
            if row.method == "LTR_full"
        )

    rhea_source_candidates = int(
        positive[
            positive["candidate_type"].eq("rhea_source_backed_reaction")
        ].drop_duplicates(["source_id", "reaction_or_pair"]).shape[0]
    )
    known_context = int(pubmed["triage_status_round14"].eq("known_existing_gold_or_supporting_context").sum())
    negative_leads = int(pubmed["triage_status_round14"].isin(["negative_evidence_lead_not_label", "context_only_not_negative_label"]).sum())

    blocker_counts = blockers["blocker_type_round14"].value_counts().to_dict()
    top_blockers = "\n".join(
        f"- `{row.priority_family}`: `{row.blocker_type_round14}`; {row.judgment_round14}"
        for row in blockers.itertuples(index=False)
    )
    positive_preview = "\n".join(
        f"- `{row.source_id}`: {row.reaction_or_pair}"
        for row in positive[positive["candidate_type"].ne("none")].head(12).itertuples(index=False)
    ) or "- No positive candidate passed source-backed candidate gates."
    negative_preview = "\n".join(
        f"- `{row.negative_route}`: allowed_now={row.allowed_for_training_now}; {row.acceptance_gate}"
        for row in negative.itertuples(index=False)
    )

    text = f"""# Round14 Reaction Sample Expansion Review

Date: 2026-06-04

## Direct Answer

Ray's diagnosis is mostly right, but the production blocker has two layers:

1. **Reaction-family coverage is incomplete.** Rhea supports some missing families, but gut-specific reactions such as urolithin, equol, isoxanthohumol, and flavanol ring-fission remain sparse in curated reaction databases.
2. **Generator recall is still the sharper blocker.** Several known gold positives already exist, but clean candidate generation does not retain the true product. Retraining cannot fix a product that never reaches the ranker.

The current high scores are not production evidence because `clean2_metrics.csv` evaluates only small panels (`{npts}`).

## What Was Run

- PubMed Round14 search: `{len(pubmed_query)}` query rows, `{pubmed_query['pubmed_ids_returned'].sum() if not pubmed_query.empty else 0}` returned PMID slots, `{pubmed['pmid'].nunique() if not pubmed.empty else 0}` unique PMID records parsed.
- Rhea Round14 search: `{rhea['query_id'].nunique() if not rhea.empty else 0}` query terms, `{rhea['rhea_id'].astype(str).ne('').sum() if not rhea.empty else 0}` returned reaction rows.
- GitHub standards checked: ECREACT/RXN biocatalysis, RetroPathRL/RetroRules, gapseq, and EnzymeMap literature.

## Main Production Diagnosis

Blocker counts: `{blocker_counts}`

{top_blockers}

## Safe Positive Expansion Candidates

Round14 found `{rhea_source_candidates}` Rhea source-backed reaction candidates and `{known_context}` PubMed rows that are already known/supporting context. **None are allowed into training immediately.**

{positive_preview}

Why not train yet:

- structures must be mapped to local SMILES/InChIKey;
- duplicate existing positives must be excluded;
- food-gut relevance and source license must be recorded;
- the deployment-faithful generator must enumerate the true product;
- split leakage must be checked.

## Negative Sample Logic

Round14 found `{negative_leads}` PubMed hits in equol/non-producer/no-conversion searches, but these are not negative labels until exact assay conditions are extracted.

{negative_preview}

Hard decoys are the only Round14 route that can be used now, and only as ranking decoys, not biological negatives.

## Files Created

- `data/curation/round14_external_model_standards.csv`
- `data/curation/round14_pubmed_query_audit.csv`
- `data/curation/round14_pubmed_candidate_triage.csv`
- `data/curation/round14_rhea_reaction_triage.csv`
- `data/curation/round14_reaction_type_blocker_diagnosis.csv`
- `data/curation/round14_positive_sample_candidates.csv`
- `data/curation/round14_negative_sample_strategy.csv`
- `data/curation/round14_collect_verify_modify_pipeline.csv`

Raw evidence:

- `data/curation/ncbi_round14_raw/*.json`
- `data/curation/rhea_round14_raw/*.json`
- `data/curation/ncbi_round14_raw/*.xml` is intentionally ignored by git.

## Round15 Direction

1. Map structures for Rhea positive candidates, starting with chlorogenate hydrolysis, daidzein glycoside hydrolysis, and bile acid deconjugation.
2. For known gold generator misses, continue generator overlay repair before importing duplicate labels.
3. For negative samples, extract explicit no-conversion assay rows only; otherwise keep rows as unlabeled or hard decoys.
4. After mapping, run generator recall before any retraining.

## Claim Boundary

Allowed after Round14:

> We have a source-backed queue of candidate positive reactions, a strict negative-sample policy, and a diagnosis separating reaction coverage gaps from generator recall failures.

Forbidden:

> Round14 has expanded the training set or proved production readiness.
"""
    OUT_DOC.write_text(textwrap.dedent(text).strip() + "\n", encoding="utf-8")


def main() -> None:
    CURATION.mkdir(parents=True, exist_ok=True)
    DOCS.mkdir(parents=True, exist_ok=True)

    standards = build_external_standards()
    pubmed_query, pmid_to_queries = load_pubmed_query_audit()
    pubmed_records = parse_pubmed_xml(NCBI_RAW / "round14_pubmed_efetch.xml")
    pubmed = build_pubmed_triage(pubmed_records, pmid_to_queries)
    rhea = build_rhea_triage()
    gaps = safe_read_csv(ROUND10_GAPS)
    blockers = build_blocker_diagnosis(gaps, pubmed_query, rhea)
    positive = build_positive_candidates(rhea, pubmed)
    negative = build_negative_strategy(pubmed)
    pipeline = build_pipeline()

    outputs = [
        (OUT_STANDARDS, standards),
        (OUT_PUBMED_QUERY, pubmed_query),
        (OUT_PUBMED_TRIAGE, pubmed),
        (OUT_RHEA, rhea),
        (OUT_BLOCKERS, blockers),
        (OUT_POSITIVE, positive),
        (OUT_NEGATIVE, negative),
        (OUT_PIPELINE, pipeline),
    ]
    for path, table in outputs:
        table.to_csv(path, index=False)
        print(f"wrote {path} rows={len(table)}")

    write_doc(standards, pubmed_query, pubmed, rhea, blockers, positive, negative)
    print(f"wrote {OUT_DOC}")
    print("round14_blockers=", blockers["blocker_type_round14"].value_counts().to_dict())
    print(
        "round14_rhea_positive_candidates=",
        int(
            positive[
                positive["candidate_type"].eq("rhea_source_backed_reaction")
            ].drop_duplicates(["source_id", "reaction_or_pair"]).shape[0]
        ),
    )
    print("round14_pubmed_known_or_support=", int(pubmed["triage_status_round14"].eq("known_existing_gold_or_supporting_context").sum()))


if __name__ == "__main__":
    main()
