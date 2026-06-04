#!/usr/bin/env python
"""Build Round10 evidence-acquisition audit tables.

Round10 turns Ray's production concern into a reproducible curation pipeline:
local generation gaps -> external evidence screens -> gated positive/negative
sample policy -> next acquisition queue.
"""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd


REPO = Path(__file__).resolve().parents[1]
CURATION = REPO / "data" / "curation"
NCBI_RAW = CURATION / "ncbi_round10_raw"
RHEA_RAW = CURATION / "rhea_round10_raw"

FULLRULE = REPO / "outputs" / "gen_gap" / "fullrule_hit_miss.parquet"
CHALLENGE = REPO / "data" / "production_challenge_panel.csv"
ROUND9_IMPORT = CURATION / "positive_sample_import_round9.csv"

OUT_GAPS = CURATION / "reaction_family_gap_priorities_round10.csv"
OUT_STANDARDS = CURATION / "external_repository_reaction_standards_round10.csv"
OUT_PUBMED = CURATION / "pubmed_candidate_screen_round10.csv"
OUT_RHEA = CURATION / "rhea_reaction_screen_round10.csv"
OUT_POLICY = CURATION / "sample_acquisition_policy_round10.csv"
OUT_PIPELINE = CURATION / "evidence_expansion_pipeline_round10.csv"
OUT_QUEUE = CURATION / "round10_next_evidence_queue.csv"
OUT_LOG = CURATION / "round10_iteration_log.csv"


PUBMED_QUERIES = {
    "caffeic_dihydrocaffeic": "gut microbiota caffeic acid dihydrocaffeic acid",
    "daidzein_equol": "gut microbiota daidzein equol",
    "catechin_valerolactone": "gut microbiota catechin valerolactone",
    "quercetin_ring_fission": "gut microbiota quercetin dihydroxyphenylacetic acid",
    "urolithin_c_a": "urolithin C urolithin A Enterocloster",
    "fucosyllactose_fucose": "2'-fucosyllactose fucose Bifidobacterium infantis",
    "pinoresinol_lariciresinol": "pinoresinol lariciresinol human intestinal microflora",
    "isoxanthohumol_8pn": "isoxanthohumol 8-prenylnaringenin human intestine",
    "assay_negative_no_conversion": "gut microbiota no conversion substrate metabolite assay negative",
}

RHEA_QUERIES = {
    "caffeic_acid": "caffeic acid",
    "dihydrocaffeic_acid": "dihydrocaffeic acid",
    "daidzein": "daidzein",
    "equol": "equol",
    "catechin": "catechin",
    "urolithin": "urolithin",
    "fucosyllactose": "fucosyllactose",
    "bile_salt_hydrolase": "bile salt hydrolase",
    "phenolic_acid_decarboxylase": "phenolic acid decarboxylase",
}


def read_json(path: Path) -> dict:
    if not path.exists():
        return {}
    return json.loads(path.read_text(encoding="utf-8"))


def compact_join(values: list[str], limit: int = 4) -> str:
    clean = [str(v).strip() for v in values if str(v).strip()]
    suffix = "" if len(clean) <= limit else f"; +{len(clean) - limit} more"
    return "; ".join(clean[:limit]) + suffix


def challenge_cases(challenge: pd.DataFrame, contains_terms: list[str]) -> tuple[int, str]:
    if challenge.empty:
        return 0, ""
    family = challenge["reaction_family"].fillna("").str.lower()
    mask = pd.Series(False, index=challenge.index)
    for term in contains_terms:
        mask = mask | family.str.contains(term.lower(), regex=False)
    hits = challenge[mask]
    return len(hits), compact_join(hits["case_id"].astype(str).tolist(), limit=8)


def local_gap_stats(fullrule: pd.DataFrame, module: str, category_terms: list[str]) -> dict[str, object]:
    if fullrule.empty:
        return {"fullrule_rows": 0, "fullrule_hits": 0, "fullrule_recall": ""}
    category = fullrule["reaction_category"].fillna("").str.lower()
    mask = fullrule["module"].eq(module)
    term_mask = pd.Series(False, index=fullrule.index)
    for term in category_terms:
        term_mask = term_mask | category.str.contains(term.lower(), regex=False)
    hits = fullrule[mask & term_mask]
    if hits.empty:
        return {"fullrule_rows": 0, "fullrule_hits": 0, "fullrule_recall": ""}
    return {
        "fullrule_rows": len(hits),
        "fullrule_hits": int(hits["hit"].sum()),
        "fullrule_recall": round(float(hits["hit"].mean()), 4),
    }


def round9_missing_keys(round9: pd.DataFrame, terms: list[str]) -> str:
    if round9.empty or "clean_candidates_full_hit_count" not in round9.columns or not terms:
        return ""
    text = (
        round9.get("substrate_name", pd.Series("", index=round9.index)).fillna("").astype(str)
        + " "
        + round9.get("product_name", pd.Series("", index=round9.index)).fillna("").astype(str)
    ).str.lower()
    term_mask = pd.Series(False, index=round9.index)
    for term in terms:
        term_mask = term_mask | text.str.contains(term.lower(), regex=False)
    mask = round9["clean_candidates_full_hit_count"].eq(0) & term_mask
    return compact_join(round9.loc[mask, "pair_key"].astype(str).tolist(), limit=8)


def build_gap_priorities() -> pd.DataFrame:
    fullrule = pd.read_parquet(FULLRULE) if FULLRULE.exists() else pd.DataFrame()
    challenge = pd.read_csv(CHALLENGE) if CHALLENGE.exists() else pd.DataFrame()
    round9 = pd.read_csv(ROUND9_IMPORT) if ROUND9_IMPORT.exists() else pd.DataFrame()

    specs = [
        {
            "priority_rank": 1,
            "module": "A",
            "priority_family": "polyphenol_ring_fission",
            "local_terms": ["bond_cleavage_ring_opening", "polyphenol_quercetin_fission"],
            "challenge_terms": ["ring_fission"],
            "round9_terms": [],
            "why": "Challenge families include catechin, epicatechin, naringenin, and quercetin ring-fission; local full-rule recall is low for ring-opening categories.",
            "external_database_route": "PubMed exact-pair extraction; Rhea is sparse for gut flavanol/flavonol ring fission.",
        },
        {
            "priority_rank": 2,
            "module": "A",
            "priority_family": "urolithin_dehydroxylation",
            "local_terms": ["polyphenol_urolithin", "functional_group_removal"],
            "challenge_terms": ["urolithin", "ellagitannin"],
            "round9_terms": ["urolithin"],
            "why": "Round9 shows urolithin C -> urolithin A is an existing gold positive but absent from clean candidate generation.",
            "external_database_route": "PubMed exact-pair and enzyme/organism extraction; Rhea query returned no urolithin reactions.",
        },
        {
            "priority_rank": 3,
            "module": "A",
            "priority_family": "hydroxycinnamate_reduction_and_hydrolysis",
            "local_terms": ["hydroxycinnamic_acid", "reduction_hydrogenation", "hydrolysis_deconjugation"],
            "challenge_terms": ["hydroxycinnamate", "phenolic_ester"],
            "round9_terms": [],
            "why": "Food phenolic acids are core gut metabolites; direct reduction evidence needs exact-pair validation, while chlorogenate hydrolysis has stronger database support.",
            "external_database_route": "Rhea for chlorogenate/caffeate reactions plus PubMed exact-pair extraction for gut anaerobic reductions.",
        },
        {
            "priority_rank": 4,
            "module": "A",
            "priority_family": "isoflavone_reductive_and_glycoside_metabolism",
            "local_terms": ["reduction_hydrogenation", "group_transfer_reconjugation_modification"],
            "challenge_terms": ["isoflavone"],
            "round9_terms": [],
            "why": "Daidzein/equol and puerarin/daidzein are high-value production challenges; Rhea covers some glycoside chemistry but not equol production.",
            "external_database_route": "Rhea for daidzein glycosides and PubMed for equol-producing strains and stepwise exact pairs.",
        },
        {
            "priority_rank": 5,
            "module": "A",
            "priority_family": "lignan_redox_and_deglycosylation",
            "local_terms": ["group_transfer_reconjugation_modification", "oxidation_dehydrogenation"],
            "challenge_terms": ["lignan"],
            "round9_terms": ["pinoresinol", "lariciresinol"],
            "why": "Round9 shows pinoresinol -> lariciresinol is an existing gold positive but absent from clean candidate generation.",
            "external_database_route": "PubMed exact-pair extraction; database coverage is likely weaker than literature for gut lignan transformations.",
        },
        {
            "priority_rank": 6,
            "module": "A",
            "priority_family": "prenylflavonoid_o_demethylation",
            "local_terms": ["functional_group_removal", "group_transfer_reconjugation_modification"],
            "challenge_terms": [],
            "round9_terms": ["isoxanthohumol", "8-prenylnaringenin"],
            "why": "Round9 shows isoxanthohumol -> 8-prenylnaringenin is an existing gold positive but absent from clean candidate generation.",
            "external_database_route": "PubMed exact-pair extraction and source-rule recovery; do not infer this transformation from generic demethylation alone.",
        },
        {
            "priority_rank": 7,
            "module": "B",
            "priority_family": "glycoside_and_hmo_hydrolysis",
            "local_terms": ["group_transfer_reconjugation_modification", "starch and sucrose", "aminosugar"],
            "challenge_terms": ["deglucosylation", "glycoside"],
            "round9_terms": ["fucosyllactose", "fucose"],
            "why": "B module shows deceptively high gold recall on n=10, but Round9 shows 2'-FL -> L-fucose is existing gold and absent from clean candidates.",
            "external_database_route": "PubMed for strain/enzyme HMO utilization; Rhea fucosyllactose hits are mostly biosynthetic fucosyltransferase reactions, not the target hydrolysis edge.",
        },
        {
            "priority_rank": 8,
            "module": "D",
            "priority_family": "bile_acid_deconjugation_and_lipid_context",
            "local_terms": ["bile", "glycerophospholipid", "fatty acid"],
            "challenge_terms": ["bile"],
            "round9_terms": [],
            "why": "D has better full-rule recall, but bile/lipid rules need organism and enzyme context before production claims.",
            "external_database_route": "Rhea and curated bacterial genome/pathway sources; use as context-aware positives or conditional negatives.",
        },
        {
            "priority_rank": 9,
            "module": "C",
            "priority_family": "amino_acid_catabolism_and_redox",
            "local_terms": ["methionine", "cysteine", "valine", "leucine", "isoleucine", "oxidation_dehydrogenation", "deamination"],
            "challenge_terms": ["amino"],
            "round9_terms": [],
            "why": "C module has several low-recall amino-acid categories; important for gut microbial metabolism but less represented in the current challenge panel.",
            "external_database_route": "Rhea/MetaNetX/ECREACT plus organism-context databases; avoid importing broad pathway reactions without food-gut relevance.",
        },
    ]

    rows = []
    for spec in specs:
        local = local_gap_stats(fullrule, spec["module"], spec["local_terms"])
        count, cases = challenge_cases(challenge, spec["challenge_terms"])
        rows.append(
            {
                "round": 10,
                "priority_rank": spec["priority_rank"],
                "module": spec["module"],
                "priority_family": spec["priority_family"],
                **local,
                "challenge_case_count": count,
                "challenge_cases": cases,
                "round9_known_gold_missing_clean_candidate_keys": round9_missing_keys(
                    round9, spec.get("round9_terms", [])
                ),
                "production_blocker": spec["why"],
                "external_database_route": spec["external_database_route"],
                "positive_sample_route": "exact substrate-product evidence only; source PMID/DOI/Rhea/MetaNetX/EC required",
                "negative_sample_route": "assay no-conversion or context-specific organism/genome absence; otherwise use hard decoy or unlabeled",
                "next_action_round10": "screen external hits, extract exact pairs, then fix generator before retraining",
            }
        )
    return pd.DataFrame(rows)


def build_external_standards() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "round": 10,
                "source_name": "ECREACT / RXN biocatalysis model",
                "source_type": "GitHub + Nature Communications paper",
                "url": "https://github.com/rxn4chemistry/biocatalysis-model/blob/main/README.md",
                "coverage_statement": "Enzymatic reactions from Rhea, BRENDA, PathBank, and MetaNetX; all 7 first-level EC classes; fields include rxn_smiles, ec, source.",
                "standard_for_this_project": "Every training reaction should carry source and EC/database provenance when possible.",
            },
            {
                "round": 10,
                "source_name": "RetroPathRL / RetroRules",
                "source_type": "GitHub rule-based biosynthesis planning",
                "url": "https://github.com/brsynth/RetroPathRL/blob/master/README.md",
                "coverage_statement": "Uses mono-component reaction rules from RetroRules for candidate-generation and pathway search.",
                "standard_for_this_project": "Rule hits are candidates, not labels; import rules only with source-reaction provenance.",
            },
            {
                "round": 10,
                "source_name": "gapseq",
                "source_type": "GitHub + Genome Biology paper",
                "url": "https://github.com/jotech/gapseq/blob/master/README.md",
                "coverage_statement": "Combines pathway prediction, transporter inference, metabolic model construction, and gap filling from curated reaction/pathway/genome resources.",
                "standard_for_this_project": "Use organism/genome/pathway context for conditional support and conditional negatives, not global reaction truth.",
            },
            {
                "round": 10,
                "source_name": "gutSMASH",
                "source_type": "GitHub + gut microbiome specialized metabolism",
                "url": "https://github.com/victoriapascal/gutsmash/blob/gutsmash/README.md",
                "coverage_statement": "Predicts known and novel anaerobic gut metabolic gene clusters validated against a curated dataset.",
                "standard_for_this_project": "Gene-cluster evidence can support biological plausibility, but exact substrate-product labels still need database or paper evidence.",
            },
        ]
    )


def build_pubmed_screen() -> pd.DataFrame:
    summary = read_json(NCBI_RAW / "pubmed_esummary_round10.json").get("result", {})
    rows = []
    for slug, term in PUBMED_QUERIES.items():
        search = read_json(NCBI_RAW / f"{slug}_esearch.json")
        ids = search.get("esearchresult", {}).get("idlist", [])
        count = search.get("esearchresult", {}).get("count", "")
        for rank, pmid in enumerate(ids, start=1):
            item = summary.get(str(pmid), {})
            title = item.get("title", "")
            doi = ""
            for article_id in item.get("articleids", []) or []:
                if article_id.get("idtype") == "doi":
                    doi = article_id.get("value", "")
                    break
            rows.append(
                {
                    "round": 10,
                    "query_slug": slug,
                    "query": term,
                    "pubmed_available_count": count,
                    "rank": rank,
                    "pmid": pmid,
                    "title": title,
                    "journal": item.get("source", ""),
                    "pubdate": item.get("pubdate", ""),
                    "doi": doi,
                    "pubmed_url": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/",
                    "screen_status_round10": classify_pubmed(slug, title),
                    "training_decision_round10": "not_training_label_until_exact_pair_extracted",
                    "required_next_check": "read abstract/full text or table for exact substrate, product, direction, organism/enzyme, assay context, and structure IDs",
                }
            )
        if not ids:
            rows.append(
                {
                    "round": 10,
                    "query_slug": slug,
                    "query": term,
                    "pubmed_available_count": count,
                    "rank": "",
                    "pmid": "",
                    "title": "",
                    "journal": "",
                    "pubdate": "",
                    "doi": "",
                    "pubmed_url": "",
                    "screen_status_round10": "no_pubmed_hits_for_query",
                    "training_decision_round10": "no_label",
                    "required_next_check": "revise query or search source databases; do not convert absence into negative",
                }
            )
    return pd.DataFrame(rows)


def classify_pubmed(slug: str, title: str) -> str:
    text = f"{slug} {title}".lower()
    if "review" in text or "overview" in text:
        return "review_or_context_only"
    if "metabolite" in text or "metabolism" in text or "biotransformation" in text:
        return "candidate_literature_screen"
    if "no conversion" in text or "negative" in text:
        return "possible_assay_negative_screen"
    return "candidate_literature_screen"


def build_rhea_screen() -> pd.DataFrame:
    rows = []
    for slug, term in RHEA_QUERIES.items():
        raw = read_json(RHEA_RAW / f"{slug}.json")
        results = raw.get("results", [])
        if not results:
            rows.append(
                {
                    "round": 10,
                    "query_slug": slug,
                    "query": term,
                    "rhea_available_count": raw.get("count", 0),
                    "rhea_id": "",
                    "equation": "",
                    "balanced": "",
                    "transport": "",
                    "rhea_url": "",
                    "screen_status_round10": "no_rhea_hit_for_query",
                    "training_decision_round10": "no_label; database absence is not a negative",
                    "required_next_check": "try exact compound synonyms, MetaNetX, BRENDA, PubMed, or organism-specific literature",
                }
            )
            continue
        for item in results:
            rid = str(item.get("id", ""))
            equation = str(item.get("equation", ""))
            rows.append(
                {
                    "round": 10,
                    "query_slug": slug,
                    "query": term,
                    "rhea_available_count": raw.get("count", ""),
                    "rhea_id": rid,
                    "equation": equation,
                    "balanced": item.get("balanced", ""),
                    "transport": item.get("transport", ""),
                    "rhea_url": f"https://www.rhea-db.org/rhea/{rid}" if rid else "",
                    "screen_status_round10": classify_rhea(slug, equation),
                    "training_decision_round10": "source_rule_candidate_not_positive_label"
                    if rid
                    else "no_label",
                    "required_next_check": "map participants to substrate/product IDs, confirm direction, EC/source, license, and whether this is exact to food-gut target",
                }
            )
    return pd.DataFrame(rows)


def classify_rhea(slug: str, equation: str) -> str:
    text = f"{slug} {equation}".lower()
    if "urolithin" in slug or "equol" in slug:
        return "gap_confirmed_in_rhea_screen"
    if "fucosyllactose" in slug and "h2o" not in text:
        return "related_biosynthesis_not_target_hydrolysis"
    if "catechin" in slug and "valerolactone" not in text:
        return "related_plant_or_redox_reaction_not_gut_ring_fission"
    if "bile_salt" in slug and ("glycocholate" in text or "tauro" in text):
        return "strong_source_rule_candidate_for_bile_deconjugation"
    if "phenolic_acid" in slug or "dihydrocaffeic" in slug or "caffeic" in slug:
        return "source_rule_candidate_for_phenolic_acid_family"
    if "daidzein" in slug:
        return "source_rule_candidate_for_isoflavone_glycoside_or_related_step"
    return "source_rule_candidate_requires_exactness_check"


def build_policy() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "round": 10,
                "sample_kind": "exact_positive",
                "can_train_as": "positive",
                "allowed_when": "Exact substrate-product edge appears in paper/table/database with stable chemical IDs and direction/context.",
                "forbidden_when": "Only family, pathway, mechanism, review, or rule-hit evidence is available.",
                "required_fields": "substrate/product name; SMILES; InChIKey; source PMID/DOI/database/version; organism/enzyme; direction; evidence tier; split guard",
                "negative_logic": "not applicable",
                "current_use": "gold/silver positive after gates",
            },
            {
                "round": 10,
                "sample_kind": "source_traceable_rule",
                "can_train_as": "candidate_generator_rule_or_feature",
                "allowed_when": "Rule maps to source reaction ID, EC/Rhea/MetaNetX/BRENDA/RetroRules provenance, and license/version.",
                "forbidden_when": "Anonymous SMARTS, MicrobeRX-only rule without recoverable source, or invented transformation.",
                "required_fields": "rule_id; SMARTS; source_reaction_id; EC; source database; direction; participant mapping; license",
                "negative_logic": "not a negative source",
                "current_use": "generator/ranking feature, not label by itself",
            },
            {
                "round": 10,
                "sample_kind": "assay_negative",
                "can_train_as": "conditional_negative_or_high_confidence_negative_under_context",
                "allowed_when": "Paper explicitly reports no conversion or no product under defined organism/enzyme/assay/detection conditions.",
                "forbidden_when": "No database hit, no rule hit, or no mention in literature.",
                "required_fields": "substrate; missing product/product class; organism/enzyme; assay time; detection method; detection limit; PMID/DOI",
                "negative_logic": "valid only under stated assay context",
                "current_use": "rare high-value negative; keep separate from global negatives",
            },
            {
                "round": 10,
                "sample_kind": "conditional_negative",
                "can_train_as": "context_negative",
                "allowed_when": "Versioned organism/genome/pathway database lacks gene/pathway under a defined condition and scope.",
                "forbidden_when": "Used as global biological impossibility.",
                "required_fields": "organism/strain; database version; pathway/gene/reaction checked; scope; source link",
                "negative_logic": "context-specific absence only",
                "current_use": "calibration or context-aware ranking, not universal negative",
            },
            {
                "round": 10,
                "sample_kind": "hard_decoy",
                "can_train_as": "ranking_negative",
                "allowed_when": "Same substrate candidate generated by valid rules but absent from curated positives after leakage-safe construction.",
                "forbidden_when": "Claimed as biologically impossible.",
                "required_fields": "substrate; candidate product; rule provenance; decoy construction method; split assignment",
                "negative_logic": "contrastive learning only",
                "current_use": "LTR negative with clear label_kind=hard_decoy",
            },
            {
                "round": 10,
                "sample_kind": "unlabeled_or_unknown",
                "can_train_as": "unlabeled",
                "allowed_when": "Evidence is absent, ambiguous, family-only, or not exact.",
                "forbidden_when": "Converted to zero/negative label.",
                "required_fields": "reason for unknown; attempted sources; date/version of search",
                "negative_logic": "absence is not negative",
                "current_use": "hold out from supervised positive/negative labels",
            },
            {
                "round": 10,
                "sample_kind": "review_or_context_evidence",
                "can_train_as": "weak_feature_or_search_pointer",
                "allowed_when": "Review or context source names a pathway/family but not exact row-level reaction.",
                "forbidden_when": "Used as exact substrate-product positive.",
                "required_fields": "source; pathway/family; downstream exact-extraction query",
                "negative_logic": "not a negative source",
                "current_use": "prioritization only",
            },
        ]
    )


def build_pipeline() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "round": 10,
                "stage_order": 1,
                "stage": "collect",
                "input": "local gap tables + challenge panel + PubMed/Rhea/GitHub queries",
                "output": "pubmed_candidate_screen_round10.csv; rhea_reaction_screen_round10.csv",
                "pass_gate": "query returns a source that can be inspected",
                "fail_action": "record no-hit; revise query; do not create a negative",
            },
            {
                "round": 10,
                "stage_order": 2,
                "stage": "verify_exactness",
                "input": "candidate papers/database reactions",
                "output": "exact-pair manifest for next round",
                "pass_gate": "exact substrate, product, direction, organism/enzyme/context are present",
                "fail_action": "mark review/context/unlabeled",
            },
            {
                "round": 10,
                "stage_order": 3,
                "stage": "map_chemistry",
                "input": "exact candidate rows",
                "output": "canonical SMILES/InChIKey mapping table",
                "pass_gate": "local and external identifiers agree or mismatch is manually resolved",
                "fail_action": "manual structure review; no training import",
            },
            {
                "round": 10,
                "stage_order": 4,
                "stage": "split_and_holdout_guard",
                "input": "mapped exact rows",
                "output": "training_allowed flag and split assignment",
                "pass_gate": "no challenge/holdout/source/scaffold leakage",
                "fail_action": "hold out or keep as evaluation challenge",
            },
            {
                "round": 10,
                "stage_order": 5,
                "stage": "generator_recall_check",
                "input": "approved exact rows + source-traceable rules",
                "output": "candidate-generation hit/miss table",
                "pass_gate": "deployment generator enumerates the known product",
                "fail_action": "recover/import/fix source-traceable rule before ranker training",
            },
            {
                "round": 10,
                "stage_order": 6,
                "stage": "freeze_and_review",
                "input": "curated CSVs and generated candidates",
                "output": "commit, diff, model-audit report",
                "pass_gate": "all labels have evidence class and source fields",
                "fail_action": "do not retrain; return to verification",
            },
            {
                "round": 10,
                "stage_order": 7,
                "stage": "retrain_after_gates",
                "input": "approved labels and generator outputs",
                "output": "leakage-controlled model metrics",
                "pass_gate": "evaluation sample size and split design support the claim",
                "fail_action": "report prototype-only claim boundary",
            },
        ]
    )


def build_next_queue() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "round": 10,
                "priority_rank": 1,
                "priority_family": "polyphenol_ring_fission",
                "queue_item": "extract exact catechin/epicatechin/flavanone/flavonol ring-fission pairs",
                "seed_sources": "PubMed screen: catechin_valerolactone, quercetin_ring_fission",
                "candidate_source_ids": "PMID:40597460; PMID:40543124; PMID:39674423; PMID:38975869; PMID:37569640",
                "expected_output_next_round": "exact_pair_screen_round11.csv",
                "acceptance_gate": "paper/table states exact substrate and valerolactone/phenylacid product with fermentation or microbial context",
                "reject_if": "review-only or metabolite trend without exact substrate-product edge",
            },
            {
                "round": 10,
                "priority_rank": 2,
                "priority_family": "urolithin_dehydroxylation",
                "queue_item": "recover source-traceable rule for known gold urolithin C -> urolithin A",
                "seed_sources": "Round8/Round9 exact pair; PubMed urolithin screen; Rhea no-hit",
                "candidate_source_ids": "PMID:39856097; PMID:41797252; pair_key:HHXMEXZVPJFAIJ__RIUPLDUFZCXCHM",
                "expected_output_next_round": "source_rule_recovery_round11.csv",
                "acceptance_gate": "source supports exact dehydroxylation direction and generator enumerates urolithin A",
                "reject_if": "urolithin family evidence without exact urolithin C -> A step",
            },
            {
                "round": 10,
                "priority_rank": 3,
                "priority_family": "hydroxycinnamate_reduction_and_hydrolysis",
                "queue_item": "separate hydrolysis positives from reduction candidates",
                "seed_sources": "Rhea caffeic/dihydrocaffeic screens; PubMed caffeic_dihydrocaffeic screen",
                "candidate_source_ids": "RHEA:20689; RHEA:79503; PMID:27977191; PMID:40528807",
                "expected_output_next_round": "hydroxycinnamate_exactness_round11.csv",
                "acceptance_gate": "exact edge and direction; caffeate reductions must be direct and gut/microbial-context supported",
                "reject_if": "rosmarinic-acid pathway support is used as direct caffeic acid -> dihydrocaffeic acid label",
            },
            {
                "round": 10,
                "priority_rank": 4,
                "priority_family": "isoflavone_reductive_and_glycoside_metabolism",
                "queue_item": "screen daidzein/equol and isoflavone glycoside steps separately",
                "seed_sources": "Rhea daidzein/equol screens; PubMed daidzein_equol screen",
                "candidate_source_ids": "RHEA:69683; RHEA:78831; RHEA:78835; PubMed query daidzein_equol",
                "expected_output_next_round": "isoflavone_exactness_round11.csv",
                "acceptance_gate": "stepwise exact pair with organism/enzyme context; equol no-hit in Rhea is not a negative",
                "reject_if": "equol production capacity study without substrate-product reaction row",
            },
            {
                "round": 10,
                "priority_rank": 5,
                "priority_family": "lignan_redox_and_deglycosylation",
                "queue_item": "recover generator rule for known gold pinoresinol -> lariciresinol",
                "seed_sources": "Round8/Round9 exact pair; PubMed pinoresinol_lariciresinol screen",
                "candidate_source_ids": "PMID:12736449; pair_key:HGXBRUKMWQGOIE__MHXCIKYXNYCMHY",
                "expected_output_next_round": "lignan_rule_recovery_round11.csv",
                "acceptance_gate": "source exactness plus deployment generator hit",
                "reject_if": "food lignan content database row without microbial transformation",
            },
            {
                "round": 10,
                "priority_rank": 6,
                "priority_family": "prenylflavonoid_o_demethylation",
                "queue_item": "recover generator rule for known gold isoxanthohumol -> 8-prenylnaringenin",
                "seed_sources": "Round8/Round9 exact pair; PubMed isoxanthohumol_8pn screen",
                "candidate_source_ids": "PMID:16772450; pair_key:YKGCBLWILMDSAV__LPEPZZAVFJPLNZ",
                "expected_output_next_round": "prenylflavonoid_rule_recovery_round11.csv",
                "acceptance_gate": "exact O-demethylation/source context and generator hit",
                "reject_if": "hop phenol pharmacokinetic association without exact conversion",
            },
            {
                "round": 10,
                "priority_rank": 7,
                "priority_family": "glycoside_and_hmo_hydrolysis",
                "queue_item": "resolve 2'-FL structure mismatch and recover hydrolysis generator rule",
                "seed_sources": "Round8/Round9 exact pair; PubMed fucosyllactose_fucose screen; Rhea fucosyllactose related biosynthesis",
                "candidate_source_ids": "PMID:31138818; PMID:33330584; pair_key:SNFSYLYCDAVZGP__SHZGCJCMOBCMKK",
                "expected_output_next_round": "hmo_mapping_and_rule_recovery_round11.csv",
                "acceptance_gate": "manual structure resolution plus exact hydrolysis/product edge and generator hit",
                "reject_if": "Rhea fucosyltransferase biosynthesis is used as hydrolysis evidence",
            },
            {
                "round": 10,
                "priority_rank": 8,
                "priority_family": "bile_acid_deconjugation_and_lipid_context",
                "queue_item": "use Rhea-supported bile deconjugation as source-rule benchmark",
                "seed_sources": "Rhea bile_salt_hydrolase screen",
                "candidate_source_ids": "RHEA:16309; RHEA:19353",
                "expected_output_next_round": "bile_rule_benchmark_round11.csv",
                "acceptance_gate": "source reaction maps to local conjugated bile acid substrate/product and organism context is captured",
                "reject_if": "bile acid pathway context is used as universal positive without exact pair",
            },
            {
                "round": 10,
                "priority_rank": 9,
                "priority_family": "amino_acid_catabolism_and_redox",
                "queue_item": "screen C-module low-recall amino-acid reactions against Rhea/ECREACT before expansion",
                "seed_sources": "local fullrule low-recall categories; future Rhea/ECREACT screen",
                "candidate_source_ids": "local categories: methionine/cysteine; valine/leucine/isoleucine; deamination; redox",
                "expected_output_next_round": "amino_acid_source_rule_screen_round11.csv",
                "acceptance_gate": "food-gut relevance plus exact source reaction and split-safe labels",
                "reject_if": "generic central-metabolism reaction is imported without food-gut scope",
            },
        ]
    )


def build_iteration_log() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "round": 10,
                "iteration": 1,
                "step": "local_generation_gap_triage",
                "evidence_checked": "fullrule_hit_miss.parquet, production_challenge_panel.csv, positive_sample_import_round9.csv",
                "finding": "Reaction family gaps are real, but Round9 shows known gold positives can also be missed by clean candidate generation.",
                "decision": "Prioritize generator recall and source-traceable rule recovery before ranker-only retraining.",
            },
            {
                "round": 10,
                "iteration": 2,
                "step": "external_model_standard_check",
                "evidence_checked": "ECREACT, RetroPathRL/RetroRules, gapseq, gutSMASH GitHub documentation",
                "finding": "Mature systems keep source/EC/reaction database provenance and treat rules as generators or context, not automatic labels.",
                "decision": "Use exact source reactions for positives; use rules as candidate generation with provenance.",
            },
            {
                "round": 10,
                "iteration": 3,
                "step": "life_science_database_screen",
                "evidence_checked": "PubMed esearch/esummary and Rhea reaction search for priority families",
                "finding": "Rhea strongly supports some core reactions such as phenolic acid/bile deconjugation families, but is sparse for urolithin, equol, HMO hydrolysis, and flavanol ring-fission targets.",
                "decision": "Use PubMed/full-text extraction for gut-specific exact pairs; do not treat Rhea no-hit as a negative.",
            },
            {
                "round": 10,
                "iteration": 4,
                "step": "positive_negative_policy",
                "evidence_checked": "Round7/Round8/Round9 outcomes and external standards",
                "finding": "Negatives are the hardest part; most safe negatives are conditional assay negatives or ranking decoys, not global biological negatives.",
                "decision": "Separate exact positives, source rules, assay negatives, conditional negatives, hard decoys, and unlabeled rows in future CSVs.",
            },
        ]
    )


def main() -> None:
    CURATION.mkdir(parents=True, exist_ok=True)
    tables = {
        OUT_GAPS: build_gap_priorities(),
        OUT_STANDARDS: build_external_standards(),
        OUT_PUBMED: build_pubmed_screen(),
        OUT_RHEA: build_rhea_screen(),
        OUT_POLICY: build_policy(),
        OUT_PIPELINE: build_pipeline(),
        OUT_QUEUE: build_next_queue(),
        OUT_LOG: build_iteration_log(),
    }
    for path, table in tables.items():
        table.to_csv(path, index=False)
        print(f"wrote {path} rows={len(table)}")
    print("round10_priority_families=", tables[OUT_GAPS]["priority_family"].tolist())


if __name__ == "__main__":
    main()
