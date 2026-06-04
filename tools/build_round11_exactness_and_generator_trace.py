#!/usr/bin/env python
"""Round11 exactness and generator-miss trace.

This round starts from the Round10 queue and checks two hard gates:
1. Does the literature/database source support an exact substrate-product edge?
2. If a known gold edge exists, why does the deployment-faithful clean generator
   fail to keep it as a candidate?
"""

from __future__ import annotations

import re
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit import RDLogger


REPO = Path(__file__).resolve().parents[1]
MODULAR = REPO / "modular"
SRC = REPO / "src"
sys.path.insert(0, str(SRC))
sys.path.insert(0, str(MODULAR))

import kio  # noqa: E402
from generate_candidates import run_reactants  # noqa: E402
from ltr_candidates import MACRO2MOD, POOL, ec_class, load_module_rules  # noqa: E402


RDLogger.DisableLog("rdApp.*")

CURATION = REPO / "data" / "curation"
NCBI_XML = CURATION / "ncbi_round11_raw" / "pubmed_exactness_round11.xml"
ROUND10_QUEUE = CURATION / "round10_next_evidence_queue.csv"
ROUND9_IMPORT = CURATION / "positive_sample_import_round9.csv"
REACTIONS = REPO / "outputs" / "modular" / "ltr" / "reactions.parquet"
CLEAN_FULL = REPO / "outputs" / "modular" / "ltr" / "clean_candidates_full.parquet"
CLEAN = REPO / "outputs" / "modular" / "ltr" / "clean_candidates.parquet"
FULLRULE = REPO / "outputs" / "gen_gap" / "fullrule_hit_miss.parquet"
ROUND5_RECOVERY = CURATION / "master_rule_recovery_candidates_round5.csv"

OUT_PUBMED_EXACT = CURATION / "pubmed_exact_pair_screen_round11.csv"
OUT_GENERATOR = CURATION / "known_gold_generator_trace_round11.csv"
OUT_DECISIONS = CURATION / "round11_training_gate_decisions.csv"
OUT_NEXT = CURATION / "round11_next_actions.csv"


TARGETS = [
    {
        "priority_family": "urolithin_dehydroxylation",
        "pair_key": "HHXMEXZVPJFAIJ__RIUPLDUFZCXCHM",
        "substrate_name": "urolithin C",
        "product_name": "urolithin A",
        "primary_pmids": ["39856097"],
        "supporting_pmids": ["41797252"],
        "exact_terms": ["urolithin c", "urolithin a"],
    },
    {
        "priority_family": "glycoside_and_hmo_hydrolysis",
        "pair_key": "SNFSYLYCDAVZGP__SHZGCJCMOBCMKK",
        "substrate_name": "2'-fucosyllactose",
        "product_name": "L-fucose",
        "primary_pmids": ["31138818"],
        "supporting_pmids": ["33330584"],
        "exact_terms": ["2'-fucosyllactose", "fucose"],
    },
    {
        "priority_family": "lignan_redox_and_deglycosylation",
        "pair_key": "HGXBRUKMWQGOIE__MHXCIKYXNYCMHY",
        "substrate_name": "pinoresinol",
        "product_name": "lariciresinol",
        "primary_pmids": ["12736449"],
        "supporting_pmids": [],
        "exact_terms": ["pinoresinol", "lariciresinol"],
    },
    {
        "priority_family": "prenylflavonoid_o_demethylation",
        "pair_key": "YKGCBLWILMDSAV__LPEPZZAVFJPLNZ",
        "substrate_name": "isoxanthohumol",
        "product_name": "8-prenylnaringenin",
        "primary_pmids": ["16772450"],
        "supporting_pmids": [],
        "exact_terms": ["isoxanthohumol", "8-prenylnaringenin"],
    },
    {
        "priority_family": "hydroxycinnamate_reduction_and_hydrolysis",
        "pair_key": "chlorogenate_or_dicaffeoylquinic_acid__caffeate",
        "substrate_name": "chlorogenate / dicaffeoylquinic acids",
        "product_name": "caffeate / caffeic acid",
        "primary_pmids": ["27977191"],
        "supporting_pmids": [],
        "exact_terms": ["hydrolysis", "intestinal microbiota"],
    },
    {
        "priority_family": "hydroxycinnamate_reduction_and_hydrolysis",
        "pair_key": "caffeic_acid__dihydrocaffeic_acid",
        "substrate_name": "caffeic acid",
        "product_name": "dihydrocaffeic acid",
        "primary_pmids": ["40528807"],
        "supporting_pmids": [],
        "exact_terms": ["caffeic acid", "dihydrocaffeic"],
    },
]


def normalize(text: object) -> str:
    return re.sub(r"\s+", " ", "" if text is None else str(text)).strip()


def b1(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    return str(value).split("-", 1)[0]


def table_with_pair_key(path: Path) -> pd.DataFrame:
    df = pd.read_parquet(path)
    df["pair_key"] = df["sb"].astype(str) + "__" + df["pb"].astype(str)
    return df


def parse_pubmed_xml(path: Path) -> dict[str, dict[str, str]]:
    if not path.exists():
        return {}
    root = ET.parse(path).getroot()
    records: dict[str, dict[str, str]] = {}
    for art in root.findall(".//PubmedArticle"):
        pmid = normalize(art.findtext(".//PMID"))
        title = normalize("".join(art.findtext(".//ArticleTitle") or ""))
        abstract_parts = []
        for node in art.findall(".//Abstract/AbstractText"):
            label = node.attrib.get("Label", "")
            text = normalize("".join(node.itertext()))
            if label and text:
                abstract_parts.append(f"{label}: {text}")
            elif text:
                abstract_parts.append(text)
        doi_values = []
        for node in art.findall("./PubmedData/ArticleIdList/ArticleId"):
            if node.attrib.get("IdType", "").lower() == "doi" and node.text:
                doi_values.append(normalize(node.text))
        for node in art.findall("./MedlineCitation/Article/ELocationID"):
            if node.attrib.get("EIdType", "").lower() == "doi" and node.text:
                doi_values.append(normalize(node.text))
        records[pmid] = {
            "pmid": pmid,
            "title": title,
            "journal": normalize(art.findtext(".//Journal/Title")),
            "year": normalize(art.findtext(".//PubDate/Year")),
            "doi": ";".join(sorted(set(doi_values))),
            "abstract": " ".join(abstract_parts),
        }
    return records


def term_hits(text: str, terms: list[str]) -> list[str]:
    lower = text.lower()
    return [term for term in terms if term.lower() in lower]


def exactness_status(target: dict[str, object], record: dict[str, str], role: str) -> str:
    text = f"{record.get('title', '')} {record.get('abstract', '')}".lower()
    pair_key = str(target["pair_key"])
    if pair_key == "caffeic_acid__dihydrocaffeic_acid":
        return "context_pathway_not_direct_pair"
    if pair_key == "chlorogenate_or_dicaffeoylquinic_acid__caffeate":
        return "hydrolysis_family_requires_specific_pair_mapping"
    terms = target["exact_terms"]
    hits = term_hits(text, list(terms))
    if len(hits) == len(terms):
        if role == "primary":
            return "abstract_or_title_exact_pair_candidate"
        return "supporting_exact_or_mechanistic_context"
    if hits:
        return "partial_context_not_exact_pair"
    return "no_exact_pair_terms_in_title_or_abstract"


def is_exact_source_status(status: object) -> bool:
    return str(status) in {
        "abstract_or_title_exact_pair_candidate",
        "supporting_exact_or_mechanistic_context",
    }


def build_pubmed_exact_screen(records: dict[str, dict[str, str]]) -> pd.DataFrame:
    rows = []
    for target in TARGETS:
        for role, pmids in [("primary", target["primary_pmids"]), ("supporting", target["supporting_pmids"])]:
            for pmid in pmids:
                rec = records.get(str(pmid), {})
                text = f"{rec.get('title', '')} {rec.get('abstract', '')}"
                hits = term_hits(text, list(target["exact_terms"]))
                rows.append(
                    {
                        "round": 11,
                        "priority_family": target["priority_family"],
                        "pair_key": target["pair_key"],
                        "substrate_name": target["substrate_name"],
                        "product_name": target["product_name"],
                        "source_role": role,
                        "pmid": pmid,
                        "doi": rec.get("doi", ""),
                        "title": rec.get("title", ""),
                        "journal": rec.get("journal", ""),
                        "year": rec.get("year", ""),
                        "matched_exact_terms": ";".join(hits),
                        "exactness_status_round11": exactness_status(target, rec, role),
                        "training_decision_round11": "not_training_until_structure_source_and_generator_gates_pass",
                        "required_next_check": "full-text/table extraction for direction, organism/enzyme/assay context, structure IDs, license, and split guard",
                        "pubmed_url": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/",
                    }
                )
    return pd.DataFrame(rows)


def load_all_pool_rules() -> dict[str, list[tuple[str, int]]]:
    d = pd.read_csv(POOL, low_memory=False, usecols=["macro_module", "reaction_rule_smarts", "enzyme_ec"])
    d["module"] = d["macro_module"].map(MACRO2MOD)
    d = d.dropna(subset=["module", "reaction_rule_smarts"])
    d = d[d["reaction_rule_smarts"].astype(str).str.contains(">>")]
    d["ecc"] = d["enzyme_ec"].map(ec_class)
    out: dict[str, list[tuple[str, int]]] = {}
    for module, group in d.groupby("module"):
        unique = group.drop_duplicates("reaction_rule_smarts")[["reaction_rule_smarts", "ecc"]]
        out[module] = list(unique.itertuples(index=False, name=None))
    return out


def generator_hit_for_rules(substrate_smiles: str, substrate_block: str, product_block: str, rules: list[tuple[str, int]]) -> dict[str, object]:
    mol = Chem.MolFromSmiles(substrate_smiles)
    if mol is None:
        return {"target_generated": False, "generated_product_count": 0, "target_rule_count": 0}
    molh = Chem.AddHs(mol)
    generated = set()
    target_rule_count = 0
    for smarts, _ecc in rules:
        for product_smiles in run_reactants(smarts, mol, molh):
            product_b1 = b1(kio.smiles_to_inchikey(product_smiles))
            if not product_b1 or product_b1 == substrate_block:
                continue
            generated.add(product_b1)
            if product_b1 == product_block:
                target_rule_count += 1
    return {
        "target_generated": product_block in generated,
        "generated_product_count": len(generated),
        "target_rule_count": target_rule_count,
    }


def build_generator_trace() -> pd.DataFrame:
    round9 = pd.read_csv(ROUND9_IMPORT)
    reactions = table_with_pair_key(REACTIONS)
    clean_full = table_with_pair_key(CLEAN_FULL)
    clean = table_with_pair_key(CLEAN)
    fullrule = pd.read_parquet(FULLRULE)
    fullrule["pair_key"] = fullrule["sb"].astype(str) + "__" + fullrule["pb"].astype(str)
    recovery = pd.read_csv(ROUND5_RECOVERY)
    recovery["pair_key"] = recovery["sb"].astype(str) + "__" + recovery["pb"].astype(str)

    sampled_rules = load_module_rules(np.random.default_rng(0), 800)
    pool_rules = load_all_pool_rules()

    rows = []
    for item in round9.itertuples(index=False):
        key = item.pair_key
        rx = reactions[reactions["pair_key"] == key]
        if rx.empty:
            continue
        r = rx.iloc[0]
        module = str(r["module"])
        sb = str(r["sb"])
        pb = str(r["pb"])
        sampled = generator_hit_for_rules(r["substrate_smiles"], sb, pb, sampled_rules.get(module, []))
        pool = generator_hit_for_rules(r["substrate_smiles"], sb, pb, pool_rules.get(module, []))
        full_hits = fullrule[fullrule["pair_key"] == key]
        rec_hits = recovery[recovery["pair_key"] == key]
        clean_full_hits = clean_full[clean_full["pair_key"] == key]
        clean_hits = clean[clean["pair_key"] == key]

        if len(clean_full_hits) == 0 and bool(full_hits["hit"].max()) and len(rec_hits) > 0:
            cause = "master_or_recovered_rule_can_generate_but_clean_generator_does_not_keep_pair"
        elif len(clean_full_hits) == 0 and pool["target_generated"] and not sampled["target_generated"]:
            cause = "module_rule_sampling_dropout"
        elif len(clean_full_hits) == 0 and not pool["target_generated"]:
            cause = "module_rule_pool_missing_target_rule"
        elif len(clean_full_hits) > 0 and len(clean_hits) == 0:
            cause = "post_generation_clean_filter_or_sampling_drop"
        else:
            cause = "not_a_current_clean_generation_miss"

        rows.append(
            {
                "round": 11,
                "pair_key": key,
                "module": module,
                "substrate_name": item.substrate_name,
                "product_name": item.product_name,
                "reaction_tier": r.get("tier", ""),
                "reaction_is_clean": bool(r.get("is_clean", False)),
                "in_clean_candidates_full": len(clean_full_hits) > 0,
                "in_clean_candidates": len(clean_hits) > 0,
                "sampled_module_rules_target_generated": sampled["target_generated"],
                "sampled_module_rules_generated_product_count": sampled["generated_product_count"],
                "sampled_module_rules_target_rule_count": sampled["target_rule_count"],
                "full_pool_module_rules_target_generated": pool["target_generated"],
                "full_pool_module_rules_generated_product_count": pool["generated_product_count"],
                "full_pool_module_rules_target_rule_count": pool["target_rule_count"],
                "full_master_rule_diagnose_hit": bool(full_hits["hit"].max()) if not full_hits.empty else False,
                "round5_recovery_candidate_count": len(rec_hits),
                "round5_recovery_sources": ";".join(sorted(set(rec_hits.get("rule_source", pd.Series(dtype=str)).dropna().astype(str)))),
                "round5_generated_target_count": int(rec_hits.get("generated_target_product", pd.Series(dtype=bool)).fillna(False).astype(bool).sum()),
                "generator_miss_cause_round11": cause,
                "recommended_action_round11": recommend_generator_action(cause),
            }
        )
    return pd.DataFrame(rows)


def recommend_generator_action(cause: str) -> str:
    if cause == "master_or_recovered_rule_can_generate_but_clean_generator_does_not_keep_pair":
        return "promote source-traceable recovered/master rule into deployment generator after provenance/license/exactness gate"
    if cause == "module_rule_sampling_dropout":
        return "remove random rule sampling for production-critical gold families or stratify by source-supported recovered rules"
    if cause == "module_rule_pool_missing_target_rule":
        return "recover/import source-traceable rule before retraining"
    if cause == "post_generation_clean_filter_or_sampling_drop":
        return "inspect clean filter and max-candidate sampling"
    return "no generator action required"


def build_training_gate_decisions(pubmed: pd.DataFrame, generator: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for target in TARGETS:
        key = target["pair_key"]
        source_rows = pubmed[pubmed["pair_key"].eq(key)]
        exact_count = int(source_rows["exactness_status_round11"].map(is_exact_source_status).sum()) if not source_rows.empty else 0
        gen_rows = generator[generator["pair_key"].eq(key)]
        generator_ready = bool(gen_rows["in_clean_candidates"].max()) if not gen_rows.empty else False
        known_gold_miss = bool((~gen_rows["in_clean_candidates_full"]).max()) if not gen_rows.empty else False
        if str(key).count("__") == 1 and key in set(generator["pair_key"]):
            decision = "do_not_import_duplicate_positive_fix_generator_first"
        elif exact_count > 0:
            decision = "candidate_positive_requires_mapping_split_and_generator_gate"
        else:
            decision = "not_a_training_label_current_evidence"
        rows.append(
            {
                "round": 11,
                "priority_family": target["priority_family"],
                "pair_key": key,
                "substrate_name": target["substrate_name"],
                "product_name": target["product_name"],
                "exact_source_candidate_count": exact_count,
                "generator_ready_in_clean_candidates": generator_ready,
                "known_gold_generation_miss": known_gold_miss,
                "training_allowed_round11": False,
                "training_decision_round11": decision,
                "why": gate_reason(decision, exact_count, generator_ready, known_gold_miss),
            }
        )
    return pd.DataFrame(rows)


def gate_reason(decision: str, exact_count: int, generator_ready: bool, known_gold_miss: bool) -> str:
    if decision == "do_not_import_duplicate_positive_fix_generator_first":
        return "Pair is already an existing gold positive; current blocker is clean candidate generator recall."
    if exact_count > 0 and not generator_ready:
        return "Source may support an exact edge, but chemistry mapping, split guard, and generator recall are not complete."
    if known_gold_miss:
        return "Known positive exists but is not generated in the deployment-faithful clean candidates."
    return "Evidence is partial, context-only, or not exact enough for supervised training."


def build_next_actions(generator: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for row in generator.itertuples(index=False):
        rows.append(
            {
                "round": 11,
                "queue_type": "generator_recall_fix",
                "priority": "high",
                "pair_key": row.pair_key,
                "substrate_name": row.substrate_name,
                "product_name": row.product_name,
                "blocking_cause": row.generator_miss_cause_round11,
                "next_action": row.recommended_action_round11,
                "must_not_do": "do not duplicate-import this as a new positive; do not use unproven rule hits as gold labels",
            }
        )
    rows.append(
        {
            "round": 11,
            "queue_type": "exact_pair_extraction",
            "priority": "high",
            "pair_key": "polyphenol_ring_fission_family",
            "substrate_name": "catechin/epicatechin/naringenin/quercetin families",
            "product_name": "valerolactones and phenylacid ring-fission products",
            "blocking_cause": "high-priority family has low local full-rule recall and mostly literature-context screens",
            "next_action": "extract exact substrate-product rows from full text/tables before any label import",
            "must_not_do": "do not train on review-only or metabolomics-trend rows",
        }
    )
    return pd.DataFrame(rows)


def main() -> None:
    records = parse_pubmed_xml(NCBI_XML)
    pubmed = build_pubmed_exact_screen(records)
    generator = build_generator_trace()
    decisions = build_training_gate_decisions(pubmed, generator)
    next_actions = build_next_actions(generator)

    for path, table in [
        (OUT_PUBMED_EXACT, pubmed),
        (OUT_GENERATOR, generator),
        (OUT_DECISIONS, decisions),
        (OUT_NEXT, next_actions),
    ]:
        table.to_csv(path, index=False)
        print(f"wrote {path} rows={len(table)}")

    print("round11_exact_status=", pubmed["exactness_status_round11"].value_counts().to_dict())
    print("round11_generator_causes=", generator["generator_miss_cause_round11"].value_counts().to_dict())


if __name__ == "__main__":
    main()
