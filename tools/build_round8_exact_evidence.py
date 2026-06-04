#!/usr/bin/env python
"""Build Round8 exact-pair evidence and label-readiness manifests.

This pass uses PubMed abstracts and Round7 database checks to separate exact
substrate-product evidence from pathway/family evidence. It intentionally keeps
training gates closed until chemical mapping, source extraction, and leakage
guards are complete.
"""

from __future__ import annotations

import re
import xml.etree.ElementTree as ET
from pathlib import Path

import pandas as pd


REPO = Path(__file__).resolve().parents[1]
CURATION = REPO / "data" / "curation"

ROUND7_QUEUE = CURATION / "exact_evidence_extraction_queue_round7.csv"
ROUND7_SOURCE = CURATION / "source_reaction_check_round7.csv"
ROUND7_RHEA = CURATION / "rhea_source_reaction_summary_round7.csv"
PUBMED_XML = CURATION / "ncbi_pubmed_efetch_round8_raw.xml"
PUBMED_EXTRA_XML = CURATION / "ncbi_pubmed_efetch_round8_extra_raw.xml"

OUT_EVIDENCE = CURATION / "exact_pair_evidence_round8.csv"
OUT_POSITIVES = CURATION / "positive_sample_expansion_round8.csv"
OUT_NEGATIVES = CURATION / "negative_evidence_candidates_round8.csv"
OUT_PIPELINE = CURATION / "round8_pipeline_queue.csv"


PMID_RE = re.compile(r"PMID:(\d+)|\b(\d{7,9})\b")
DOI_RE = re.compile(r"10\.\d{4,9}/[^\s;]+", re.I)


MANUAL_ASSESSMENTS = {
    ("QAIPRVGONGVQAS__DZAUWHJDUNRCTF", "40528807"): {
        "abstract_evidence_class": "pathway_product_support_not_exact_pair",
        "exact_pair_supported_by_abstract": False,
        "source_strength": "supports rosmarinic-acid pathway producing DHCA, not caffeic acid to DHCA directly",
        "round8_decision": "do_not_train_positive_current_evidence",
        "next_action": "extract full reaction sequence and decide whether CA->DHCA is a valid intermediate edge",
    },
    ("ATJVZXXHKSYELS__KSEBMYQBYZTDHS", "19502437"): {
        "abstract_evidence_class": "substrate_product_family_support_table_needed",
        "exact_pair_supported_by_abstract": False,
        "source_strength": "ethyl ferulate substrate and ferulic-acid esterase activity are present, but exact product mapping needs table/full text",
        "round8_decision": "fulltext_table_required_before_positive",
        "next_action": "extract enzyme substrate table for ethyl ferulate hydrolysis product and EC/protein context",
    },
    ("HGXBRUKMWQGOIE__MHXCIKYXNYCMHY", "12736449"): {
        "abstract_evidence_class": "abstract_exact_pair_strong",
        "exact_pair_supported_by_abstract": True,
        "source_strength": "abstract explicitly reports Enterococcus faecalis PDG-1 transformation of pinoresinol to lariciresinol",
        "round8_decision": "exact_positive_candidate_pending_mapping_holdout",
        "next_action": "map stereochemistry and source row; keep Round2 holdout guard before training",
    },
    ("MHXCIKYXNYCMHY__PUETUDUXMCLALY", "11453749"): {
        "abstract_evidence_class": "pathway_scheme_fulltext_needed",
        "exact_pair_supported_by_abstract": False,
        "source_strength": "abstract says a metabolic scheme is presented but does not state lariciresinol to secoisolariciresinol directly",
        "round8_decision": "fulltext_scheme_required_before_positive",
        "next_action": "extract pathway scheme/table and reconcile with Rhea direction and stereochemistry",
    },
    ("RIUPLDUFZCXCHM__WXUQMTRHPNOXBV", "29583112"): {
        "abstract_evidence_class": "organism_product_context_not_exact_pair",
        "exact_pair_supported_by_abstract": False,
        "source_strength": "abstract supports isourolithin A production from ellagic acid, not urolithin A to urolithin B",
        "round8_decision": "do_not_train_positive_current_evidence",
        "next_action": "search exact urolithin A to B source or mark as candidate only",
    },
    ("HHXMEXZVPJFAIJ__RIUPLDUFZCXCHM", "39856097"): {
        "abstract_evidence_class": "abstract_exact_pair_strong",
        "exact_pair_supported_by_abstract": True,
        "source_strength": "abstract identifies urolithin C precursor and Enterocloster ucd-dependent production of urolithin A",
        "round8_decision": "exact_positive_candidate_pending_mapping_holdout",
        "next_action": "extract Ucd assay/ex vivo details, position-specific dehydroxylation, and organism context",
    },
    ("HHXMEXZVPJFAIJ__WDGSXHQNUPZEHA", "29583112"): {
        "abstract_evidence_class": "organism_product_context_not_exact_pair",
        "exact_pair_supported_by_abstract": False,
        "source_strength": "abstract supports isourolithin A production from ellagic acid, not urolithin C to isourolithin A",
        "round8_decision": "do_not_train_positive_current_evidence",
        "next_action": "search exact intermediate route or keep as candidate only",
    },
    ("SNFSYLYCDAVZGP__SHZGCJCMOBCMKK", "31138818"): {
        "abstract_evidence_class": "abstract_exact_pair_strong",
        "exact_pair_supported_by_abstract": True,
        "source_strength": "abstract reports fucose cleavage from 2'-fucosyllactose by B. infantis Bi-26 metabolism",
        "round8_decision": "exact_positive_candidate_pending_mapping_holdout",
        "next_action": "confirm L-fucose chemical mapping and whether product edge should include lactose/monomers",
    },
    ("SNFSYLYCDAVZGP__SHZGCJCMOBCMKK", "19520709"): {
        "abstract_evidence_class": "enzyme_family_support_with_conditional_negative",
        "exact_pair_supported_by_abstract": False,
        "source_strength": "abstract characterizes fucosidases and includes context-specific non-activity for AfcB on alpha1,2-fucosyl substrates",
        "round8_decision": "context_evidence_not_global_positive",
        "next_action": "use for enzyme-context annotation; do not duplicate positive row without exact 2-FL assay",
    },
    ("YKGCBLWILMDSAV__LPEPZZAVFJPLNZ", "20397197"): {
        "abstract_evidence_class": "abstract_exact_pair_moderate",
        "exact_pair_supported_by_abstract": True,
        "source_strength": "abstract confirms microbial formation of 8-prenylnaringenin after incubation of xanthohumol and isoxanthohumol with fecal slurries",
        "round8_decision": "exact_positive_candidate_better_source_available",
        "next_action": "prefer PMID:16772450 for direct isoxanthohumol conversion evidence",
    },
    ("YKGCBLWILMDSAV__LPEPZZAVFJPLNZ", "41571201"): {
        "abstract_evidence_class": "health_effect_context_not_reaction_assay",
        "exact_pair_supported_by_abstract": False,
        "source_strength": "abstract calls 8-PN an isoxanthohumol metabolite but focuses on endothelial function",
        "round8_decision": "do_not_use_as_primary_training_source",
        "next_action": "keep as biological context only; use reaction-focused source for positive label",
    },
    ("YKGCBLWILMDSAV__LPEPZZAVFJPLNZ", "16772450"): {
        "abstract_evidence_class": "abstract_exact_pair_strong",
        "exact_pair_supported_by_abstract": True,
        "source_strength": "abstract reports intestinal activation/conversion of isoxanthohumol into 8-prenylnaringenin in vitro and in humans",
        "round8_decision": "exact_positive_candidate_pending_mapping_holdout",
        "next_action": "use this reaction-focused source preferentially; extract SHIME/fecal/human context and mapping",
    },
    ("REFJWTPEDVJJIY__CXQWRCVTCMQVQX", "36432010"): {
        "abstract_evidence_class": "mechanism_not_target_pair",
        "exact_pair_supported_by_abstract": False,
        "source_strength": "abstract studies taxifolin binding and conversion to alphitonin, not quercetin to taxifolin",
        "round8_decision": "reject_as_positive_current_evidence",
        "next_action": "remove from positive expansion until exact quercetin-to-taxifolin source exists",
    },
    ("IBHWREHFNDMRPR__QCDYQQDYXPDABM", ""): {
        "abstract_evidence_class": "no_pubmed_or_rhea_evidence_found",
        "exact_pair_supported_by_abstract": False,
        "source_strength": "Round8 PubMed and Rhea narrow searches returned zero records",
        "round8_decision": "do_not_train_positive_current_evidence",
        "next_action": "search broader decarboxylase/pathway names or keep as external curation gap",
    },
}


EXTRA_SOURCES = [
    {
        "pair_key": "YKGCBLWILMDSAV__LPEPZZAVFJPLNZ",
        "module": "A",
        "gap_family": "prenylflavonoid_o_demethylation",
        "substrate_name": "isoxanthohumol",
        "product_name": "8-prenylnaringenin",
        "source_db": "PubMed",
        "source_id": "PMID:16772450",
        "source_title": "The prenylflavonoid isoxanthohumol from hops is activated into 8-prenylnaringenin in vitro and in the human intestine.",
        "evidence_type": "in_vitro_and_human_intestinal_conversion",
        "source_reaction_check_status": "external_exact_literature_source_added_round8",
        "round7_import_recommendation": "reaction_focused_replacement_source_for_exact_pair",
    }
]


def pmids_from_text(value: object) -> list[str]:
    if value is None or pd.isna(value):
        return []
    out: list[str] = []
    for match in PMID_RE.finditer(str(value)):
        pmid = match.group(1) or match.group(2)
        if pmid and pmid not in out:
            out.append(pmid)
    return out


def dois_from_text(value: object) -> list[str]:
    if value is None or pd.isna(value):
        return []
    return [x.rstrip(".") for x in DOI_RE.findall(str(value))]


def normalize_text(text: str) -> str:
    return re.sub(r"\s+", " ", text or "").strip()


def parse_pubmed_xml(paths: list[Path]) -> dict[str, dict[str, str]]:
    records: dict[str, dict[str, str]] = {}
    for path in paths:
        if not path.exists():
            continue
        root = ET.parse(path).getroot()
        for art in root.findall(".//PubmedArticle"):
            pmid = art.findtext(".//PMID") or ""
            title_node = art.find(".//ArticleTitle")
            title = normalize_text("".join(title_node.itertext())) if title_node is not None else ""
            abstract = normalize_text(
                " ".join("".join(node.itertext()) for node in art.findall(".//AbstractText"))
            )
            journal = art.findtext(".//Article/Journal/Title") or art.findtext(".//MedlineTA") or ""
            year = (
                art.findtext(".//Article/Journal/JournalIssue/PubDate/Year")
                or art.findtext(".//Article/ArticleDate/Year")
                or ""
            )
            doi_values = []
            for node in art.findall("./PubmedData/ArticleIdList/ArticleId"):
                if node.attrib.get("IdType", "").lower() == "doi" and node.text:
                    doi_values.append(node.text)
            for node in art.findall("./MedlineCitation/Article/ELocationID"):
                if node.attrib.get("EIdType", "").lower() == "doi" and node.text:
                    doi_values.append(node.text)
            records[pmid] = {
                "pmid": pmid,
                "title": title,
                "abstract": abstract,
                "journal": journal,
                "year": year,
                "doi": ";".join(sorted(set(doi_values))),
            }
    return records


def abstract_term_flags(abstract: str, substrate: str, product: str) -> dict[str, bool | str]:
    a = abstract.lower()
    sub = str(substrate).lower()
    prod = str(product).lower()
    gut_terms = ["gut", "fecal", "faecal", "intestinal", "microbiota", "microflora", "bifidobacterium", "enterocloster", "enterococcus", "lactobacillus"]
    enzyme_terms = ["enzyme", "esterase", "fucosidase", "dehydroxylase", "reductase", "hydrolase", "operon"]
    negative_terms = ["did not act", "remained stable", "no conversion", "not act", "no activity", "only samples"]
    return {
        "abstract_contains_substrate_name": sub in a,
        "abstract_contains_product_name": prod in a,
        "abstract_contains_gut_or_microbe_context": any(term in a for term in gut_terms),
        "abstract_contains_enzyme_or_gene_context": any(term in a for term in enzyme_terms),
        "abstract_contains_negative_or_context_limit": any(term in a for term in negative_terms),
        "abstract_character_count": len(abstract),
    }


def assessment_for(pair_key: str, pmid: str) -> dict[str, object]:
    if (pair_key, pmid) in MANUAL_ASSESSMENTS:
        return MANUAL_ASSESSMENTS[(pair_key, pmid)]
    if (pair_key, "") in MANUAL_ASSESSMENTS:
        return MANUAL_ASSESSMENTS[(pair_key, "")]
    return {
        "abstract_evidence_class": "not_assessed",
        "exact_pair_supported_by_abstract": False,
        "source_strength": "No manual Round8 assessment available.",
        "round8_decision": "manual_review_required",
        "next_action": "assess exact pair evidence",
    }


def build_evidence_rows() -> pd.DataFrame:
    queue = pd.read_csv(ROUND7_QUEUE)
    extras = pd.DataFrame(EXTRA_SOURCES)
    rows_in = pd.concat([queue, extras], ignore_index=True, sort=False)
    pubmed = parse_pubmed_xml([PUBMED_XML, PUBMED_EXTRA_XML])
    source = pd.read_csv(ROUND7_SOURCE)
    rhea = pd.read_csv(ROUND7_RHEA)

    source_status = {
        (row.substrate_name, row.product_name): row.source_reaction_check_status
        for row in source.itertuples(index=False)
    }
    rhea_exact = {
        (row.substrate_name, row.product_name): row.rhea_pair_exactness_by_name
        for row in rhea.itertuples(index=False)
        if str(row.rhea_pair_exactness_by_name).startswith("rhea_equation")
    }

    rows: list[dict[str, object]] = []
    for item in rows_in.itertuples(index=False):
        pair_key = getattr(item, "pair_key", "")
        source_id = getattr(item, "source_id", "")
        pmids = pmids_from_text(source_id)
        if not pmids:
            pmids = [""]
        for pmid in pmids:
            rec = pubmed.get(pmid, {})
            abstract = rec.get("abstract", "")
            assessment = assessment_for(pair_key, pmid)
            flags = abstract_term_flags(abstract, item.substrate_name, item.product_name)
            rows.append(
                {
                    "round": 8,
                    "pair_key": pair_key,
                    "module": getattr(item, "module", ""),
                    "gap_family": getattr(item, "gap_family", ""),
                    "substrate_name": item.substrate_name,
                    "product_name": item.product_name,
                    "source_db": getattr(item, "source_db", ""),
                    "source_id": source_id,
                    "pmid": pmid,
                    "doi": rec.get("doi", ";".join(dois_from_text(source_id))),
                    "year": rec.get("year", ""),
                    "journal": rec.get("journal", ""),
                    "title": rec.get("title", getattr(item, "source_title", "")),
                    "pubmed_xml_found": bool(rec),
                    "source_reaction_check_status": getattr(
                        item,
                        "source_reaction_check_status",
                        source_status.get((item.substrate_name, item.product_name), ""),
                    ),
                    "rhea_pair_status": rhea_exact.get((item.substrate_name, item.product_name), ""),
                    **flags,
                    "abstract_evidence_class": assessment["abstract_evidence_class"],
                    "exact_pair_supported_by_abstract": assessment["exact_pair_supported_by_abstract"],
                    "source_strength_summary": assessment["source_strength"],
                    "round8_decision": assessment["round8_decision"],
                    "training_allowed_round8": False,
                    "why_training_still_blocked": (
                        "Round8 uses abstract-level evidence screening. Training still requires exact source extraction, "
                        "chemical identifier mapping, direction/stereochemistry check, license note, and split/holdout guard."
                    ),
                    "next_action": assessment["next_action"],
                }
            )
    return pd.DataFrame(rows)


def build_positive_candidates(evidence: pd.DataFrame) -> pd.DataFrame:
    candidates = evidence[
        evidence["round8_decision"].isin(
            [
                "exact_positive_candidate_pending_mapping_holdout",
            ]
        )
    ].copy()
    if candidates.empty:
        return pd.DataFrame()
    candidates["round8_positive_status"] = "verified_exact_positive_candidate_not_training_yet"
    candidates["required_before_import"] = (
        "substrate/product InChIKey and SMILES mapping; exact source extraction; direction; stereochemistry; "
        "organism/enzyme/assay context; source license; holdout/split guard"
    )
    candidates["label_type_after_all_gates"] = "positive"
    return candidates[
        [
            "round",
            "pair_key",
            "module",
            "gap_family",
            "substrate_name",
            "product_name",
            "pmid",
            "doi",
            "title",
            "abstract_evidence_class",
            "source_strength_summary",
            "round8_positive_status",
            "training_allowed_round8",
            "required_before_import",
            "label_type_after_all_gates",
            "next_action",
        ]
    ]


def build_negative_candidates(evidence: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for item in evidence.itertuples(index=False):
        if not item.abstract_contains_negative_or_context_limit:
            continue
        if item.pmid == "19520709":
            rows.append(
                {
                    "round": 8,
                    "negative_candidate_type": "conditional_assay_negative",
                    "pair_key": item.pair_key,
                    "pmid": item.pmid,
                    "title": item.title,
                    "context": "AfcB alpha-L-fucosidase non-activity on alpha1,2-fucosyl substrates; not a global 2-FL negative",
                    "can_train_global_negative": False,
                    "can_train_context_negative": True,
                    "required_before_use": "extract enzyme name, substrate class, exact assay substrate, conditions, and detection statement",
                }
            )
        elif item.pmid == "40528807":
            rows.append(
                {
                    "round": 8,
                    "negative_candidate_type": "conditional_context_limit",
                    "pair_key": item.pair_key,
                    "pmid": item.pmid,
                    "title": item.title,
                    "context": "Rosmarinic acid stability in upper gut and antibiotic-treated chicken cecum; not a caffeic-acid negative",
                    "can_train_global_negative": False,
                    "can_train_context_negative": True,
                    "required_before_use": "extract exact compartment/antibiotic context and substrate tested",
                }
            )
    if not rows:
        return pd.DataFrame(
            columns=[
                "round",
                "negative_candidate_type",
                "pair_key",
                "pmid",
                "title",
                "context",
                "can_train_global_negative",
                "can_train_context_negative",
                "required_before_use",
            ]
        )
    return pd.DataFrame(rows)


def build_pipeline(evidence: pd.DataFrame, positives: pd.DataFrame, negatives: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for item in positives.itertuples(index=False):
        rows.append(
            {
                "round": 8,
                "queue_type": "positive_exact_extraction",
                "priority": "high",
                "pair_key": item.pair_key,
                "substrate_name": item.substrate_name,
                "product_name": item.product_name,
                "source_id": f"PMID:{item.pmid}" if item.pmid else "",
                "next_action": item.next_action,
            }
        )
    for item in evidence.itertuples(index=False):
        if str(item.round8_decision).startswith("fulltext") or item.round8_decision == "fulltext_table_required_before_positive":
            rows.append(
                {
                    "round": 8,
                    "queue_type": "fulltext_table_extraction",
                    "priority": "medium",
                    "pair_key": item.pair_key,
                    "substrate_name": item.substrate_name,
                    "product_name": item.product_name,
                    "source_id": f"PMID:{item.pmid}" if item.pmid else "",
                    "next_action": item.next_action,
                }
            )
    for item in negatives.itertuples(index=False):
        rows.append(
            {
                "round": 8,
                "queue_type": "conditional_negative_extraction",
                "priority": "medium",
                "pair_key": item.pair_key,
                "substrate_name": "",
                "product_name": "",
                "source_id": f"PMID:{item.pmid}" if item.pmid else "",
                "next_action": item.required_before_use,
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    evidence = build_evidence_rows()
    positives = build_positive_candidates(evidence)
    negatives = build_negative_candidates(evidence)
    pipeline = build_pipeline(evidence, positives, negatives)

    evidence.to_csv(OUT_EVIDENCE, index=False)
    positives.to_csv(OUT_POSITIVES, index=False)
    negatives.to_csv(OUT_NEGATIVES, index=False)
    pipeline.to_csv(OUT_PIPELINE, index=False)

    print(f"wrote {OUT_EVIDENCE} rows={len(evidence)}")
    print(evidence["round8_decision"].value_counts(dropna=False).to_string())
    print(f"wrote {OUT_POSITIVES} rows={len(positives)}")
    print(f"wrote {OUT_NEGATIVES} rows={len(negatives)}")
    print(f"wrote {OUT_PIPELINE} rows={len(pipeline)}")


if __name__ == "__main__":
    main()
