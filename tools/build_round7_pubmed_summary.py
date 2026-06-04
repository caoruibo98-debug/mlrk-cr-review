#!/usr/bin/env python
"""Summarize Round7 PubMed identity checks into a reviewable CSV."""

from __future__ import annotations

import json
import re
from pathlib import Path

import pandas as pd


REPO = Path(__file__).resolve().parents[1]
CURATION = REPO / "data" / "curation"
RAW = CURATION / "ncbi_pubmed_esummary_round7_raw.json"
QUEUE = CURATION / "exact_evidence_extraction_queue_round7.csv"
OUT = CURATION / "literature_pubmed_summary_round7.csv"


PMID_RE = re.compile(r"PMID:(\d+)|\b(\d{7,9})\b")
DOI_RE = re.compile(r"10\.\d{4,9}/[^\s;]+", re.I)


def extract_pmids(source_id: object) -> list[str]:
    if source_id is None or pd.isna(source_id):
        return []
    pmids = []
    for match in PMID_RE.finditer(str(source_id)):
        pmid = match.group(1) or match.group(2)
        if pmid and pmid not in pmids:
            pmids.append(pmid)
    return pmids


def extract_dois(source_id: object) -> list[str]:
    if source_id is None or pd.isna(source_id):
        return []
    return [x.rstrip(".") for x in DOI_RE.findall(str(source_id))]


def title_signal(title: str, substrate: str, product: str) -> str:
    title_l = title.lower()
    sub_l = str(substrate).lower()
    prod_l = str(product).lower()
    has_sub = sub_l in title_l
    has_prod = prod_l in title_l
    if has_sub and has_prod:
        return "title_names_both_compounds"
    if has_sub:
        return "title_names_substrate_only"
    if has_prod:
        return "title_names_product_only"
    if any(word in title_l for word in ["metabolism", "biotransformation", "enzyme", "microbiota", "microflora"]):
        return "title_supports_relevant_metabolism_context"
    return "title_identity_only"


def load_pubmed() -> dict[str, dict[str, object]]:
    raw = json.loads(RAW.read_text(encoding="utf-8"))
    result = raw.get("result", {})
    return {
        str(uid): result.get(str(uid), {})
        for uid in result.get("uids", [])
        if isinstance(result.get(str(uid), {}), dict)
    }


def main() -> None:
    queue = pd.read_csv(QUEUE)
    pubmed = load_pubmed()
    rows: list[dict[str, object]] = []

    for item in queue.itertuples(index=False):
        pmids = extract_pmids(item.source_id)
        dois = extract_dois(item.source_id)
        if not pmids:
            rows.append(
                {
                    "round": 7,
                    "pair_key": item.pair_key,
                    "module": item.module,
                    "substrate_name": item.substrate_name,
                    "product_name": item.product_name,
                    "source_id": item.source_id,
                    "pmid": "",
                    "doi": ";".join(dois),
                    "pubmed_found": False,
                    "pubdate": "",
                    "journal": "",
                    "title": "",
                    "title_signal": "no_pmid_to_check",
                    "identity_check_status": "non_pubmed_or_search_placeholder",
                    "exact_pair_proven": False,
                    "why_not_training_yet": "Need exact source extraction and compound mapping before this can become a positive label.",
                }
            )
            continue

        for pmid in pmids:
            rec = pubmed.get(pmid, {})
            title = str(rec.get("title", ""))
            article_dois = []
            for article_id in rec.get("articleids", []) or []:
                if isinstance(article_id, dict) and str(article_id.get("idtype", "")).lower() == "doi":
                    article_dois.append(str(article_id.get("value", "")))
            rows.append(
                {
                    "round": 7,
                    "pair_key": item.pair_key,
                    "module": item.module,
                    "substrate_name": item.substrate_name,
                    "product_name": item.product_name,
                    "source_id": item.source_id,
                    "pmid": pmid,
                    "doi": ";".join(dois or article_dois),
                    "pubmed_found": bool(rec),
                    "pubdate": rec.get("pubdate", ""),
                    "journal": rec.get("fulljournalname", rec.get("source", "")),
                    "title": title,
                    "title_signal": title_signal(title, item.substrate_name, item.product_name),
                    "identity_check_status": "pubmed_identity_verified" if rec else "pubmed_missing_from_round7_esummary",
                    "exact_pair_proven": False,
                    "why_not_training_yet": (
                        "PubMed identity check is not enough. Extract exact substrate/product, organism or enzyme, "
                        "direction, stereochemistry, and assay context from the paper/table before y=1."
                    ),
                }
            )

    out = pd.DataFrame(rows)
    out.to_csv(OUT, index=False)
    print(f"wrote {OUT} rows={len(out)}")
    print(out["identity_check_status"].value_counts(dropna=False).to_string())
    print(out["title_signal"].value_counts(dropna=False).to_string())


if __name__ == "__main__":
    main()
