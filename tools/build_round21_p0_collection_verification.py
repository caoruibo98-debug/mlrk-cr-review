#!/usr/bin/env python
"""Round21 P0 collection and verification worklist.

Audit-only. This script turns Round20 reaction-type gaps into a first concrete
collection batch. It does not modify training data or model code.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Any

import pandas as pd
from rdkit import RDLogger


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "curation"
DOC = ROOT / "docs" / "reviews" / "round21_p0_collection_verification.md"
LTR_OUT = ROOT / "outputs" / "modular" / "ltr"

sys.path.insert(0, str(ROOT / "modular"))
sys.path.insert(0, str(ROOT / "src"))
import ltr_build  # noqa: E402


RDLogger.DisableLog("rdApp.*")

RHEA_RAW = OUT / "rhea_round21_raw"
NCBI_RAW = OUT / "ncbi_round21_raw"


def read_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {}
    return json.loads(path.read_text(encoding="utf-8"))


def pair_key(df: pd.DataFrame) -> pd.Series:
    return df["module"].astype(str) + "__" + df["sb"].astype(str) + "__" + df["pb"].astype(str)


def clean_text(value: Any) -> str:
    if pd.isna(value):
        return ""
    return str(value).strip()


def pubmed_title_map() -> dict[str, dict[str, str]]:
    mapping: dict[str, dict[str, str]] = {}
    for path in [
        NCBI_RAW / "gold_dehydroxylation_pmids.json",
        NCBI_RAW / "urolithin_dehydroxylation_esummary.json",
        NCBI_RAW / "legacy_gold_pmids_esummary.json",
    ]:
        payload = read_json(path)
        result = payload.get("summary", {}).get("result", payload.get("result", {}))
        for pmid in result.get("uids", []):
            rec = result.get(str(pmid), {})
            mapping[str(pmid)] = {
                "title": clean_text(rec.get("title")),
                "journal": clean_text(rec.get("source")),
                "pubdate": clean_text(rec.get("pubdate")),
            }
    return mapping


def pubmed_title_looks_relevant(row: pd.Series, title: str) -> bool:
    if not title:
        return True
    text = " ".join(
        [
            title,
            clean_text(row.get("reaction_type")),
            clean_text(row.get("microbe_or_strain")),
        ]
    ).lower()
    title_only = title.lower()
    terms = [
        clean_text(row.get("substrate_name")).lower(),
        clean_text(row.get("product_name")).lower(),
        "urolithin",
        "ellagic",
        "catechol",
        "dehydroxyl",
        "gut",
        "intestinal",
        "bacter",
        "microb",
        "corticosteroid",
    ]
    useful_terms = [term for term in terms if len(term) >= 4]
    if any(term in title_only for term in useful_terms):
        return True
    if "photoreceptor" in title_only or "vision" in title_only:
        return False
    return any(term in text for term in useful_terms)


def rhea_records(path: Path) -> list[dict[str, Any]]:
    payload = read_json(path)
    if "records" in payload:
        return payload.get("records", [])
    return payload.get("results", [])


def build_local_candidate_pool() -> pd.DataFrame:
    pool = ltr_build.load_pool()
    reactions = pd.read_parquet(LTR_OUT / "reactions.parquet")
    clean_full = pd.read_parquet(LTR_OUT / "clean_candidates_full.parquet")
    clean = pd.read_parquet(LTR_OUT / "clean_candidates.parquet")
    for df in [pool, reactions, clean_full, clean]:
        df["_pair_key"] = pair_key(df)
    clean_full_pos = set(clean_full.loc[clean_full["y"].eq(1), "_pair_key"])
    clean_pos = set(clean.loc[clean["y"].eq(1), "_pair_key"])
    rxn_keys = set(reactions["_pair_key"])
    pool = pool.copy()
    pool["tier"] = pool["tier_rank"].map({2: "gold", 1: "silver", 0: "weak"}).fillna("weak")
    pool["in_reactions"] = pool["_pair_key"].isin(rxn_keys)
    pool["in_clean_full_pos"] = pool["_pair_key"].isin(clean_full_pos)
    pool["in_clean_pos"] = pool["_pair_key"].isin(clean_pos)
    pool["_evidence_score"] = pool[
        ["tier_rank", "ev_has_pmid", "ev_has_ec", "ev_has_rule", "ev_microbe_n", "ev_gene_n"]
    ].sum(axis=1)
    pool = pool.sort_values("_evidence_score", ascending=False).drop_duplicates("_pair_key")
    return pool


def p0_selector(pool: pd.DataFrame) -> pd.DataFrame:
    route_lost_gold = pool[
        pool["tier"].eq("gold")
        & pool["in_reactions"]
        & ~pool["in_clean_pos"]
        & (
            pool["reaction_type"].astype(str).str.contains("urolithin|catechol|dehydroxylase|dehydroxylating", case=False, na=False)
            | pool["reaction_category"].astype(str).str.contains("functional_group_removal", case=False, na=False)
        )
    ].copy()
    route_lost_gold = route_lost_gold.sort_values(
        ["in_clean_full_pos", "ev_has_pmid", "pubmed_year"],
        ascending=[True, False, False],
    ).head(14)

    dopamine = pool[
        pool["substrate_name"].astype(str).str.lower().eq("dopamine")
        & pool["product_name"].astype(str).str.lower().isin(["m-tyramine", "3-tyramine"])
    ].head(1)

    rows = pd.concat([route_lost_gold, dopamine], ignore_index=True).drop_duplicates("_pair_key")
    return rows


def candidate_verification(row: pd.Series, pmids: dict[str, dict[str, str]]) -> dict[str, Any]:
    substrate = clean_text(row["substrate_name"]).lower()
    product = clean_text(row["product_name"]).lower()
    pmid = clean_text(row.get("pmid"))
    pmid = pmid[:-2] if pmid.endswith(".0") else pmid
    pubmed = pmids.get(pmid, {}) if pmid else {}
    evidence_sources = []
    verification_status = "blocked_pending_external_verification"
    rhea_id = ""
    rhea_equation = ""
    rhea_status = ""
    if substrate == "dopamine" and product in {"m-tyramine", "3-tyramine"}:
        recs = rhea_records(RHEA_RAW / "RHEA_61520.json")
        if recs:
            rhea_id = f"RHEA:{recs[0].get('id')}"
            rhea_equation = clean_text(recs[0].get("equation"))
            rhea_status = clean_text(recs[0].get("status"))
            evidence_sources.append("Rhea exact approved reaction")
            verification_status = "verified_exact_reaction_rhea_plus_literature"
    elif "urolithin" in substrate or "urolithin" in product:
        urolithin_payload = read_json(RHEA_RAW / "urolithin.json")
        no_rhea = (
            urolithin_payload.get("record_count_available") == 0
            or urolithin_payload.get("count") == 0
        )
        if pubmed:
            evidence_sources.append("PubMed curated gut urolithin literature")
            evidence_sources.append("Rhea no urolithin exact hit" if no_rhea else "Rhea search needs manual resolution")
            verification_status = "literature_verified_rhea_absent_not_training_yet"
        else:
            evidence_sources.append("Rhea no urolithin exact hit" if no_rhea else "Rhea search needs manual resolution")
            verification_status = "manual_literature_check_required"
    elif "hydrocaffeic" in substrate or "tyramine" in product or "catechol" in clean_text(row["reaction_type"]).lower():
        if pubmed or str(pmid) == "32067637":
            evidence_sources.append("PubMed catechol dehydroxylase literature")
            verification_status = "literature_verified_rhea_partial_or_absent"
        else:
            verification_status = "manual_catechol_literature_check_required"

    pubmed_relevance = "not_checked"
    if pmid and pubmed:
        pubmed_relevance = "title_context_pass" if pubmed_title_looks_relevant(row, pubmed.get("title", "")) else "title_context_mismatch"
        if pubmed_relevance == "title_context_mismatch":
            verification_status = "blocked_pubmed_title_mismatch_requires_recuration"
            evidence_sources.append("PubMed title mismatch")

    return {
        "verification_status_round21": verification_status,
        "pubmed_relevance_round21": pubmed_relevance,
        "rhea_id": rhea_id,
        "rhea_equation": rhea_equation,
        "rhea_status": rhea_status,
        "pubmed_pmid": pmid,
        "pubmed_title": pubmed.get("title", ""),
        "pubmed_journal": pubmed.get("journal", ""),
        "pubmed_pubdate": pubmed.get("pubdate", ""),
        "evidence_sources_round21": ";".join(evidence_sources),
    }


def build_collection_worklist() -> pd.DataFrame:
    pmids = pubmed_title_map()
    pool = build_local_candidate_pool()
    chosen = p0_selector(pool)
    rows = []
    for idx, row in chosen.reset_index(drop=True).iterrows():
        ver = candidate_verification(row, pmids)
        route_loss = []
        if bool(row["in_reactions"]) and not bool(row["in_clean_full_pos"]):
            route_loss.append("reactions_to_clean_full_loss")
        if bool(row["in_clean_full_pos"]) and not bool(row["in_clean_pos"]):
            route_loss.append("clean_full_to_final_loss")
        if bool(row["in_reactions"]) and not bool(row["in_clean_pos"]):
            route_loss.append("final_clean_absent")
        rows.append(
            {
                "round": 21,
                "candidate_id": f"R21-P0-{idx + 1:03d}",
                "priority": "P0",
                "module": row["module"],
                "reaction_family": clean_text(row["reaction_category"]) or "missing_reaction_category",
                "reaction_type": clean_text(row["reaction_type"]),
                "substrate_name": clean_text(row["substrate_name"]),
                "product_name": clean_text(row["product_name"]),
                "substrate_smiles": clean_text(row["substrate_smiles"]),
                "product_smiles": clean_text(row["product_smiles"]),
                "sb": clean_text(row["sb"]),
                "pb": clean_text(row["pb"]),
                "pair_key": clean_text(row["_pair_key"]),
                "tier": clean_text(row["tier"]),
                "source_dataset": clean_text(row["source_dataset"]),
                "source_origin_type": clean_text(row["source_origin_type"]),
                "training_use_recommendation": clean_text(row["training_use_recommendation"]),
                "enzyme_ec": clean_text(row.get("enzyme_ec")),
                "microbe_or_strain": clean_text(row.get("microbe_or_strain")),
                "supporting_sentence": clean_text(row.get("supporting_sentence")),
                "in_reactions": bool(row["in_reactions"]),
                "in_clean_full_pos": bool(row["in_clean_full_pos"]),
                "in_clean_pos": bool(row["in_clean_pos"]),
                "route_gap_round21": ";".join(route_loss) if route_loss else "no_route_gap_observed",
                "decision_round21": "route_repair_or_manual_import_gate_not_training",
                "training_allowed_round21": False,
                **ver,
            }
        )
    return pd.DataFrame(rows)


def build_source_queries() -> pd.DataFrame:
    rows = [
        {
            "round": 21,
            "query_id": "R21-Q-RHEA-UROLITHIN",
            "source": "Rhea",
            "query": "urolithin",
            "raw_output": "data/curation/rhea_round21_raw/urolithin.json",
            "result_summary": "0 returned records; use PubMed/literature as primary evidence for gut urolithin dehydroxylation.",
            "next_action": "manual exact participant mapping against PubMed-curated urolithin paper tables",
        },
        {
            "round": 21,
            "query_id": "R21-Q-RHEA-DOPAMINE",
            "source": "Rhea",
            "query": "RHEA:61520; dopamine m-tyramine",
            "raw_output": "data/curation/rhea_round21_raw/RHEA_61520.json",
            "result_summary": "Approved Rhea reaction: dopamine + AH2 = 3-tyramine + A + H2O.",
            "next_action": "use as exact positive verification and generator target test seed",
        },
        {
            "round": 21,
            "query_id": "R21-Q-RHEA-CAFFEATE",
            "source": "Rhea",
            "query": "caffeate dehydroxylase",
            "raw_output": "data/curation/rhea_round21_raw/caffeate_dehydroxylase.json",
            "result_summary": "0 returned records; hydrocaffeic/caffeic dehydroxylation remains literature-led.",
            "next_action": "use PubMed PMID 32067637 and exact structure check before route repair",
        },
        {
            "round": 21,
            "query_id": "R21-Q-PUBMED-UROLITHIN",
            "source": "PubMed",
            "query": "urolithin dehydroxylase Enterocloster gut microbiota",
            "raw_output": "data/curation/ncbi_round21_raw/urolithin_dehydroxylase_esearch.json",
            "result_summary": "Returned PMIDs 41298472, 39856097, 37494568.",
            "next_action": "extract exact substrate/product rows and sentence evidence",
        },
        {
            "round": 21,
            "query_id": "R21-Q-PUBMED-CATECHOL",
            "source": "PubMed",
            "query": "A widely distributed metalloenzyme class catechols",
            "raw_output": "data/curation/ncbi_round21_raw/catechol_metalloenzyme_esearch.json",
            "result_summary": "Returned PMID 32067637.",
            "next_action": "extract substrate panel and exact products for catechol dehydroxylation family",
        },
    ]
    return pd.DataFrame(rows)


def build_route_actions(worklist: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for family, group in worklist.groupby("reaction_family", dropna=False):
        route_losses = sorted(set(";".join(group["route_gap_round21"]).split(";")))
        rows.append(
            {
                "round": 21,
                "reaction_family": family,
                "candidate_count": len(group),
                "verified_or_literature_supported_count": int(
                    group["verification_status_round21"].str.contains("verified|literature", na=False).sum()
                ),
                "route_gap_types": ";".join([x for x in route_losses if x]),
                "next_generator_gate": "exact target route dry-run before any training import",
                "next_negative_gate": "generated decoys must be screened against Rhea/ECReact/EnzymeMap/current positive pool",
                "training_allowed_round21": False,
            }
        )
    return pd.DataFrame(rows)


def write_doc(worklist: pd.DataFrame, queries: pd.DataFrame, actions: pd.DataFrame) -> None:
    DOC.parent.mkdir(parents=True, exist_ok=True)
    text = "# Round21 P0 Collection And Verification\n\n"
    text += "## Working conclusion\n\n"
    text += (
        "Round21 confirms that the next audit should not simply add more rows to training. "
        "The P0 batch contains source-backed positives that are already in `reactions.parquet` but absent from final clean candidates. "
        "The right next move is exact evidence extraction plus route dry-run, then negative screening.\n\n"
    )
    text += "## Candidate worklist summary\n\n```text\n"
    cols = [
        "candidate_id",
        "substrate_name",
        "product_name",
        "tier",
        "verification_status_round21",
        "route_gap_round21",
        "training_allowed_round21",
    ]
    text += worklist[cols].to_string(index=False)
    text += "\n```\n\n"
    text += "## Source queries\n\n```text\n"
    text += queries[["query_id", "source", "result_summary", "next_action"]].to_string(index=False)
    text += "\n```\n\n"
    text += "## Family route actions\n\n```text\n"
    text += actions.to_string(index=False)
    text += "\n```\n\n"
    text += "## Written artifacts\n\n"
    text += "- `data/curation/round21_p0_collection_worklist.csv`\n"
    text += "- `data/curation/round21_source_query_manifest.csv`\n"
    text += "- `data/curation/round21_route_repair_actions.csv`\n"
    DOC.write_text(text, encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    worklist = build_collection_worklist()
    queries = build_source_queries()
    actions = build_route_actions(worklist)
    outputs = {
        "round21_p0_collection_worklist.csv": worklist,
        "round21_source_query_manifest.csv": queries,
        "round21_route_repair_actions.csv": actions,
    }
    for name, df in outputs.items():
        path = OUT / name
        df.to_csv(path, index=False)
        print(f"wrote {path.relative_to(ROOT)} rows={len(df)}")
    write_doc(worklist, queries, actions)
    print(f"wrote {DOC.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
