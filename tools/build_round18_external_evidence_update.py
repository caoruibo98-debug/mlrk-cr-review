from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "curation"
DOC = ROOT / "docs" / "reviews" / "training_gate_round18.md"
RHEA_RAW = OUT / "rhea_round18_raw"
NCBI_RAW = OUT / "ncbi_round18_raw"


EXTERNAL_SOURCE_CONTRACTS = [
    {
        "round": 18,
        "source_lane": "github_web",
        "source_name": "ECREACT / rxn-biocatalysis-tools",
        "url": "https://pypi.org/project/rxn-biocatalysis-tools/1.0.0/",
        "evidence_checked": "Reaction records are represented as rxn_smiles, ec, source and are aggregated from Rhea, BRENDA, PathBank, and MetaNetX.",
        "coverage_or_validation_signal": "All EC classes are represented; preprocessing exports train/valid/test source and target files.",
        "foodgut_training_contract": "Every imported positive needs reaction SMILES or substrate/product structures, EC/source provenance, split trace, and source database label.",
        "gap_if_missing": "A positive imported only from free-text pathway context is not comparable to ECREACT-style training data.",
    },
    {
        "round": 18,
        "source_lane": "github_web",
        "source_name": "EnzymeMap",
        "url": "https://github.com/hesther/enzymemap",
        "evidence_checked": "The project atom-maps, corrects, and suggests enzymatic reactions; its reproduction path resolves trivial names to SMILES and creates processed reaction tables.",
        "coverage_or_validation_signal": "Processed reaction tables and EC-specific processing are first-class outputs.",
        "foodgut_training_contract": "Exact names must be resolved to structures before a label is admitted; generic reaction text is rule context only.",
        "gap_if_missing": "Without name-to-structure resolution, Chinese/English/food aliases are user-facing metadata, not machine labels.",
    },
    {
        "round": 18,
        "source_lane": "github_web_and_paper",
        "source_name": "gapseq",
        "url": "https://link.springer.com/article/10.1186/s13059-021-02295-1",
        "evidence_checked": "gapseq uses curated biochemical reaction databases and validates microbial models against enzyme activity, carbon use, fermentation products, and community interactions.",
        "coverage_or_validation_signal": "Reported benchmark spans 14,931 bacterial phenotypes and a curated reaction database.",
        "foodgut_training_contract": "For gut-microbe deployment claims, reaction evidence should carry organism/sequence/phenotype or literature support when available.",
        "gap_if_missing": "A chemistry-only rule hit is not enough to claim microbial biological feasibility.",
    },
    {
        "round": 18,
        "source_lane": "database",
        "source_name": "Rhea",
        "url": "https://www.rhea-db.org/",
        "evidence_checked": "Rhea is expert-curated, uses ChEBI participants, and exposes reaction IDs, equations, status, balance, and references.",
        "coverage_or_validation_signal": "Current release 140 lists 18,343 reactions and 15,125 unique compounds.",
        "foodgut_training_contract": "Rhea exact reactions can seed candidate positives only after ChEBI structure mapping, local de-duplication, generator-route check, and leakage-safe split assignment.",
        "gap_if_missing": "Generic Rhea equations such as 'a glucuronoside' are family templates, not exact training positives.",
    },
]


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def rhea_rows() -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for path in sorted(RHEA_RAW.glob("*_rerun.json")):
        data = load_json(path)
        results = data.get("results", [])
        query_slug = path.stem.replace("_rerun", "")
        if not results:
            rows.append(
                {
                    "round": 18,
                    "source_type": "Rhea",
                    "query_slug": query_slug,
                    "source_id": "NO_RHEA_RESULT",
                    "title_or_equation": "",
                    "status": "no_result_for_query",
                    "candidate_gate": "not_a_sample",
                    "training_allowed_round18": False,
                    "why": "No result for this broad query; rerun exact substrate names or known Rhea IDs before declaring absence.",
                    "raw_file": str(path.relative_to(ROOT)),
                }
            )
            continue
        for rec in results:
            equation = rec.get("equation", "")
            status = rec.get("status", "")
            is_generic = equation.lower().startswith(("a ", "an ")) or " a " in equation.lower()
            rows.append(
                {
                    "round": 18,
                    "source_type": "Rhea",
                    "query_slug": query_slug,
                    "source_id": f"RHEA:{rec.get('id', '')}",
                    "title_or_equation": equation,
                    "status": status,
                    "candidate_gate": "rule_context_only" if is_generic else "exact_reaction_partial_pass",
                    "training_allowed_round18": False,
                    "why": (
                        "Generic participant language prevents direct positive-label admission."
                        if is_generic
                        else "Approved/specific Rhea equation still needs ChEBI structures, local de-dup, generator route, and split gate."
                    ),
                    "raw_file": str(path.relative_to(ROOT)),
                }
            )
    return rows


def pubmed_rows() -> list[dict[str, Any]]:
    summary_path = NCBI_RAW / "gut_choline_tma_cutc_esummary.json"
    if not summary_path.exists():
        return []
    data = load_json(summary_path)
    result = data.get("result", {})
    rows: list[dict[str, Any]] = []
    for uid in result.get("uids", []):
        rec = result.get(uid, {})
        rows.append(
            {
                "round": 18,
                "source_type": "PubMed",
                "query_slug": "gut_choline_tma_cutc",
                "source_id": f"PMID:{uid}",
                "title_or_equation": rec.get("title", ""),
                "status": "literature_lead",
                "candidate_gate": "manual_exact_pair_extraction_required",
                "training_allowed_round18": False,
                "why": "Literature lead must be manually extracted for substrate, product, organism/strain/community, assay condition, and analytical method.",
                "raw_file": str(summary_path.relative_to(ROOT)),
            }
        )
    return rows


def append_doc(source_contract: pd.DataFrame, live: pd.DataFrame) -> None:
    block = "\n\n## Round18 external evidence update\n\n"
    block += "External-source contract rows:\n\n"
    block += "```text\n"
    block += source_contract[["source_name", "foodgut_training_contract"]].to_string(index=False)
    block += "\n```\n\n"
    block += "Live Rhea/PubMed rerun gates:\n\n"
    block += "```text\n"
    if live.empty:
        block += "no live rerun rows\n"
    else:
        block += live[["source_type", "source_id", "candidate_gate", "training_allowed_round18"]].to_string(index=False)
        block += "\n"
    block += "```\n"
    block += "\nAdditional written artifacts:\n\n"
    block += "- `data/curation/round18_external_source_contract_update.csv`\n"
    block += "- `data/curation/round18_live_rhea_pubmed_rerun_triage.csv`\n"
    text = DOC.read_text(encoding="utf-8") if DOC.exists() else ""
    marker = "## Round18 external evidence update"
    if marker in text:
        text = text[: text.index(marker)].rstrip() + block
    else:
        text = text.rstrip() + block
    DOC.write_text(text + "\n", encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    source_contract = pd.DataFrame(EXTERNAL_SOURCE_CONTRACTS)
    live = pd.DataFrame(rhea_rows() + pubmed_rows())
    source_path = OUT / "round18_external_source_contract_update.csv"
    live_path = OUT / "round18_live_rhea_pubmed_rerun_triage.csv"
    source_contract.to_csv(source_path, index=False)
    live.to_csv(live_path, index=False)
    append_doc(source_contract, live)
    print(f"wrote {source_path.relative_to(ROOT)} rows={len(source_contract)}")
    print(f"wrote {live_path.relative_to(ROOT)} rows={len(live)}")
    print(f"updated {DOC.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
