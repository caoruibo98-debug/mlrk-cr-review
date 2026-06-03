from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "src"))

import kio  # noqa: E402


BENCHMARK_DIR = ROOT / "outputs" / "external_benchmarks"
EXPECTED = BENCHMARK_DIR / "benchmark_panel_unified.csv"
DEFAULT_PRODUCT_RESULTS = BENCHMARK_DIR / "external_product_results.csv"
DEFAULT_ENZYME_RESULTS = BENCHMARK_DIR / "external_enzyme_results.csv"
DEFAULT_OUT = BENCHMARK_DIR / "external_result_scorecard.json"
TEMPLATES_DIR = BENCHMARK_DIR / "result_templates"


PRODUCT_COLUMNS = [
    "tool",
    "case_id",
    "rank",
    "predicted_product_name",
    "predicted_product_smiles",
    "predicted_product_inchikey",
    "reaction_id",
    "enzyme_ec",
    "evidence_reference",
    "raw_source_file",
]

ENZYME_COLUMNS = [
    "tool",
    "case_id",
    "substrate_pubchem_cid",
    "predicted_ec",
    "predicted_enzyme",
    "predicted_microbe_or_strain",
    "score",
    "raw_source_file",
]


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))


def write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def normalize_inchikey(value: str | None) -> str:
    return "" if not value else str(value).strip().upper()


def inchikey_from_smiles(smiles: str | None) -> str:
    if not smiles:
        return ""
    try:
        return normalize_inchikey(kio.smiles_to_inchikey(smiles))
    except Exception:
        return ""


def block1(value: str | None) -> str:
    value = normalize_inchikey(value)
    return kio.inchikey_block1(value) if value else ""


def rank_value(row: dict[str, str], fallback: int) -> int:
    raw = row.get("rank") or row.get("Rank") or row.get("prediction_rank") or ""
    try:
        return int(float(raw))
    except (TypeError, ValueError):
        return fallback


def row_tool(row: dict[str, str], fallback: str) -> str:
    return (row.get("tool") or row.get("Tool") or fallback).strip() or fallback


def row_case_id(row: dict[str, str]) -> str:
    return (row.get("case_id") or row.get("Case ID") or row.get("id") or row.get("query_id") or "").strip()


def row_product_smiles(row: dict[str, str]) -> str:
    candidates = [
        "predicted_product_smiles",
        "product_smiles",
        "Product SMILES",
        "product",
        "SMILES",
        "smiles",
    ]
    for col in candidates:
        value = row.get(col)
        if value:
            return value.strip()
    return ""


def row_product_inchikey(row: dict[str, str]) -> str:
    candidates = [
        "predicted_product_inchikey",
        "product_inchikey",
        "Product InChIKey",
        "inchikey",
        "InChIKey",
    ]
    for col in candidates:
        value = normalize_inchikey(row.get(col))
        if value:
            return value
    return inchikey_from_smiles(row_product_smiles(row))


def expected_map() -> dict[str, dict[str, str]]:
    rows = read_csv(EXPECTED)
    return {row["case_id"]: row for row in rows}


def load_product_results(paths: list[Path]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for path in paths:
        if not path.exists():
            continue
        source_rows = read_csv(path)
        fallback_tool = path.stem.replace("_results", "")
        for idx, row in enumerate(source_rows, start=1):
            case_id = row_case_id(row)
            if not case_id:
                continue
            product_inchikey = row_product_inchikey(row)
            rows.append(
                {
                    "tool": row_tool(row, fallback_tool),
                    "case_id": case_id,
                    "rank": str(rank_value(row, idx)),
                    "predicted_product_name": row.get("predicted_product_name") or row.get("product_name") or row.get("Product Name") or "",
                    "predicted_product_smiles": row_product_smiles(row),
                    "predicted_product_inchikey": product_inchikey,
                    "predicted_product_block1": block1(product_inchikey),
                    "reaction_id": row.get("reaction_id") or row.get("Reaction ID") or row.get("reaction") or "",
                    "enzyme_ec": row.get("enzyme_ec") or row.get("EC") or row.get("ec") or "",
                    "evidence_reference": row.get("evidence_reference") or row.get("PubMed") or row.get("pubmed") or row.get("reference") or "",
                    "raw_source_file": str(path.relative_to(ROOT)).replace("\\", "/") if path.is_absolute() else str(path),
                }
            )
    return rows


def score_product_results(product_rows: list[dict[str, str]], expected: dict[str, dict[str, str]]) -> dict[str, Any]:
    grouped: dict[str, dict[str, list[dict[str, str]]]] = defaultdict(lambda: defaultdict(list))
    for row in product_rows:
        if row["case_id"] in expected:
            grouped[row["tool"]][row["case_id"]].append(row)

    summaries: dict[str, Any] = {}
    case_results: list[dict[str, Any]] = []
    for tool, by_case in sorted(grouped.items()):
        top1 = top5 = any_hit = full = 0
        for case_id, exp in expected.items():
            rows = sorted(by_case.get(case_id, []), key=lambda r: int(r["rank"]))
            exp_block = exp["expected_product_block1"]
            exp_full = normalize_inchikey(exp["expected_product_inchikey"])
            hit_rank = None
            full_hit_rank = None
            for row in rows:
                predicted = normalize_inchikey(row.get("predicted_product_inchikey"))
                if not predicted:
                    continue
                if predicted == exp_full and full_hit_rank is None:
                    full_hit_rank = int(row["rank"])
                if block1(predicted) == exp_block and hit_rank is None:
                    hit_rank = int(row["rank"])
            if hit_rank == 1:
                top1 += 1
            if hit_rank is not None and hit_rank <= 5:
                top5 += 1
            if hit_rank is not None:
                any_hit += 1
            if full_hit_rank is not None:
                full += 1
            case_results.append(
                {
                    "tool": tool,
                    "case_id": case_id,
                    "panel": exp["panel"],
                    "expected_product_name": exp["expected_product_name"],
                    "expected_product_block1": exp_block,
                    "expected_product_inchikey": exp_full,
                    "predicted_count": len(rows),
                    "block1_hit_rank": hit_rank,
                    "full_inchikey_hit_rank": full_hit_rank,
                    "top5_block1_hit": hit_rank is not None and hit_rank <= 5,
                }
            )
        denom = max(1, len(expected))
        summaries[tool] = {
            "case_count": len(expected),
            "cases_with_predictions": len(by_case),
            "top1_block1_recall": round(top1 / denom, 3),
            "top5_block1_recall": round(top5 / denom, 3),
            "any_rank_block1_recall": round(any_hit / denom, 3),
            "any_rank_full_inchikey_recall": round(full / denom, 3),
        }
    return {"tool_summaries": summaries, "case_results": case_results}


def load_enzyme_results(paths: list[Path]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for path in paths:
        if not path.exists():
            continue
        fallback_tool = path.stem.replace("_results", "")
        for row in read_csv(path):
            case_id = row_case_id(row)
            if not case_id:
                continue
            rows.append(
                {
                    "tool": row_tool(row, fallback_tool),
                    "case_id": case_id,
                    "substrate_pubchem_cid": row.get("substrate_pubchem_cid") or row.get("PubChem CID") or row.get("cid") or "",
                    "predicted_ec": row.get("predicted_ec") or row.get("EC") or row.get("ec") or "",
                    "predicted_enzyme": row.get("predicted_enzyme") or row.get("enzyme") or row.get("enzyme_name") or "",
                    "predicted_microbe_or_strain": row.get("predicted_microbe_or_strain") or row.get("microbe") or row.get("strain") or "",
                    "score": row.get("score") or row.get("probability") or "",
                    "raw_source_file": str(path.relative_to(ROOT)).replace("\\", "/") if path.is_absolute() else str(path),
                }
            )
    return rows


def summarize_enzyme_results(rows: list[dict[str, str]], expected: dict[str, dict[str, str]]) -> dict[str, Any]:
    by_tool: dict[str, set[str]] = defaultdict(set)
    for row in rows:
        if row["case_id"] in expected:
            by_tool[row["tool"]].add(row["case_id"])
    return {
        tool: {
            "case_count": len(expected),
            "cases_with_enzyme_or_ec_predictions": len(cases),
            "scoring_note": "No expected EC labels are curated yet; these rows support evidence plausibility review, not product recall.",
        }
        for tool, cases in sorted(by_tool.items())
    }


def write_templates(expected: dict[str, dict[str, str]]) -> None:
    sample_rows = []
    for case_id, row in list(expected.items())[:3]:
        sample_rows.append(
            {
                "tool": "example_product_tool",
                "case_id": case_id,
                "rank": 1,
                "predicted_product_name": row["expected_product_name"],
                "predicted_product_smiles": row["expected_product_smiles"],
                "predicted_product_inchikey": row["expected_product_inchikey"],
                "reaction_id": "",
                "enzyme_ec": "",
                "evidence_reference": row.get("reference", ""),
                "raw_source_file": "",
            }
        )
    write_csv(TEMPLATES_DIR / "external_product_results_template.csv", sample_rows, PRODUCT_COLUMNS)
    enzyme_rows = [
        {
            "tool": "example_ec_tool",
            "case_id": case_id,
            "substrate_pubchem_cid": row.get("substrate_pubchem_cid", ""),
            "predicted_ec": "",
            "predicted_enzyme": "",
            "predicted_microbe_or_strain": "",
            "score": "",
            "raw_source_file": "",
        }
        for case_id, row in list(expected.items())[:3]
    ]
    write_csv(TEMPLATES_DIR / "external_enzyme_results_template.csv", enzyme_rows, ENZYME_COLUMNS)


def resolve_paths(values: list[str] | None, default: Path) -> list[Path]:
    if values:
        return [Path(v) if Path(v).is_absolute() else ROOT / v for v in values]
    return [default]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--product-results", action="append", default=None)
    parser.add_argument("--enzyme-results", action="append", default=None)
    parser.add_argument("--out", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    expected = expected_map()
    write_templates(expected)
    product_paths = resolve_paths(args.product_results, DEFAULT_PRODUCT_RESULTS)
    enzyme_paths = resolve_paths(args.enzyme_results, DEFAULT_ENZYME_RESULTS)
    product_rows = load_product_results(product_paths)
    enzyme_rows = load_enzyme_results(enzyme_paths)
    product_scores = score_product_results(product_rows, expected)
    enzyme_summary = summarize_enzyme_results(enzyme_rows, expected)

    status = "scored_external_outputs" if product_rows or enzyme_rows else "awaiting_external_outputs"
    report = {
        "status": status,
        "case_count": len(expected),
        "product_result_files": [str(p.relative_to(ROOT)).replace("\\", "/") for p in product_paths if p.exists()],
        "enzyme_result_files": [str(p.relative_to(ROOT)).replace("\\", "/") for p in enzyme_paths if p.exists()],
        "template_files": [
            str((TEMPLATES_DIR / "external_product_results_template.csv").relative_to(ROOT)).replace("\\", "/"),
            str((TEMPLATES_DIR / "external_enzyme_results_template.csv").relative_to(ROOT)).replace("\\", "/"),
        ],
        "product_tool_summaries": product_scores["tool_summaries"],
        "product_case_results": product_scores["case_results"],
        "enzyme_tool_summaries": enzyme_summary,
        "claim_boundary": "External validation is not complete unless real external tool outputs are imported and scored.",
    }
    out = Path(args.out)
    if not out.is_absolute():
        out = ROOT / out
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps({"status": status, "out": str(out), "product_tools": sorted(report["product_tool_summaries"])}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
