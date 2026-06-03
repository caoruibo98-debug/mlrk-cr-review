from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]

CORE_REPORT = ROOT / "outputs" / "appraisal" / "life_science_appraisal.json"
CHALLENGE_REPORT = ROOT / "outputs" / "appraisal" / "challenge_appraisal.json"
DEFAULT_JSON_OUT = ROOT / "outputs" / "appraisal" / "reaction_family_kpis.json"
DEFAULT_CSV_OUT = ROOT / "outputs" / "appraisal" / "reaction_family_kpis.csv"

PANEL_REPORTS = [
    ("core", "core_food_glycoside_panel", CORE_REPORT),
    ("challenge", "production_challenge_panel", CHALLENGE_REPORT),
]

CSV_FIELDS = [
    "panel",
    "reaction_family",
    "case_count",
    "expected_generated_count",
    "expected_in_candidate_pool_count",
    "candidate_pool_recall",
    "connectivity_top5_hit_count",
    "stereo_mismatch_top5_count",
    "top5_hit_count",
    "top5_hit_rate",
    "benchmark_traceable_count",
    "benchmark_traceability_rate",
    "model_evidence_count",
    "model_evidence_rate",
    "top_failure_type",
    "production_gap",
    "recommended_next_action",
]


def read_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        rel = str(path.relative_to(ROOT)).replace("\\", "/")
        raise FileNotFoundError(f"Missing appraisal report: {rel}")
    return json.loads(path.read_text(encoding="utf-8"))


def as_rank(value: Any) -> int | None:
    if value in (None, ""):
        return None
    return int(value)


def is_connectivity_top5(case: dict[str, Any]) -> bool:
    rank = as_rank(case.get("expected_rank"))
    return rank is not None and rank <= 5


def is_strict_top5(case: dict[str, Any]) -> bool:
    return is_connectivity_top5(case) and case.get("expected_full_inchikey_match") is not False


def evidence_source(case: dict[str, Any]) -> str:
    value = str(case.get("expected_evidence_source") or "none").strip()
    return value or "none"


def has_benchmark_traceability(case: dict[str, Any]) -> bool:
    return evidence_source(case) != "none"


def has_model_output_evidence(case: dict[str, Any]) -> bool:
    return evidence_source(case) not in {"none", "benchmark_panel"}


def rounded_rate(numerator: int, denominator: int) -> float:
    return round(numerator / denominator, 3) if denominator else 0.0


def production_gap(row: dict[str, Any]) -> str:
    if row["case_count"] == 0:
        return "no_cases"
    if row["candidate_pool_recall"] == 0:
        return "candidate_generation_blocked"
    if row["top5_hit_rate"] == 0 and row["expected_generated_count"] > 0:
        return "ranking_or_identity_blocked"
    if row["top5_hit_rate"] < row["candidate_pool_recall"]:
        return "ranking_or_topn_blocked"
    if row["top5_hit_rate"] >= 0.8 and row["model_evidence_rate"] == 0:
        return "evidence_integration_blocked"
    if row["top5_hit_rate"] >= 0.8:
        return "externally_unvalidated_strength"
    return "partial_coverage"


def recommended_next_action(row: dict[str, Any]) -> str:
    gap = row["production_gap"]
    if gap == "candidate_generation_blocked":
        return "Add curated reaction-family candidate generation before tuning the ranker."
    if gap == "ranking_or_identity_blocked":
        return "Inspect generated expected products for rank, identity, and stereo matching errors."
    if gap == "ranking_or_topn_blocked":
        return "Calibrate top-k ranking after confirming the expected product is in the candidate pool."
    if gap == "evidence_integration_blocked":
        return "Attach independent enzyme, microbe, strain, or literature evidence to matched products."
    if gap == "externally_unvalidated_strength":
        return "Run matched external-tool and wet-lab-facing validation before stronger production claims."
    if gap == "partial_coverage":
        return "Expand cases and split remaining failures into generation, ranking, and evidence work."
    return "Add enough benchmark cases to make this family measurable."


def failure_examples(cases: list[dict[str, Any]], limit: int = 3) -> list[dict[str, Any]]:
    failures = [case for case in cases if not is_strict_top5(case)]
    examples = []
    for case in failures[:limit]:
        examples.append(
            {
                "case_id": case.get("case_id"),
                "substrate": case.get("substrate_name"),
                "expected_product": case.get("expected_product_name"),
                "failure_type": case.get("failure_type"),
                "expected_rank": case.get("expected_rank"),
                "top_product": case.get("top_product"),
            }
        )
    return examples


def summarize_cases(panel: str, family: str, cases: list[dict[str, Any]]) -> dict[str, Any]:
    case_count = len(cases)
    expected_generated = [case for case in cases if as_rank(case.get("expected_rank")) is not None]
    expected_in_pool = [case for case in cases if case.get("expected_in_candidate_pool") is True]
    connectivity_top5 = [case for case in cases if is_connectivity_top5(case)]
    strict_top5 = [case for case in cases if is_strict_top5(case)]
    stereo_mismatch = [
        case
        for case in connectivity_top5
        if case.get("expected_full_inchikey_match") is False
    ]
    benchmark_traceable = [case for case in cases if has_benchmark_traceability(case)]
    model_evidence = [case for case in cases if has_model_output_evidence(case)]
    failures = Counter(str(case.get("failure_type") or "not_recorded") for case in cases)
    evidence = Counter(evidence_source(case) for case in cases)
    top_failure = failures.most_common(1)[0][0] if failures else "not_recorded"

    row: dict[str, Any] = {
        "panel": panel,
        "reaction_family": family,
        "case_count": case_count,
        "expected_generated_count": len(expected_generated),
        "expected_in_candidate_pool_count": len(expected_in_pool),
        "candidate_pool_recall": rounded_rate(len(expected_in_pool), case_count),
        "connectivity_top5_hit_count": len(connectivity_top5),
        "stereo_mismatch_top5_count": len(stereo_mismatch),
        "top5_hit_count": len(strict_top5),
        "top5_hit_rate": rounded_rate(len(strict_top5), case_count),
        "benchmark_traceable_count": len(benchmark_traceable),
        "benchmark_traceability_rate": rounded_rate(len(benchmark_traceable), case_count),
        "model_evidence_count": len(model_evidence),
        "model_evidence_rate": rounded_rate(len(model_evidence), case_count),
        "failure_type_counts": dict(sorted(failures.items())),
        "evidence_source_counts": dict(sorted(evidence.items())),
        "top_failure_type": top_failure,
        "example_failure_cases": failure_examples(cases),
    }
    row["production_gap"] = production_gap(row)
    row["recommended_next_action"] = recommended_next_action(row)
    return row


def group_by_family(panel: str, cases: list[dict[str, Any]]) -> list[dict[str, Any]]:
    families: dict[str, list[dict[str, Any]]] = {}
    for case in cases:
        family = str(case.get("reaction_family") or "unknown")
        families.setdefault(family, []).append(case)
    return [
        summarize_cases(panel, family, families[family])
        for family in sorted(families)
    ]


def build_report() -> dict[str, Any]:
    panel_reports: dict[str, dict[str, Any]] = {}
    all_cases: list[dict[str, Any]] = []
    family_rows: list[dict[str, Any]] = []

    for panel, label, path in PANEL_REPORTS:
        appraisal = read_json(path)
        cases = list(appraisal.get("cases", []))
        all_cases.extend(cases)
        rows = group_by_family(panel, cases)
        family_rows.extend(rows)
        panel_reports[panel] = {
            "label": label,
            "panel_file": appraisal.get("panel_file"),
            "status": appraisal.get("status"),
            "score": appraisal.get("score", {}),
            "case_count": len(cases),
            "family_count": len(rows),
            "top5_hit_rate": rounded_rate(sum(1 for case in cases if is_strict_top5(case)), len(cases)),
            "candidate_pool_recall": rounded_rate(
                sum(1 for case in cases if case.get("expected_in_candidate_pool") is True),
                len(cases),
            ),
        }

    combined_rows = group_by_family("all", all_cases)
    production_priorities = sorted(
        combined_rows,
        key=lambda row: (
            row["top5_hit_rate"],
            row["candidate_pool_recall"],
            -row["case_count"],
            row["reaction_family"],
        ),
    )

    return {
        "status": "ready",
        "schema_version": 1,
        "source_reports": {
            "core": str(CORE_REPORT.relative_to(ROOT)).replace("\\", "/"),
            "challenge": str(CHALLENGE_REPORT.relative_to(ROOT)).replace("\\", "/"),
        },
        "metric_definitions": {
            "top5_hit_rate": "Strict expected-product top-5 recall; connectivity hits with full-InChIKey mismatches are not counted.",
            "candidate_pool_recall": "Fraction of cases where the expected product appears anywhere in the generated candidate pool.",
            "benchmark_traceability_rate": "Fraction of cases where the expected product has any traceable benchmark or model evidence source.",
            "model_evidence_rate": "Fraction of cases where matched output evidence is not limited to benchmark-panel labels.",
        },
        "panels": panel_reports,
        "family_kpis": family_rows,
        "combined_family_kpis": combined_rows,
        "production_priority_queue": production_priorities,
        "interpretation": {
            "main_use": "Locate production blockers by reaction family before adding rules, labels, or ranker features.",
            "claim_boundary": "Family KPI strength is internal benchmark evidence only; it is not external validation or wet-lab proof.",
        },
    }


def csv_row(row: dict[str, Any]) -> dict[str, Any]:
    output = {field: row.get(field) for field in CSV_FIELDS}
    return output


def write_outputs(report: dict[str, Any], json_out: Path, csv_out: Path) -> None:
    json_out.parent.mkdir(parents=True, exist_ok=True)
    json_out.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    csv_out.parent.mkdir(parents=True, exist_ok=True)
    rows = [*report["family_kpis"], *report["combined_family_kpis"]]
    with csv_out.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=CSV_FIELDS)
        writer.writeheader()
        for row in rows:
            writer.writerow(csv_row(row))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--json-out", default=str(DEFAULT_JSON_OUT.relative_to(ROOT)).replace("\\", "/"))
    parser.add_argument("--csv-out", default=str(DEFAULT_CSV_OUT.relative_to(ROOT)).replace("\\", "/"))
    args = parser.parse_args()

    report = build_report()
    json_out = ROOT / args.json_out
    csv_out = ROOT / args.csv_out
    write_outputs(report, json_out, csv_out)
    print(
        json.dumps(
            {
                "status": report["status"],
                "json_out": str(json_out.relative_to(ROOT)).replace("\\", "/"),
                "csv_out": str(csv_out.relative_to(ROOT)).replace("\\", "/"),
                "families": len(report["combined_family_kpis"]),
                "highest_priority_gap": report["production_priority_queue"][0]["production_gap"],
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
