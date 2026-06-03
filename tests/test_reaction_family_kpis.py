from __future__ import annotations

import csv
import json
from pathlib import Path

from tools.reaction_family_kpis import build_report


ROOT = Path(__file__).resolve().parents[1]
KPI_JSON = ROOT / "outputs" / "appraisal" / "reaction_family_kpis.json"
KPI_CSV = ROOT / "outputs" / "appraisal" / "reaction_family_kpis.csv"


def family(report: dict, panel: str, reaction_family: str) -> dict:
    for row in report["family_kpis"]:
        if row["panel"] == panel and row["reaction_family"] == reaction_family:
            return row
    raise AssertionError(f"missing {panel}/{reaction_family}")


def combined_family(report: dict, reaction_family: str) -> dict:
    for row in report["combined_family_kpis"]:
        if row["reaction_family"] == reaction_family:
            return row
    raise AssertionError(f"missing combined {reaction_family}")


def test_reaction_family_kpi_report_is_generated_and_schema_ready() -> None:
    report = json.loads(KPI_JSON.read_text(encoding="utf-8"))
    assert report["status"] == "ready"
    assert report["schema_version"] == 1
    assert report["panels"]["core"]["case_count"] == 6
    assert report["panels"]["challenge"]["case_count"] == 22
    assert report["production_priority_queue"]


def test_reaction_family_kpis_separate_generation_from_evidence_gaps() -> None:
    report = build_report()
    hydroxy = family(report, "challenge", "hydroxycinnamate_reduction")
    assert hydroxy["case_count"] == 3
    assert hydroxy["top5_hit_rate"] == 1.0
    assert hydroxy["production_gap"] == "evidence_integration_blocked"
    assert "evidence" in hydroxy["recommended_next_action"].lower()

    urolithin = family(report, "challenge", "ellagitannin_urolithin_multistep")
    assert urolithin["candidate_pool_recall"] == 0.0
    assert urolithin["production_gap"] == "candidate_generation_blocked"


def test_combined_family_kpis_preserve_core_strength_claim_boundary() -> None:
    report = build_report()
    glycoside = combined_family(report, "glycoside_hydrolysis")
    assert glycoside["case_count"] >= 6
    assert glycoside["top5_hit_count"] >= 6
    assert glycoside["model_evidence_rate"] == 0.0
    assert glycoside["production_gap"] == "evidence_integration_blocked"


def test_reaction_family_kpi_csv_has_actionable_rows() -> None:
    rows = list(csv.DictReader(KPI_CSV.open("r", encoding="utf-8", newline="")))
    assert rows
    assert any(row["panel"] == "all" for row in rows)
    assert all(row["recommended_next_action"] for row in rows)
