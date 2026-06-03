from __future__ import annotations

import json
from pathlib import Path

from tools.fresh_clone_report import build_report


ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "outputs" / "appraisal" / "fresh_clone_report.json"


def test_fresh_clone_report_records_release_candidate_pass() -> None:
    report = json.loads(REPORT.read_text(encoding="utf-8"))
    assert report["status"] == "passed"
    assert report["source_ref"] == "generation_18"
    assert report["commit"].startswith("fda7f6c")
    assert report["metrics"]["core_top5"] == 1.0
    assert report["metrics"]["challenge_top5"] == 0.364
    assert all(check["status"] == "passed" for check in report["checks"])


def test_fresh_clone_report_builder_keeps_claim_boundary() -> None:
    report = build_report(
        source_url="https://example.invalid/repo.git",
        source_ref="generation_x",
        commit="abc123",
        clone_path="D:/tmp/clone",
        contract_tests="PASSED 1 contract tests",
        scorecard_status="internal_mvp_only",
        core_top5=1.0,
        challenge_top5=0.1,
        readiness_status="internal_mvp_only",
        repo_doctor_status="ready",
        repo_doctor_checks=1,
        model_card_status="ready",
    )
    assert report["status"] == "passed"
    assert "not external biological validation" in report["claim_boundary"]
