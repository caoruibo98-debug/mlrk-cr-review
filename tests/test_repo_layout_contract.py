from __future__ import annotations

import json
from pathlib import Path

from tools.repo_doctor import SCRIPT_TOOL_PAIRS, build_report


ROOT = Path(__file__).resolve().parents[1]
REPO_DOCTOR = ROOT / "outputs" / "appraisal" / "repo_doctor.json"


def test_repo_doctor_report_is_ready() -> None:
    report = json.loads(REPO_DOCTOR.read_text(encoding="utf-8"))
    assert report["status"] == "ready"
    assert report["summary"]["fail_count"] == 0


def test_repo_doctor_declares_human_reviewer_entrypoints() -> None:
    report = build_report()
    command_ids = {row["id"] for row in report["canonical_commands"]}
    for expected in {"test", "readiness", "scorecard", "family_kpis", "model_card", "repo_doctor", "predict", "serve"}:
        assert expected in command_ids


def test_scripts_and_tools_are_paired_without_moving_legacy_code() -> None:
    report = build_report()
    paired = {row["id"].split(":", 1)[1] for row in report["checks"] if row["category"] == "script_tool_pair"}
    assert set(SCRIPT_TOOL_PAIRS) <= paired
    assert (ROOT / "modular" / "predict_substrate_clean.py").is_file()
    assert (ROOT / "src").is_dir()


def test_repository_guide_documents_legacy_boundary() -> None:
    text = (ROOT / "docs" / "REPOSITORY_GUIDE.md").read_text(encoding="utf-8")
    assert "Keep legacy scientific code in place" in text
    assert "modular/" in text
    assert "scripts/repo_doctor.py" in text


def test_readme_is_reviewer_first_and_claim_limited() -> None:
    text = (ROOT / "README.md").read_text(encoding="utf-8")
    for heading in (
        "## Reviewer Summary",
        "## Fast Review Path",
        "## Current Evidence",
        "## Important Outputs",
        "## Known Limitations",
        "## Reviewer Documents",
    ):
        assert heading in text
    opening = text[:1600]
    assert "internal_mvp_only" in opening
    assert "Current safe claim" in opening
    assert "Current unsafe claim" in opening
    for command in (
        "python scripts/run_contract_tests.py",
        "python scripts/validate_production_readiness.py",
        "python scripts/repo_doctor.py",
        "python scripts/production_scorecard.py",
    ):
        assert command in text
