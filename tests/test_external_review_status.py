from __future__ import annotations

import json
from pathlib import Path

from tools.external_review_status import build_report


ROOT = Path(__file__).resolve().parents[1]
STATUS = ROOT / "outputs" / "appraisal" / "external_review_status.json"


def test_external_review_timeout_is_not_marked_as_passed() -> None:
    report = json.loads(STATUS.read_text(encoding="utf-8"))
    assert report["tool"] == "CodeRabbit"
    assert report["status"] == "timed_out"
    assert report["review_claim_allowed"] is False
    assert report["issue_count"] is None


def test_external_review_status_builder_only_allows_claim_after_completion() -> None:
    failed = build_report(
        tool="CodeRabbit",
        status="failed",
        command="coderabbit review --agent",
        duration_ms=10,
        issue_count=0,
        note="example",
    )
    assert failed["review_claim_allowed"] is False
    assert failed["issue_count"] is None

    completed = build_report(
        tool="CodeRabbit",
        status="completed",
        command="coderabbit review --agent",
        duration_ms=10,
        issue_count=0,
        note="example",
    )
    assert completed["review_claim_allowed"] is True
    assert completed["issue_count"] == 0
