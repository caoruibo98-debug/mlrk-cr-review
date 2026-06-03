from __future__ import annotations

import json
from pathlib import Path

from tools.release_summary import build_summary


ROOT = Path(__file__).resolve().parents[1]
SUMMARY = ROOT / "outputs" / "appraisal" / "release_summary.json"
MARKDOWN = ROOT / "docs" / "RELEASE_SUMMARY.md"
SCORECARD = ROOT / "outputs" / "appraisal" / "production_scorecard.json"


def test_release_summary_covers_required_evaluation_lanes() -> None:
    summary = json.loads(SUMMARY.read_text(encoding="utf-8"))
    assert summary["status"] == "ready"
    assert summary["production_position"] == "internal_research_mvp"
    assert summary["public_or_clinical_ready"] is False
    assert summary["internal_metrics"]["anti_cheat_pass"] is True
    assert summary["external_benchmark_status"]["result_status"] == "awaiting_external_outputs"
    assert summary["engineering_reproducibility"]["fresh_clone_status"] == "passed"


def test_release_summary_markdown_keeps_final_claim_boundary() -> None:
    text = MARKDOWN.read_text(encoding="utf-8")
    assert "Internal And Ablation Metrics" in text
    assert "Cross-Family Coverage" in text
    assert "External Status" in text
    assert "Do not present it as an externally validated predictor" in text


def test_release_summary_builder_reads_scorecard_values() -> None:
    scorecard = json.loads(SCORECARD.read_text(encoding="utf-8"))
    summary = build_summary(scorecard)
    assert summary["core_panel"]["top5_hit_rate"] == 1.0
    assert summary["challenge_panel"]["top5_hit_rate"] == 0.364
    assert summary["cross_family_coverage"]["reaction_family_count"] == 19
