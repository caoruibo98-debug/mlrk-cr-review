from __future__ import annotations

import json
from pathlib import Path

from tools.model_card import build_summary, read_json


ROOT = Path(__file__).resolve().parents[1]
SUMMARY = ROOT / "outputs" / "appraisal" / "model_card_summary.json"
MODEL_CARD = ROOT / "docs" / "MODEL_CARD.md"
SCORECARD = ROOT / "outputs" / "appraisal" / "production_scorecard.json"


def test_model_card_summary_is_scorecard_backed() -> None:
    summary = json.loads(SUMMARY.read_text(encoding="utf-8"))
    assert summary["status"] == "ready"
    assert summary["source_scorecard"] == "outputs/appraisal/production_scorecard.json"
    assert summary["readiness_status"] == "internal_mvp_only"
    assert summary["challenge_panel"]["top5_hit_count"] == 8
    assert summary["challenge_panel"]["case_count"] == 22


def test_model_card_markdown_keeps_claim_boundary_visible() -> None:
    text = MODEL_CARD.read_text(encoding="utf-8")
    assert "Generated from `outputs/appraisal/production_scorecard.json`" in text
    assert "not wet-lab probabilities" in text
    assert "Production challenge panel | 22 | 8 | 0.364" in text
    assert "External result status: `awaiting_external_outputs`" in text


def test_model_card_builder_matches_current_scorecard() -> None:
    scorecard = read_json(SCORECARD)
    summary = build_summary(scorecard)
    assert summary["core_panel"]["top5_hit_rate"] == 1.0
    assert summary["external_result_status"] == "awaiting_external_outputs"
    assert summary["external_ai_review_status"] == "timed_out"
    assert summary["external_ai_review_claim_allowed"] is False
    assert summary["generalization_evidence_level"]["level"] == 3
