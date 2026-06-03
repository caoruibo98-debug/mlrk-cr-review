from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
SCORECARD = ROOT / "outputs" / "appraisal" / "production_scorecard.json"
DEFAULT_SUMMARY_OUT = ROOT / "outputs" / "appraisal" / "model_card_summary.json"
DEFAULT_MODEL_CARD_OUT = ROOT / "docs" / "MODEL_CARD.md"


def read_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        rel = str(path.relative_to(ROOT)).replace("\\", "/")
        raise FileNotFoundError(f"Missing required model-card source: {rel}")
    return json.loads(path.read_text(encoding="utf-8"))


def panel_summary(panel: dict[str, Any]) -> dict[str, Any]:
    return {
        "case_count": int(panel.get("case_count", 0)),
        "top5_hit_count": int(panel.get("top5_hit_count", 0)),
        "top5_hit_rate": float(panel.get("top5_hit_rate", 0.0)),
        "score_0_to_5": panel.get("score", {}).get("score_0_to_5"),
        "failure_type_counts": panel.get("failure_type_counts", {}),
        "evidence_source_counts": panel.get("evidence_source_counts", {}),
    }


def build_summary(scorecard: dict[str, Any]) -> dict[str, Any]:
    reaction_family_kpis = scorecard.get("reaction_family_kpis", {})
    combined_families = reaction_family_kpis.get("combined_family_kpis", [])
    external_results = scorecard.get("external_result_scorecard", {})
    external_review = scorecard.get("external_review_status", {})
    repo_doctor = scorecard.get("repository_doctor", {})
    fresh_clone = scorecard.get("fresh_clone_report", {})
    return {
        "status": "ready",
        "source_scorecard": "outputs/appraisal/production_scorecard.json",
        "readiness_status": scorecard.get("status"),
        "system_type": [
            "rule_based_candidate_generation",
            "chemistry_only_learned_to_rank_scoring",
            "evidence_overlay",
            "contract_checked_internal_api",
        ],
        "generalization_evidence_level": {
            "level": 3,
            "label": "internal benchmark evidence only; no external tool outputs or wet-lab validation imported",
        },
        "allowed_claim": scorecard.get("current_claim"),
        "not_allowed_claims": [
            "Broad prediction of all food-derived gut microbial metabolites.",
            "Strain-aware or genome-aware gut microbiome metabolism prediction.",
            "Wet-lab occurrence probability.",
            "Consumer health, clinical, or treatment recommendation.",
        ],
        "core_panel": panel_summary(scorecard.get("core_panel", {})),
        "challenge_panel": panel_summary(scorecard.get("challenge_panel", {})),
        "internal_comparison": scorecard.get("internal_comparison", {}),
        "reaction_family_count": len(combined_families),
        "external_result_status": external_results.get("status", "unknown"),
        "external_ai_review_status": external_review.get("status", "unknown"),
        "external_ai_review_claim_allowed": external_review.get("review_claim_allowed", False),
        "repository_doctor_status": repo_doctor.get("status", "unknown"),
        "fresh_clone_status": fresh_clone.get("status", "unknown"),
        "readiness_issues": scorecard.get("readiness_issues", []),
        "remaining_gaps": scorecard.get("remaining_gap_to_full_food_microbiome_metabolite_prediction", []),
    }


def pct(rate: float) -> str:
    return f"{rate:.3f}"


def render_model_card(summary: dict[str, Any]) -> str:
    core = summary["core_panel"]
    challenge = summary["challenge_panel"]
    internal = summary.get("internal_comparison", {})
    method_means = internal.get("r5_mean_by_method", {})
    issues = summary.get("readiness_issues", [])
    gaps = summary.get("remaining_gaps", [])
    review_claim_allowed = str(summary["external_ai_review_claim_allowed"]).lower()
    not_allowed = "\n".join(f"- {item}" for item in summary["not_allowed_claims"])
    issue_lines = "\n".join(f"- `{item.get('code')}`: {item.get('message')}" for item in issues)
    gap_lines = "\n".join(f"- {item}" for item in gaps)

    return f"""# Model Card

Generated from `{summary['source_scorecard']}`.

## Status

- Readiness status: `{summary['readiness_status']}`
- Generalization evidence level: `{summary['generalization_evidence_level']['level']}` ({summary['generalization_evidence_level']['label']})
- External result status: `{summary['external_result_status']}`
- External AI/code review status: `{summary['external_ai_review_status']}` (claim allowed: `{review_claim_allowed}`)
- Repository doctor status: `{summary['repository_doctor_status']}`
- Fresh-clone release-candidate status: `{summary['fresh_clone_status']}`

## System Type

This is a hybrid internal research system:

- rule-based candidate generation,
- chemistry-only learned-to-rank scoring,
- evidence overlay,
- contract-checked internal API.

It is not an externally validated biological occurrence predictor.

## Intended Use

{summary['allowed_claim']}

## Not Intended Use

{not_allowed}

## Evaluation Snapshot

| Evaluation | Cases | Strict top-5 hits | Strict top-5 rate | Score |
| --- | ---: | ---: | ---: | ---: |
| Core food-glycoside panel | {core['case_count']} | {core['top5_hit_count']} | {pct(core['top5_hit_rate'])} | {core['score_0_to_5']} |
| Production challenge panel | {challenge['case_count']} | {challenge['top5_hit_count']} | {pct(challenge['top5_hit_rate'])} | {challenge['score_0_to_5']} |

Internal LTR_chem mean recall@5 is `{method_means.get('LTR_chem')}` versus random `{method_means.get('random')}`, EC-only `{method_means.get('ec_only')}`, and Tanimoto `{method_means.get('tanimoto')}`.

Reaction-family KPI count: `{summary['reaction_family_count']}`.

## Readiness Issues

{issue_lines}

## Remaining Gaps

{gap_lines}

## Claim Boundary

Scores are ranking and prioritization signals over generated candidates. They are not wet-lab probabilities, occurrence probabilities, clinical signals, or consumer health recommendations.
"""


def write_outputs(summary: dict[str, Any], summary_out: Path, model_card_out: Path) -> None:
    summary_out.parent.mkdir(parents=True, exist_ok=True)
    summary_out.write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    model_card_out.parent.mkdir(parents=True, exist_ok=True)
    model_card_out.write_text(render_model_card(summary), encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--summary-out", default=str(DEFAULT_SUMMARY_OUT.relative_to(ROOT)).replace("\\", "/"))
    parser.add_argument("--model-card-out", default=str(DEFAULT_MODEL_CARD_OUT.relative_to(ROOT)).replace("\\", "/"))
    args = parser.parse_args()

    scorecard = read_json(SCORECARD)
    summary = build_summary(scorecard)
    summary_out = ROOT / args.summary_out
    model_card_out = ROOT / args.model_card_out
    write_outputs(summary, summary_out, model_card_out)
    print(
        json.dumps(
            {
                "status": summary["status"],
                "readiness_status": summary["readiness_status"],
                "summary_out": str(summary_out.relative_to(ROOT)).replace("\\", "/"),
                "model_card_out": str(model_card_out.relative_to(ROOT)).replace("\\", "/"),
                "challenge_top5": summary["challenge_panel"]["top5_hit_rate"],
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
