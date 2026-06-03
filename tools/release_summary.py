from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from mlrk_prod.io_utils import atomic_write_json, atomic_write_text  # noqa: E402

SCORECARD = ROOT / "outputs" / "appraisal" / "production_scorecard.json"
DEFAULT_JSON_OUT = ROOT / "outputs" / "appraisal" / "release_summary.json"
DEFAULT_MD_OUT = ROOT / "docs" / "RELEASE_SUMMARY.md"


def read_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"status": "missing", "path": str(path.relative_to(ROOT)).replace("\\", "/")}
    return json.loads(path.read_text(encoding="utf-8"))


def build_summary(scorecard: dict[str, Any]) -> dict[str, Any]:
    internal = scorecard.get("internal_comparison", {})
    core = scorecard.get("core_panel", {})
    challenge = scorecard.get("challenge_panel", {})
    reaction = scorecard.get("reaction_family_kpis", {})
    external_results = scorecard.get("external_result_scorecard", {})
    external_review = scorecard.get("external_review_status", {})
    fresh_clone = scorecard.get("fresh_clone_report", {})
    repo_doctor = scorecard.get("repository_doctor", {})
    return {
        "status": "ready",
        "readiness_status": scorecard.get("status"),
        "production_position": "internal_research_mvp",
        "public_or_clinical_ready": False,
        "internal_metrics": {
            "r5_mean_by_method": internal.get("r5_mean_by_method", {}),
            "ablation_deltas": internal.get("ablation_deltas", {}),
            "anti_cheat_pass": internal.get("anti_cheat_pass"),
        },
        "core_panel": {
            "case_count": core.get("case_count"),
            "top5_hit_count": core.get("top5_hit_count"),
            "top5_hit_rate": core.get("top5_hit_rate"),
            "score_0_to_5": core.get("score", {}).get("score_0_to_5"),
        },
        "challenge_panel": {
            "case_count": challenge.get("case_count"),
            "top5_hit_count": challenge.get("top5_hit_count"),
            "top5_hit_rate": challenge.get("top5_hit_rate"),
            "score_0_to_5": challenge.get("score", {}).get("score_0_to_5"),
            "failure_type_counts": challenge.get("failure_type_counts", {}),
        },
        "cross_family_coverage": {
            "reaction_family_count": len(reaction.get("combined_family_kpis", [])),
            "status": reaction.get("status"),
        },
        "external_benchmark_status": {
            "result_status": external_results.get("status"),
            "tool_summaries": external_results.get("product_tool_summaries", {}),
            "claim": "external inputs and import harness are ready; real external outputs are not imported",
        },
        "external_ai_review": {
            "tool": external_review.get("tool"),
            "status": external_review.get("status"),
            "review_claim_allowed": external_review.get("review_claim_allowed", False),
        },
        "engineering_reproducibility": {
            "fresh_clone_status": fresh_clone.get("status"),
            "fresh_clone_commit": fresh_clone.get("commit"),
            "repo_doctor_status": repo_doctor.get("status"),
        },
        "remaining_blockers": scorecard.get("remaining_gap_to_full_food_microbiome_metabolite_prediction", []),
        "readiness_issues": scorecard.get("readiness_issues", []),
    }


def render_markdown(summary: dict[str, Any]) -> str:
    methods = summary["internal_metrics"]["r5_mean_by_method"]
    deltas = summary["internal_metrics"]["ablation_deltas"]
    challenge_failures = summary["challenge_panel"]["failure_type_counts"]
    blockers = "\n".join(f"- {item}" for item in summary["remaining_blockers"])
    issues = "\n".join(f"- `{item.get('code')}`: {item.get('message')}" for item in summary["readiness_issues"])
    return f"""# Release Summary

## Final Position

- Production position: `{summary['production_position']}`
- Readiness status: `{summary['readiness_status']}`
- Public or clinical ready: `{str(summary['public_or_clinical_ready']).lower()}`

The release candidate is suitable for internal research prioritization of rule-generated food-polyphenol metabolite candidates. It is not ready for public consumer, clinical, wet-lab probability, or broad gut microbiome metabolism claims.

## Internal And Ablation Metrics

- LTR_chem mean recall@5: `{methods.get('LTR_chem')}`
- Random mean recall@5: `{methods.get('random')}`
- EC-only mean recall@5: `{methods.get('ec_only')}`
- Tanimoto mean recall@5: `{methods.get('tanimoto')}`
- LTR_chem minus random recall@5: `{deltas.get('LTR_chem_minus_random_r5')}`
- LTR_chem minus EC-only recall@5: `{deltas.get('LTR_chem_minus_ec_only_r5')}`
- LTR_chem minus Tanimoto recall@5: `{deltas.get('LTR_chem_minus_tanimoto_r5')}`
- Anti-cheat pass: `{summary['internal_metrics']['anti_cheat_pass']}`

## Panel Results

| Panel | Cases | Strict top-5 hits | Strict top-5 rate | Score |
| --- | ---: | ---: | ---: | ---: |
| Core | {summary['core_panel']['case_count']} | {summary['core_panel']['top5_hit_count']} | {summary['core_panel']['top5_hit_rate']} | {summary['core_panel']['score_0_to_5']} |
| Challenge | {summary['challenge_panel']['case_count']} | {summary['challenge_panel']['top5_hit_count']} | {summary['challenge_panel']['top5_hit_rate']} | {summary['challenge_panel']['score_0_to_5']} |

Challenge failure types: `{challenge_failures}`.

## Cross-Family Coverage

- Reaction-family KPI status: `{summary['cross_family_coverage']['status']}`
- Reaction families reported: `{summary['cross_family_coverage']['reaction_family_count']}`

## External Status

- External benchmark result status: `{summary['external_benchmark_status']['result_status']}`
- External benchmark claim: {summary['external_benchmark_status']['claim']}
- External AI/code review tool: `{summary['external_ai_review']['tool']}`
- External AI/code review status: `{summary['external_ai_review']['status']}`
- External AI/code review claim allowed: `{str(summary['external_ai_review']['review_claim_allowed']).lower()}`

## Engineering Reproducibility

- Fresh-clone status: `{summary['engineering_reproducibility']['fresh_clone_status']}`
- Fresh-clone commit: `{summary['engineering_reproducibility']['fresh_clone_commit']}`
- Repository doctor status: `{summary['engineering_reproducibility']['repo_doctor_status']}`

## Remaining Blockers

{blockers}

## Readiness Issues

{issues}

## Final Claim Boundary

Use this release as an internal research MVP. Do not present it as an externally validated predictor, wet-lab occurrence model, strain-aware gut microbiome metabolism engine, clinical tool, or consumer health recommendation product.
"""


def write_outputs(summary: dict[str, Any], json_out: Path, md_out: Path) -> None:
    atomic_write_json(json_out, summary)
    atomic_write_text(md_out, render_markdown(summary))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--json-out", default=str(DEFAULT_JSON_OUT.relative_to(ROOT)).replace("\\", "/"))
    parser.add_argument("--md-out", default=str(DEFAULT_MD_OUT.relative_to(ROOT)).replace("\\", "/"))
    args = parser.parse_args()

    summary = build_summary(read_json(SCORECARD))
    json_out = ROOT / args.json_out
    md_out = ROOT / args.md_out
    write_outputs(summary, json_out, md_out)
    print(
        json.dumps(
            {
                "status": summary["status"],
                "readiness_status": summary["readiness_status"],
                "json_out": str(json_out.relative_to(ROOT)).replace("\\", "/"),
                "md_out": str(md_out.relative_to(ROOT)).replace("\\", "/"),
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
