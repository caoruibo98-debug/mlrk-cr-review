from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from mlrk_prod.io_utils import atomic_write_json  # noqa: E402

DEFAULT_OUT = ROOT / "outputs" / "appraisal" / "repo_doctor.json"

SCRIPT_TOOL_PAIRS = [
    "export_external_benchmarks.py",
    "external_review_status.py",
    "freeze_generation.py",
    "fresh_clone_report.py",
    "life_science_appraisal.py",
    "model_card.py",
    "production_scorecard.py",
    "reaction_family_kpis.py",
    "release_summary.py",
    "repo_doctor.py",
    "run_contract_tests.py",
    "score_external_results.py",
    "validate_production_readiness.py",
]

REQUIRED_PATHS = {
    "root_docs": [
        "README.md",
        "README_PRODUCTION_CANDIDATE.md",
        "production_artifact_manifest.json",
        "pyproject.toml",
        "requirements.txt",
    ],
    "production_docs": [
        "docs/API_CONTRACT.md",
        "docs/EVALUATION_PROTOCOL.md",
        "docs/MODEL_CARD.md",
        "docs/PRODUCTION_READINESS_REVIEW.md",
        "docs/REPOSITORY_GUIDE.md",
        "docs/RELEASE_SUMMARY.md",
        "docs/SECURITY_AND_DEPLOYMENT_BOUNDARIES.md",
        "docs/WEB_APP_MVP_SPEC.md",
        "docs/production/ITERATION_LEDGER.md",
        "docs/production/SCIENTIFIC_POSITIONING.md",
    ],
    "legacy_boundaries": [
        "modular/predict_substrate_clean.py",
        "modular/ltr_deploy.py",
        "src",
    ],
    "versioned_outputs": [
        "outputs/appraisal/life_science_appraisal.json",
        "outputs/appraisal/challenge_appraisal.json",
        "outputs/appraisal/production_scorecard.json",
        "outputs/appraisal/reaction_family_kpis.json",
        "outputs/appraisal/release_summary.json",
        "outputs/appraisal/model_card_summary.json",
        "outputs/appraisal/external_review_status.json",
        "outputs/appraisal/fresh_clone_report.json",
        "outputs/external_benchmarks/manifest.json",
        "outputs/external_benchmarks/external_result_scorecard.json",
        "outputs/modular/ltr/clean2_metrics.csv",
        "outputs/modular/ltr/models_clean/deploy_meta.json",
    ],
}

CANONICAL_COMMANDS = [
    {"id": "test", "command": "python scripts/run_contract_tests.py"},
    {"id": "readiness", "command": "python scripts/validate_production_readiness.py"},
    {"id": "scorecard", "command": "python scripts/production_scorecard.py"},
    {"id": "family_kpis", "command": "python scripts/reaction_family_kpis.py"},
    {"id": "model_card", "command": "python scripts/model_card.py"},
    {"id": "release_summary", "command": "python scripts/release_summary.py"},
    {"id": "external_export", "command": "python scripts/export_external_benchmarks.py"},
    {"id": "external_score", "command": "python scripts/score_external_results.py"},
    {"id": "external_review_status", "command": "python scripts/external_review_status.py"},
    {"id": "fresh_clone_report", "command": "python scripts/fresh_clone_report.py"},
    {"id": "repo_doctor", "command": "python scripts/repo_doctor.py"},
    {"id": "predict", "command": "python -m mlrk_prod.cli predict --name rutin --topn 10"},
    {"id": "serve", "command": "python -m mlrk_prod.cli serve --host 127.0.0.1 --port 8765"},
]


def rel(path: Path) -> str:
    return str(path.relative_to(ROOT)).replace("\\", "/")


def add_path_check(checks: list[dict[str, Any]], category: str, path: str) -> None:
    full = ROOT / path
    checks.append(
        {
            "id": f"{category}:{path}",
            "category": category,
            "path": path,
            "status": "pass" if full.exists() else "fail",
            "message": "exists" if full.exists() else "missing required repository path",
        }
    )


def wrapper_mentions_tool(script: Path, tool_name: str) -> bool:
    if not script.exists():
        return False
    text = script.read_text(encoding="utf-8")
    return "runpy.run_path" in text and tool_name in text


def add_script_pair_check(checks: list[dict[str, Any]], name: str) -> None:
    script = ROOT / "scripts" / name
    tool = ROOT / "tools" / name
    status = "pass" if script.exists() and tool.exists() and wrapper_mentions_tool(script, name) else "fail"
    checks.append(
        {
            "id": f"script_tool_pair:{name}",
            "category": "script_tool_pair",
            "path": f"scripts/{name}",
            "tool_path": f"tools/{name}",
            "status": status,
            "message": "script wrapper and tool implementation are paired" if status == "pass" else "missing or unpaired script/tool entrypoint",
        }
    )


def build_report() -> dict[str, Any]:
    checks: list[dict[str, Any]] = []
    for category, paths in REQUIRED_PATHS.items():
        for path in paths:
            add_path_check(checks, category, path)
    for name in SCRIPT_TOOL_PAIRS:
        add_script_pair_check(checks, name)

    failures = [check for check in checks if check["status"] != "pass"]
    return {
        "status": "ready" if not failures else "needs_attention",
        "schema_version": 1,
        "summary": {
            "check_count": len(checks),
            "pass_count": len(checks) - len(failures),
            "fail_count": len(failures),
        },
        "canonical_commands": CANONICAL_COMMANDS,
        "legacy_policy": {
            "status": "preserved",
            "message": "Production wrappers live beside legacy code; legacy modular and src paths are not moved during standardization.",
        },
        "checks": checks,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", default=str(DEFAULT_OUT.relative_to(ROOT)).replace("\\", "/"))
    args = parser.parse_args()

    report = build_report()
    out = ROOT / args.out
    atomic_write_json(out, report)
    print(
        json.dumps(
            {
                "status": report["status"],
                "out": rel(out),
                "checks": report["summary"]["check_count"],
                "failures": report["summary"]["fail_count"],
            },
            indent=2,
        )
    )
    return 0 if report["status"] == "ready" else 1


if __name__ == "__main__":
    raise SystemExit(main())
