from __future__ import annotations

from mlrk_prod.manifest import readiness_status, validate_readiness


def test_artifact_manifest_is_not_blocked_by_missing_files() -> None:
    issues = validate_readiness()
    missing = [i for i in issues if i.code == "missing_artifact"]
    assert missing == []
    assert readiness_status(issues) in {"internal_mvp_only", "ready_for_limited_beta"}


def test_external_validation_gate_is_explicit() -> None:
    issues = validate_readiness()
    assert any(i.code == "external_benchmark_results_missing" for i in issues)
