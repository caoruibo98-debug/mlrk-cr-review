from __future__ import annotations

import csv
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from .paths import resolve_kernel_path


MODULES = ("A", "B", "C", "D")
BASELINE_METHODS = {"random", "ec_only", "tanimoto", "LTR_chem"}


@dataclass(frozen=True)
class ReadinessIssue:
    severity: str
    code: str
    message: str


def load_manifest(path: str | Path = "production_artifact_manifest.json") -> dict[str, Any]:
    manifest_path = resolve_kernel_path(path)
    with manifest_path.open("r", encoding="utf-8") as f:
        return json.load(f)


def _parse_mean(value: str) -> float:
    head = str(value).split("+/-", 1)[0].split("±", 1)[0].strip()
    return float(head)


def _read_metrics(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))


def validate_readiness(manifest: dict[str, Any] | None = None) -> list[ReadinessIssue]:
    data = manifest or load_manifest()
    issues: list[ReadinessIssue] = []

    for artifact in data.get("required_artifacts", []):
        p = resolve_kernel_path(artifact["path"])
        if not p.exists():
            issues.append(ReadinessIssue("error", "missing_artifact", f"Missing {artifact['id']}: {p}"))

    deploy_meta_path = resolve_kernel_path("outputs/modular/ltr/models_clean/deploy_meta.json")
    if deploy_meta_path.exists():
        with deploy_meta_path.open("r", encoding="utf-8") as f:
            deploy_meta = json.load(f)
        missing = [m for m in MODULES if m not in deploy_meta]
        if missing:
            issues.append(ReadinessIssue("error", "missing_module_meta", f"Missing deploy meta modules: {missing}"))
        for module, meta in deploy_meta.items():
            if int(meta.get("feat_dim", 0)) != 1536:
                issues.append(ReadinessIssue("warning", "unexpected_feature_dim", f"{module} feat_dim is {meta.get('feat_dim')}"))

    metrics_path = resolve_kernel_path("outputs/modular/ltr/clean2_metrics.csv")
    if metrics_path.exists():
        rows = _read_metrics(metrics_path)
        methods_by_module: dict[str, set[str]] = {}
        r5_by_module_method: dict[tuple[str, str], float] = {}
        for row in rows:
            module = row.get("module", "")
            method = row.get("method", "")
            methods_by_module.setdefault(module, set()).add(method)
            if "r@5" in row:
                r5_by_module_method[(module, method)] = _parse_mean(row["r@5"])
        for module in MODULES:
            missing_methods = BASELINE_METHODS - methods_by_module.get(module, set())
            if missing_methods:
                issues.append(ReadinessIssue("error", "missing_baseline", f"{module} missing methods {sorted(missing_methods)}"))
            ltr = r5_by_module_method.get((module, "LTR_chem"))
            if ltr is None:
                continue
            if ltr < float(data["readiness_gates"].get("ltr_chem_r5_minimum", 0.0)):
                issues.append(ReadinessIssue("warning", "low_ltr_r5", f"{module} LTR_chem r@5={ltr:.3f}"))
            ec = r5_by_module_method.get((module, "ec_only"))
            if ec is not None and ec >= ltr:
                issues.append(ReadinessIssue("error", "ec_probe_not_controlled", f"{module} ec_only r@5={ec:.3f} >= LTR_chem r@5={ltr:.3f}"))

    gates = data.get("readiness_gates", {})
    if not gates.get("external_head_to_head_complete", False):
        issues.append(ReadinessIssue("warning", "external_benchmark_missing", "External BioTransformer/MicrobeRX/GutBug-style benchmark is not complete."))
    if not gates.get("wet_lab_validation_complete", False):
        issues.append(ReadinessIssue("info", "wet_lab_missing", "Wet-lab validation is not complete; public biological claims must stay limited."))
    if not gates.get("public_web_app_allowed", False):
        issues.append(ReadinessIssue("info", "public_web_not_allowed", "Only internal research MVP deployment is allowed by this manifest."))

    return issues


def readiness_status(issues: list[ReadinessIssue]) -> str:
    if any(i.severity == "error" for i in issues):
        return "blocked"
    if any(i.severity == "warning" for i in issues):
        return "internal_mvp_only"
    return "ready_for_limited_beta"

