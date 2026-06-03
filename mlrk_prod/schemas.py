from __future__ import annotations

from dataclasses import dataclass
from typing import Any


@dataclass(frozen=True)
class PredictionContractIssue:
    field: str
    message: str


def validate_prediction_payload(payload: dict[str, Any]) -> list[PredictionContractIssue]:
    issues: list[PredictionContractIssue] = []
    for field in ("input", "module", "module_name", "n_rule_candidates", "honest_note", "top"):
        if field not in payload:
            issues.append(PredictionContractIssue(field, "Missing required prediction field."))
    top = payload.get("top")
    if top is not None and not isinstance(top, list):
        issues.append(PredictionContractIssue("top", "Expected a list of ranked candidates."))
    if isinstance(top, list):
        for idx, row in enumerate(top):
            for field in ("rank", "product_name", "product_smiles", "score", "evidence"):
                if field not in row:
                    issues.append(PredictionContractIssue(f"top[{idx}].{field}", "Missing candidate field."))
    return issues


def public_claims_allowed(payload: dict[str, Any]) -> bool:
    note = str(payload.get("honest_note", "")).lower()
    return "non-wet-lab" in note or "non wet-lab" in note or "非湿实验概率" in note

