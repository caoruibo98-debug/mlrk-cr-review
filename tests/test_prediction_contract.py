from __future__ import annotations

import json
from pathlib import Path

from mlrk_prod.schemas import validate_prediction_payload
from mlrk_prod.biosanity import annotate_prediction_payload


def test_clean_rutin_prediction_payload_contract() -> None:
    path = Path("outputs/modular/predictions/clean_rutin.json")
    payload = json.loads(path.read_text(encoding="utf-8"))
    payload = annotate_prediction_payload(payload)
    issues = validate_prediction_payload(payload)
    assert issues == []
    assert payload["module"] == "B"
    assert payload["top"]
    assert payload["top"][0]["biochem_quality"]["tier"] == "high"
