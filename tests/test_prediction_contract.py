from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

from mlrk_prod.schemas import validate_prediction_payload
from mlrk_prod.biosanity import annotate_prediction_payload


ROOT = Path(__file__).resolve().parents[1]


def ensure_clean_rutin_prediction(path: Path) -> None:
    if path.exists():
        return
    cmd = [sys.executable, "-m", "mlrk_prod.cli", "predict", "--name", "rutin", "--topn", "10"]
    proc = subprocess.run(
        cmd,
        cwd=str(ROOT),
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
        check=False,
        timeout=240,
    )
    if proc.returncode != 0:
        raise RuntimeError(proc.stderr or proc.stdout)


def test_clean_rutin_prediction_payload_contract() -> None:
    path = ROOT / "outputs" / "modular" / "predictions" / "clean_rutin.json"
    ensure_clean_rutin_prediction(path)
    payload = json.loads(path.read_text(encoding="utf-8"))
    payload = annotate_prediction_payload(payload)
    issues = validate_prediction_payload(payload)
    assert issues == []
    assert payload["module"] == "B"
    assert payload["top"]
    assert payload["top"][0]["biochem_quality"]["tier"] == "high"


def test_no_candidate_prediction_writes_contract_payload() -> None:
    path = ROOT / "outputs" / "modular" / "predictions" / "clean_ellagic_acid.json"
    cmd = [sys.executable, "-m", "mlrk_prod.cli", "predict", "--name", "ellagic acid", "--topn", "10"]
    proc = subprocess.run(
        cmd,
        cwd=str(ROOT),
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
        check=False,
        timeout=240,
    )
    if proc.returncode != 0:
        raise RuntimeError(proc.stderr or proc.stdout)
    if not path.is_file():
        raise RuntimeError(f"prediction completed but output not written: {path}")
    payload = json.loads(path.read_text(encoding="utf-8"))
    payload = annotate_prediction_payload(payload)
    assert validate_prediction_payload(payload) == []
    assert payload["candidate_generation_status"] == "no_candidates"
    assert payload["n_rule_candidates"] == 0
    assert payload["top"] == []
