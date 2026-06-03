from __future__ import annotations

import json
import os
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path

from .biosanity import annotate_prediction_payload
from .paths import kernel_root
from .schemas import validate_prediction_payload


MAX_NAME_LENGTH = 120
MAX_SMILES_LENGTH = 1500
DEFAULT_TIMEOUT_SECONDS = 180


@dataclass(frozen=True)
class PredictionRequest:
    name: str | None = None
    smiles: str | None = None
    topn: int = 12


@dataclass(frozen=True)
class PredictionRun:
    payload: dict
    stdout: str
    stderr: str
    output_path: Path
    elapsed_seconds: float


def _safe_label(label: str) -> str:
    return "".join(c if c.isalnum() or c in "-_" else "_" for c in label)


def validate_request(req: PredictionRequest) -> None:
    if bool(req.name) == bool(req.smiles):
        raise ValueError("Provide exactly one of name or smiles.")
    if req.name and len(req.name) > MAX_NAME_LENGTH:
        raise ValueError(f"name is too long; max {MAX_NAME_LENGTH} characters.")
    if req.smiles and len(req.smiles) > MAX_SMILES_LENGTH:
        raise ValueError(f"smiles is too long; max {MAX_SMILES_LENGTH} characters.")
    if not 1 <= int(req.topn) <= 50:
        raise ValueError("topn must be between 1 and 50.")


def run_prediction(req: PredictionRequest, timeout_seconds: int = DEFAULT_TIMEOUT_SECONDS) -> PredictionRun:
    validate_request(req)
    root = kernel_root()
    script = root / "modular" / "predict_substrate_clean.py"
    if not script.exists():
        raise FileNotFoundError(script)

    label = req.name if req.name else "query"
    expected_output = root / "outputs" / "modular" / "predictions" / f"clean_{_safe_label(label)}.json"
    start = time.time()
    before_mtime = expected_output.stat().st_mtime if expected_output.exists() else None

    cmd = [sys.executable, str(script), "--topn", str(int(req.topn))]
    if req.name:
        cmd.extend(["--name", req.name])
    else:
        cmd.extend(["--smiles", req.smiles or ""])

    env = os.environ.copy()
    env["PYTHONIOENCODING"] = "utf-8"
    proc = subprocess.run(
        cmd,
        cwd=str(root),
        env=env,
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
        timeout=timeout_seconds,
        check=False,
    )
    elapsed = time.time() - start
    if proc.returncode != 0:
        raise RuntimeError(f"prediction failed with exit code {proc.returncode}: {proc.stderr or proc.stdout}")
    if not expected_output.exists():
        raise FileNotFoundError(f"prediction completed but output was not written: {expected_output}")
    if before_mtime is not None and expected_output.stat().st_mtime <= before_mtime:
        raise RuntimeError(f"prediction output was not refreshed: {expected_output}")

    payload = annotate_prediction_payload(json.loads(expected_output.read_text(encoding="utf-8")))
    issues = validate_prediction_payload(payload)
    if issues:
        detail = "; ".join(f"{i.field}: {i.message}" for i in issues)
        raise RuntimeError(f"prediction payload failed contract: {detail}")
    return PredictionRun(
        payload=payload,
        stdout=proc.stdout,
        stderr=proc.stderr,
        output_path=expected_output,
        elapsed_seconds=elapsed,
    )
