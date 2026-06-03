from __future__ import annotations

import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def main() -> int:
    cmd = [
        sys.executable,
        str(ROOT / "modular" / "ltr_deploy.py"),
        "--candidates",
        "clean_candidates.parquet",
        "--force",
    ]
    proc = subprocess.run(cmd, cwd=str(ROOT), text=True, encoding="utf-8", errors="replace")
    if proc.returncode != 0:
        return proc.returncode
    validate = [sys.executable, str(ROOT / "tools" / "validate_production_readiness.py")]
    return subprocess.run(validate, cwd=str(ROOT), text=True, encoding="utf-8", errors="replace").returncode


if __name__ == "__main__":
    raise SystemExit(main())

