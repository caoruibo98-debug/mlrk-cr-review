from __future__ import annotations

import json
import sys
from dataclasses import asdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from mlrk_prod.manifest import readiness_status, validate_readiness  # noqa: E402


def main() -> int:
    issues = validate_readiness()
    status = readiness_status(issues)
    print(json.dumps({"status": status, "issues": [asdict(i) for i in issues]}, indent=2))
    return 1 if status == "blocked" else 0


if __name__ == "__main__":
    raise SystemExit(main())

