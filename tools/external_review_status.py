from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from mlrk_prod.io_utils import atomic_write_json  # noqa: E402

DEFAULT_OUT = ROOT / "outputs" / "appraisal" / "external_review_status.json"

FINAL_STATUSES = {"completed", "timed_out", "failed", "skipped"}


def build_report(
    *,
    tool: str,
    status: str,
    command: str,
    duration_ms: int | None,
    issue_count: int | None,
    note: str,
) -> dict:
    if status not in FINAL_STATUSES:
        raise ValueError(f"status must be one of {sorted(FINAL_STATUSES)}")
    completed = status == "completed"
    return {
        "status": status,
        "tool": tool,
        "command": command,
        "duration_ms": duration_ms,
        "issue_count": issue_count if completed else None,
        "review_claim_allowed": completed,
        "note": note,
        "next_step": (
            "Summarize and address returned review issues."
            if completed
            else "Retry external AI/code review later or use a different review service; do not claim external review passed."
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--tool", required=True)
    parser.add_argument("--status", required=True, choices=sorted(FINAL_STATUSES))
    parser.add_argument("--command", required=True)
    parser.add_argument("--duration-ms", type=int, default=None)
    parser.add_argument("--issue-count", type=int, default=None)
    parser.add_argument("--note", default="")
    parser.add_argument("--out", default=str(DEFAULT_OUT.relative_to(ROOT)).replace("\\", "/"))
    args = parser.parse_args()

    report = build_report(
        tool=args.tool,
        status=args.status,
        command=args.command,
        duration_ms=args.duration_ms,
        issue_count=args.issue_count,
        note=args.note,
    )
    out = ROOT / args.out
    atomic_write_json(out, report)
    print(json.dumps({"status": report["status"], "tool": report["tool"], "out": str(out.relative_to(ROOT)).replace("\\", "/")}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
