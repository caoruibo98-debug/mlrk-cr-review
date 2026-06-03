from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


TRACKED = [
    "production_artifact_manifest.json",
    "pyproject.toml",
    "mlrk_prod",
    "tools",
    "tests",
    "docs",
    "data",
]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def iter_files() -> list[Path]:
    files: list[Path] = []
    for item in TRACKED:
        p = ROOT / item
        if p.is_file():
            files.append(p)
        elif p.is_dir():
            files.extend(x for x in p.rglob("*") if x.is_file() and "__pycache__" not in x.parts)
    return sorted(files)


def git_commit() -> str | None:
    proc = subprocess.run(["git", "rev-parse", "HEAD"], cwd=str(ROOT), capture_output=True, text=True, check=False)
    if proc.returncode != 0:
        return None
    return proc.stdout.strip()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--generation", required=True)
    parser.add_argument("--review", default="")
    parser.add_argument("--score-file", default="outputs/appraisal/life_science_appraisal.json")
    args = parser.parse_args()

    freeze_dir = ROOT / "freezes" / args.generation
    freeze_dir.mkdir(parents=True, exist_ok=True)
    score_path = ROOT / args.score_file
    files = iter_files()
    manifest = {
        "generation": args.generation,
        "git_commit": git_commit(),
        "review": args.review,
        "file_count": len(files),
        "files": [{"path": str(p.relative_to(ROOT)).replace("\\", "/"), "sha256": sha256(p)} for p in files],
        "score_file": str(score_path.relative_to(ROOT)).replace("\\", "/") if score_path.exists() else None,
    }
    if score_path.exists():
        manifest["score"] = json.loads(score_path.read_text(encoding="utf-8")).get("score")
    (freeze_dir / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(f"froze {args.generation} with {len(files)} tracked files")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

