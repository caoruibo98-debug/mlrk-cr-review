from __future__ import annotations

import json
import tempfile
from pathlib import Path

from mlrk_prod.io_utils import atomic_write_json, atomic_write_text


def test_atomic_write_json_writes_complete_json() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp) / "report.json"
        atomic_write_json(path, {"status": "ready", "value": 1})
        assert json.loads(path.read_text(encoding="utf-8")) == {"status": "ready", "value": 1}
        leftovers = list(Path(tmp).glob("*.tmp"))
        assert leftovers == []


def test_atomic_write_text_replaces_existing_file() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp) / "report.md"
        atomic_write_text(path, "old")
        atomic_write_text(path, "new")
        assert path.read_text(encoding="utf-8") == "new"
