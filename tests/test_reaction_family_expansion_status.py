from __future__ import annotations

import csv
import json
from pathlib import Path

from tools.reaction_family_expansion_status import REQUIRED_SEED_FIELDS, build_report, read_seed


ROOT = Path(__file__).resolve().parents[1]
SEED = ROOT / "data" / "reaction_family_gap_seed_manifest.tsv"
KPIS = ROOT / "outputs" / "appraisal" / "reaction_family_kpis.json"
STATUS = ROOT / "outputs" / "appraisal" / "reaction_family_expansion_status.json"
STATUS_MD = ROOT / "docs" / "reviews" / "reaction_family_expansion_status.md"


def test_gap_seed_manifest_has_required_fields_and_machine_identifiers() -> None:
    with SEED.open("r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        assert reader.fieldnames is not None
        for field in REQUIRED_SEED_FIELDS:
            assert field in reader.fieldnames
        rows = list(reader)
    assert len(rows) >= 14
    for row in rows:
        assert row["substrate_smiles"]
        assert row["product_smiles"]
        assert row["substrate_inchikey"]
        assert row["product_inchikey"]
        assert row["substrate_name_en"]
        assert row["substrate_name_zh"]
        assert row["product_name_en"]
        assert row["product_name_zh"]


def test_gap_seed_manifest_keeps_challenge_cases_out_of_training() -> None:
    rows = read_seed(SEED)
    assert rows
    assert all(row["split_role"] == "benchmark_holdout" for row in rows)
    assert all(row["training_allowed"].lower() == "false" for row in rows)
    assert all("do_not_train" in row["holdout_policy"] for row in rows)


def test_expansion_status_covers_all_generation_blocked_families() -> None:
    report = build_report(read_seed(SEED), json.loads(KPIS.read_text(encoding="utf-8")))
    assert report["status"] == "ready"
    assert report["blocked_family_count"] >= 10
    assert report["blocked_without_seed"] == []
    assert report["ready_for_retraining_families"] == []
    assert report["policy"]["challenge_holdouts_are_not_training_data"] is True


def test_expansion_status_artifacts_are_written() -> None:
    report = json.loads(STATUS.read_text(encoding="utf-8"))
    text = STATUS_MD.read_text(encoding="utf-8")
    assert report["status"] == "ready"
    assert report["seed_summary"]["holdout_rows"] >= 14
    assert "Challenge holdouts are training data: `false`" in text
