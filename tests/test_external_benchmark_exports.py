from __future__ import annotations

import csv
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs" / "external_benchmarks"


def rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f))


def test_external_benchmark_manifest_and_files_exist() -> None:
    manifest = json.loads((OUT / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["case_count"] == 28
    assert manifest["panel_counts"] == {"challenge": 22, "core": 6}
    for rel in manifest["files"]:
        assert (OUT / rel).is_file(), rel


def test_external_benchmark_unified_table_has_expected_identity_fields() -> None:
    data = rows(OUT / "benchmark_panel_unified.csv")
    assert len(data) == 28
    first = data[0]
    for field in (
        "case_id",
        "substrate_smiles",
        "substrate_inchikey",
        "expected_product_smiles",
        "expected_product_inchikey",
        "reaction_family",
    ):
        assert first[field]


def test_biotransformer_tsv_preserves_case_ids() -> None:
    lines = (OUT / "biotransformer_input.tsv").read_text(encoding="utf-8").splitlines()
    assert len(lines) == 28
    case_id, smiles = lines[0].split("\t")
    assert case_id
    assert smiles


def test_microberx_and_gutbug_exports_have_distinct_scopes() -> None:
    microberx = rows(OUT / "microberx_query_table.csv")
    gutbug = rows(OUT / "gutbug_query_table.csv")
    assert len(microberx) == len(gutbug) == 28
    assert "expected_product_smiles" in microberx[0]
    assert "substrate_pubchem_cid" in gutbug[0]
