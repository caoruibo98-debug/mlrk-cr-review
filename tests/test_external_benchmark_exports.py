from __future__ import annotations

import csv
import json
from pathlib import Path

from tools.score_external_results import (
    expected_map,
    load_enzyme_results,
    load_product_results,
    score_product_results,
    summarize_enzyme_results,
)


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs" / "external_benchmarks"
FIXTURES = ROOT / "tests" / "fixtures"


def rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f))


def test_external_benchmark_manifest_and_files_exist() -> None:
    manifest = json.loads((OUT / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["case_count"] == 28
    assert manifest["panel_counts"] == {"challenge": 22, "core": 6}
    for rel in manifest["files"]:
        assert (OUT / rel).is_file(), rel


def test_external_result_scorecard_waits_for_real_outputs() -> None:
    scorecard = json.loads((OUT / "external_result_scorecard.json").read_text(encoding="utf-8"))
    assert scorecard["status"] == "awaiting_external_outputs"
    assert scorecard["product_tool_summaries"] == {}
    for rel in scorecard["template_files"]:
        assert (ROOT / rel).is_file(), rel


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


def test_external_product_result_scoring_uses_case_id_and_inchikey() -> None:
    expected = expected_map()
    product_rows = load_product_results([FIXTURES / "external_product_results_sample.csv"])
    scored = score_product_results(product_rows, expected)
    summary = scored["tool_summaries"]["sample_product_tool"]
    assert summary["case_count"] == 28
    assert summary["top1_block1_recall"] == round(1 / 28, 3)
    assert summary["top5_block1_recall"] == round(1 / 28, 3)
    assert summary["any_rank_block1_recall"] == round(2 / 28, 3)
    rutin = [
        row
        for row in scored["case_results"]
        if row["tool"] == "sample_product_tool" and row["case_id"] == "rutin_quercetin"
    ][0]
    assert rutin["top5_block1_hit"] is True


def test_external_enzyme_result_summary_is_not_product_recall() -> None:
    expected = expected_map()
    enzyme_rows = load_enzyme_results([FIXTURES / "external_enzyme_results_sample.csv"])
    summary = summarize_enzyme_results(enzyme_rows, expected)["sample_ec_tool"]
    assert summary["case_count"] == 28
    assert summary["cases_with_enzyme_or_ec_predictions"] == 2
    assert "not product recall" in summary["scoring_note"]
