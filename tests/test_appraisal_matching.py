from __future__ import annotations

from tools.life_science_appraisal import block1_from_smiles, classify_expected_outcome, expected_in_candidate_pool, find_expected_row


def test_expected_match_uses_connectivity_not_display_name() -> None:
    naringenin = "O=C1C[C@@H](c2ccc(O)cc2)Oc2cc(O)cc(O)c21"
    expected_block1 = block1_from_smiles(naringenin)
    payload = {
        "top": [
            {
                "rank": 1,
                "product_name": "(R)-naringenin",
                "product_smiles": naringenin,
            }
        ]
    }
    rank, row, match_type = find_expected_row(payload, "naringenin", expected_block1)
    assert rank == 1
    assert row["product_name"] == "(R)-naringenin"
    assert match_type == "inchikey_block1"


def test_expected_block1_prevents_name_only_false_positive() -> None:
    naringenin = "O=C1C[C@@H](c2ccc(O)cc2)Oc2cc(O)cc(O)c21"
    expected_block1 = block1_from_smiles(naringenin)
    payload = {
        "top": [
            {
                "rank": 1,
                "product_name": "naringenin",
                "product_smiles": "CCO",
            }
        ]
    }
    rank, row, match_type = find_expected_row(payload, "naringenin", expected_block1)
    assert rank is None
    assert row == {}
    assert match_type == "miss"


def test_block1_from_smiles_tolerates_invalid_smiles() -> None:
    assert block1_from_smiles("not a smiles") is None


def test_expected_in_candidate_pool_uses_candidate_summary() -> None:
    payload = {"candidate_summary": {"all_product_blocks": ["AAAA", "BBBB"]}}
    assert expected_in_candidate_pool(payload, "BBBB") is True
    assert expected_in_candidate_pool(payload, "CCCC") is False


def test_failure_taxonomy_distinguishes_generation_from_ranking() -> None:
    payload = {
        "n_rule_candidates": 4,
        "candidate_summary": {"all_product_blocks": ["EXPECTED"]},
    }
    assert (
        classify_expected_outcome(payload, None, {}, "EXPECTED", "none")
        == "expected_generated_not_returned_topn"
    )
    assert (
        classify_expected_outcome(payload, None, {}, "MISSING", "none")
        == "expected_product_not_generated"
    )
    assert (
        classify_expected_outcome(payload, 2, {"biochem_quality": {"interpretation_allowed": True}}, "EXPECTED", "benchmark_panel")
        == "hit_top5_benchmark_only"
    )
