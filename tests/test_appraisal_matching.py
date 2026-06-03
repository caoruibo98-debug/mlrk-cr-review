from __future__ import annotations

from tools.life_science_appraisal import block1_from_smiles, find_expected_row


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
