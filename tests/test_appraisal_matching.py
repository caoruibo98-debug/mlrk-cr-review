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

