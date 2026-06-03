from __future__ import annotations

from mlrk_prod.biosanity import candidate_quality


def test_rutin_to_quercetin_quality_is_high() -> None:
    rutin = "C[C@@H]1O[C@@H](OC[C@H]2O[C@@H](Oc3c(-c4ccc(O)c(O)c4)oc4cc(O)cc(O)c4c3=O)[C@H](O)[C@@H](O)[C@@H]2O)[C@H](O)[C@H](O)[C@H]1O"
    row = {
        "product_smiles": "O=c1c(O)c(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12",
        "evidence": "enzyme=alpha-L-rhamnosidase + beta-glucosidase;pmid=9875509",
    }
    quality = candidate_quality(rutin, row, module_name="carb_glycoside")
    assert quality["tier"] == "high"
    assert any(flag["code"] == "plausible_glycoside_mass_loss" for flag in quality["flags"])


def test_phosphorus_introduction_is_rejected_for_flavonoid() -> None:
    naringenin = "O=C1C[C@@H](c2ccc(O)cc2)Oc2cc(O)cc(O)c21"
    row = {
        "product_smiles": "COP(=O)(OC)Oc1ccc([C@@H]2CC(=O)c3c(O)cc(O)cc3O2)cc1",
        "evidence": "-",
    }
    quality = candidate_quality(naringenin, row, module_name="small_molecule")
    assert quality["tier"] == "reject"
    assert any(flag["code"] == "phosphorus_introduced" for flag in quality["flags"])

