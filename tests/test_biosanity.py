from __future__ import annotations

from mlrk_prod.biosanity import candidate_quality
from mlrk_prod.biosanity import annotate_prediction_payload
from mlrk_prod.glycoside_rescue import aromatic_o_glycoside_rescue_candidates
from mlrk_prod.microbial_rescue import microbial_rescue_candidates
from rdkit import Chem


def canonical(smiles: str) -> str:
    return Chem.MolToSmiles(Chem.MolFromSmiles(smiles))


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
    assert quality["interpretation_allowed"] is False
    assert any(flag["code"] == "phosphorus_introduced" for flag in quality["flags"])


def test_annotated_payload_excludes_rejected_candidates_from_interpretation_ready_view() -> None:
    payload = {
        "input": {"smiles": "O=C1C[C@@H](c2ccc(O)cc2)Oc2cc(O)cc(O)c21"},
        "module_name": "small_molecule",
        "top": [
            {
                "rank": 1,
                "product_name": "bad phosphate",
                "product_smiles": "COP(=O)(OC)Oc1ccc([C@@H]2CC(=O)c3c(O)cc(O)cc3O2)cc1",
                "score": 1.0,
                "evidence": "-",
            },
            {
                "rank": 2,
                "product_name": "naringenin",
                "product_smiles": "O=C1C[C@@H](c2ccc(O)cc2)Oc2cc(O)cc(O)c21",
                "score": 0.5,
                "evidence": "pmid=example",
            },
        ],
    }
    out = annotate_prediction_payload(payload)
    assert [row["rank"] for row in out["interpretation_ready_top"]] == [2]
    assert out["quality_summary"]["rejected_ranks"] == [1]


def test_annotated_payload_skips_malformed_rejected_rank() -> None:
    payload = {
        "input": {"smiles": "O=C1C[C@@H](c2ccc(O)cc2)Oc2cc(O)cc(O)c21"},
        "module_name": "small_molecule",
        "top": [
            {
                "rank": "not-a-rank",
                "product_name": "bad phosphate",
                "product_smiles": "COP(=O)(OC)Oc1ccc([C@@H]2CC(=O)c3c(O)cc(O)cc3O2)cc1",
                "score": 1.0,
                "evidence": "-",
            }
        ],
    }
    out = annotate_prediction_payload(payload)
    assert out["interpretation_ready_top"] == []
    assert out["quality_summary"]["rejected_ranks"] == []


def test_aromatic_glycoside_rescue_finds_quercitrin_aglycone() -> None:
    quercitrin = "C[C@@H]1O[C@@H](Oc2c(-c3ccc(O)c(O)c3)oc3cc(O)cc(O)c3c2=O)[C@H](O)[C@H](O)[C@H]1O"
    products = aromatic_o_glycoside_rescue_candidates(quercitrin)
    assert "O=c1c(O)c(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12" in products


def test_aromatic_glycoside_rescue_does_not_fire_on_aglycone_or_caffeine() -> None:
    quercetin = "O=c1c(O)c(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12"
    caffeine = "Cn1c(=O)c2c(ncn2C)n(C)c1=O"
    assert aromatic_o_glycoside_rescue_candidates(quercetin) == []
    assert aromatic_o_glycoside_rescue_candidates(caffeine) == []


def test_microbial_rescue_generates_hydroxycinnamate_reduction() -> None:
    caffeic_acid = "O=C(O)/C=C/c1ccc(O)c(O)c1"
    expected = canonical("O=C(O)CCc1ccc(O)c(O)c1")
    products = {product for product, _source in microbial_rescue_candidates(caffeic_acid)}
    assert expected in products


def test_microbial_rescue_generates_decarboxylation() -> None:
    gallic_acid = "O=C(O)c1cc(O)c(O)c(O)c1"
    expected = canonical("Oc1cccc(O)c1O")
    products = {product for product, _source in microbial_rescue_candidates(gallic_acid)}
    assert expected in products


def test_microbial_rescue_generates_two_step_lunularin() -> None:
    resveratrol = "Oc1ccc(/C=C/c2cc(O)cc(O)c2)cc1"
    expected = canonical("Oc1ccc(CCc2cccc(O)c2)cc1")
    products = {product for product, _source in microbial_rescue_candidates(resveratrol)}
    assert expected in products
