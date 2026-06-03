from __future__ import annotations

import copy
from dataclasses import dataclass
from typing import Any

from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors


SUGAR_LOSS_DELTAS = (146.0, 162.0, 176.0, 308.0)
EMPTY_EVIDENCE = {"", "-", "--", "---", "—", "鈥?", "бк"}


@dataclass(frozen=True)
class QualityFlag:
    code: str
    severity: str
    message: str


def _mol(smiles: str | None) -> Chem.Mol | None:
    if not smiles:
        return None
    return Chem.MolFromSmiles(str(smiles))


def _atom_count(mol: Chem.Mol, atomic_num: int) -> int:
    return sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == atomic_num)


def _evidence_is_empty(value: Any) -> bool:
    text = str(value or "").strip()
    if text in EMPTY_EVIDENCE:
        return True
    return not any(token in text.lower() for token in ("enzyme=", "pmid=", "microbe=", "ec=", "known:"))


def _rank_as_int(value: Any) -> int | None:
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def _tier(score: float) -> str:
    if score >= 0.75:
        return "high"
    if score >= 0.55:
        return "medium"
    if score >= 0.35:
        return "low"
    return "reject"


def candidate_quality(input_smiles: str, candidate: dict[str, Any], module_name: str = "") -> dict[str, Any]:
    sub = _mol(input_smiles)
    prod = _mol(candidate.get("product_smiles"))
    flags: list[QualityFlag] = []
    score = 1.0

    if sub is None:
        flags.append(QualityFlag("invalid_input_smiles", "error", "Input SMILES cannot be parsed by RDKit."))
        score = 0.0
    if prod is None:
        flags.append(QualityFlag("invalid_product_smiles", "error", "Product SMILES cannot be parsed by RDKit."))
        score = 0.0

    result: dict[str, Any] = {
        "valid_smiles": prod is not None,
        "score": 0.0,
        "tier": "reject",
        "interpretation_allowed": False,
        "flags": [flag.__dict__ for flag in flags],
    }
    if sub is None or prod is None:
        return result

    sub_mass = Descriptors.ExactMolWt(sub)
    prod_mass = Descriptors.ExactMolWt(prod)
    mass_delta = prod_mass - sub_mass
    sub_c = _atom_count(sub, 6)
    prod_c = _atom_count(prod, 6)
    sub_p = _atom_count(sub, 15)
    prod_p = _atom_count(prod, 15)
    formal_charge = sum(atom.GetFormalCharge() for atom in prod.GetAtoms())
    heavy_atoms = prod.GetNumHeavyAtoms()
    evidence_missing = _evidence_is_empty(candidate.get("evidence"))

    if evidence_missing:
        flags.append(QualityFlag("no_prior_evidence", "warning", "No enzyme, microbe, EC, or PMID evidence is attached."))
        score -= 0.25
    if heavy_atoms < 3 or prod_mass < 45:
        flags.append(QualityFlag("small_byproduct", "error", "Product is too small to present as the main metabolite target."))
        score -= 0.55
    if formal_charge != 0:
        flags.append(QualityFlag("charged_candidate", "warning", "Product has a non-zero formal charge; neutral form/adduct handling needs review."))
        score -= 0.15
    if prod_p > sub_p:
        flags.append(QualityFlag("phosphorus_introduced", "error", "Product introduces phosphorus absent from the substrate."))
        score -= 0.45
    if prod_c > sub_c + 6:
        flags.append(QualityFlag("large_carbon_gain", "error", "Product gains many carbon atoms relative to the substrate."))
        score -= 0.50
    if mass_delta > 80:
        flags.append(QualityFlag("large_mass_gain", "error", "Product mass gain is large for an unconstrained gut microbial candidate."))
        score -= 0.40
    if abs(mass_delta) > 500:
        flags.append(QualityFlag("extreme_mass_delta", "warning", "Mass delta is extreme and needs manual biochemical review."))
        score -= 0.20

    if "glycoside" in module_name.lower() and mass_delta < -80:
        if any(abs(abs(mass_delta) - expected) <= 15 for expected in SUGAR_LOSS_DELTAS):
            flags.append(QualityFlag("plausible_glycoside_mass_loss", "info", "Mass delta is consistent with common sugar loss."))
            score += 0.10

    score = max(0.0, min(1.0, score))
    tier = _tier(score)
    has_error_flag = any(flag.severity == "error" for flag in flags)
    result.update(
        {
            "valid_smiles": True,
            "score": round(score, 3),
            "tier": tier,
            "interpretation_allowed": tier in {"high", "medium"} and not has_error_flag,
            "exact_mass": round(prod_mass, 5),
            "formula": rdMolDescriptors.CalcMolFormula(prod),
            "mass_delta_from_input": round(mass_delta, 5),
            "formal_charge": formal_charge,
            "flags": [flag.__dict__ for flag in flags],
        }
    )
    return result


def annotate_prediction_payload(payload: dict[str, Any]) -> dict[str, Any]:
    out = copy.deepcopy(payload)
    input_smiles = str(out.get("input", {}).get("smiles", ""))
    module_name = str(out.get("module_name", ""))
    tier_counts: dict[str, int] = {}
    rejected_ranks: list[int] = []
    for row in out.get("top", []):
        row["biochem_quality"] = candidate_quality(input_smiles, row, module_name=module_name)
        tier = row["biochem_quality"]["tier"]
        tier_counts[tier] = tier_counts.get(tier, 0) + 1
        if not row["biochem_quality"]["interpretation_allowed"]:
            rank = _rank_as_int(row.get("rank"))
            if rank is not None:
                rejected_ranks.append(rank)
    out["interpretation_ready_top"] = [
        row for row in out.get("top", []) if row.get("biochem_quality", {}).get("interpretation_allowed")
    ]
    out["quality_summary"] = {
        "returned_candidates": len(out.get("top", [])),
        "interpretation_ready_candidates": len(out["interpretation_ready_top"]),
        "tier_counts": tier_counts,
        "rejected_ranks": rejected_ranks,
    }
    out["application_boundary"] = {
        "score_meaning": "Ranking score is not a wet-lab probability.",
        "quality_meaning": "Biochemical quality flags support review; they do not prove occurrence.",
        "allowed_use": "Internal research prioritization and LC-MS target triage.",
        "forbidden_use": "Consumer, clinical, or causal microbiome health claims.",
    }
    return out
