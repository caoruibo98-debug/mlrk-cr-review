from __future__ import annotations

from contextlib import contextmanager
from dataclasses import dataclass

from rdkit import Chem
from rdkit import RDLogger
from rdkit.Chem import AllChem


@dataclass(frozen=True)
class RescueRule:
    name: str
    smarts: str


MICROBIAL_RESCUE_RULES = [
    RescueRule(
        "hydroxycinnamate_side_chain_reduction",
        "[c:1][C:2]=[C:3][C:4](=O)[O:5]>>[c:1][C:2][C:3][C:4](=O)[O:5]",
    ),
    RescueRule(
        "stilbene_double_bond_reduction",
        "[c:1][C:2]=[C:3][c:4]>>[c:1][C:2][C:3][c:4]",
    ),
    RescueRule(
        "aromatic_carboxyl_decarboxylation",
        "[c:1][C:2](=O)[O:3]>>[c:1]",
    ),
    RescueRule(
        "phenol_dehydroxylation",
        "[c:1][O;H1:2]>>[c:1]",
    ),
]


@contextmanager
def _temporary_rdkit_silence():
    RDLogger.DisableLog("rdApp.*")
    try:
        yield
    finally:
        RDLogger.EnableLog("rdApp.*")


def _canonical(smiles: str | None) -> str | None:
    mol = Chem.MolFromSmiles(smiles or "")
    if mol is None:
        return None
    return Chem.MolToSmiles(mol)


def _compile_rules():
    compiled = []
    with _temporary_rdkit_silence():
        for rule in MICROBIAL_RESCUE_RULES:
            rxn = AllChem.ReactionFromSmarts(rule.smarts)
            if rxn is None:
                continue
            rxn.Initialize()
            compiled.append((rule, rxn))
    return compiled


_COMPILED_RULES = _compile_rules()


def _run_rule(mol: Chem.Mol, rxn) -> set[str]:
    out: set[str] = set()
    with _temporary_rdkit_silence():
        product_sets = rxn.RunReactants((mol,))
    for pset in product_sets:
        for product in pset:
            try:
                with _temporary_rdkit_silence():
                    Chem.SanitizeMol(product)
                    smi = Chem.MolToSmiles(product)
            except Exception:
                continue
            if smi:
                out.add(smi)
    return out


def microbial_rescue_candidates(smiles: str, max_steps: int = 2, max_products: int = 128) -> list[tuple[str, str]]:
    start = _canonical(smiles)
    if not start:
        return []
    seen = {start}
    frontier = {start}
    products: dict[str, str] = {}
    for step in range(1, max_steps + 1):
        next_frontier: set[str] = set()
        for current in sorted(frontier):
            mol = Chem.MolFromSmiles(current)
            if mol is None:
                continue
            for rule, rxn in _COMPILED_RULES:
                for product in _run_rule(mol, rxn):
                    canonical = _canonical(product)
                    if not canonical or canonical in seen:
                        continue
                    seen.add(canonical)
                    next_frontier.add(canonical)
                    products.setdefault(canonical, f"microbial_rescue:{rule.name}:step{step}")
                    if len(products) >= max_products:
                        return sorted(products.items())
        frontier = next_frontier
        if not frontier:
            break
    return sorted(products.items())
