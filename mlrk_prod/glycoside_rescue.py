from __future__ import annotations

from contextlib import contextmanager

from rdkit import Chem
from rdkit import RDLogger
from rdkit.Chem import AllChem


AROMATIC_O_GLYCOSIDE_RESCUE_SMARTS = "[c:1]-[O:2]-[C;R;$(C(O)O):3]>>[c:1]-[O:2]"


@contextmanager
def _temporary_rdkit_silence():
    RDLogger.DisableLog("rdApp.*")
    try:
        yield
    finally:
        RDLogger.EnableLog("rdApp.*")


def _build_rescue_reaction():
    with _temporary_rdkit_silence():
        rxn = AllChem.ReactionFromSmarts(AROMATIC_O_GLYCOSIDE_RESCUE_SMARTS)
        if rxn is not None:
            rxn.Initialize()
    return rxn


_RESCUE_RXN = _build_rescue_reaction()


def aromatic_o_glycoside_rescue_candidates(smiles: str) -> list[str]:
    """Return conservative aglycone candidates for aromatic O-glycosides."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return []
    if _RESCUE_RXN is None:
        return []
    out: set[str] = set()
    try:
        with _temporary_rdkit_silence():
            product_sets = _RESCUE_RXN.RunReactants((mol,))
    except Exception:
        return []
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
    return sorted(out)
