from __future__ import annotations

from rdkit import Chem
from rdkit import RDLogger
from rdkit.Chem import AllChem


AROMATIC_O_GLYCOSIDE_RESCUE_SMARTS = "[c:1]-[O:2]-[C;R;$(C(O)O):3]>>[c:1]-[O:2]"

RDLogger.DisableLog("rdApp.*")


def aromatic_o_glycoside_rescue_candidates(smiles: str) -> list[str]:
    """Return conservative aglycone candidates for aromatic O-glycosides."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return []
    rxn = AllChem.ReactionFromSmarts(AROMATIC_O_GLYCOSIDE_RESCUE_SMARTS)
    if rxn is None:
        return []
    rxn.Initialize()
    out: set[str] = set()
    try:
        product_sets = rxn.RunReactants((mol,))
    except Exception:
        return []
    for pset in product_sets:
        for product in pset:
            try:
                Chem.SanitizeMol(product)
                smi = Chem.MolToSmiles(product)
            except Exception:
                continue
            if smi:
                out.add(smi)
    return sorted(out)
