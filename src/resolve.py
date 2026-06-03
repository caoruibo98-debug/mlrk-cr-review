"""
resolve.py — 名称 → SMILES 解析（给 predict_substrate 用名字代替 SMILES）

顺序：
  1) PubChem PUG-REST 在线（名字/同义词/缩写都能解，如 "2'-FL" → 2'-fucosyllactose）
  2) 本地 ChEBI / FooDB 全名（离线兜底）
  3) 小同义词表（少数缩写补充）
返回 (canonical_smiles, source) 或 (None, reason)。
"""
from __future__ import annotations

import json
import sys
import urllib.parse
import urllib.request
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import kio

SYNONYMS = {  # 缩写 → 全名（PubChem 多数能直接解，这里只作额外保险）
    "2'-fl": "2'-fucosyllactose", "3-fl": "3-fucosyllactose",
    "3'-sl": "3'-sialyllactose", "6'-sl": "6'-sialyllactose",
    "lnt": "lacto-N-tetraose", "lnnt": "lacto-N-neotetraose",
    "dfl": "difucosyllactose",
}


def _pubchem(name: str) -> str | None:
    url = (f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/"
           f"{urllib.parse.quote(name)}/property/IsomericSMILES/JSON")
    try:
        with urllib.request.urlopen(url, timeout=20) as r:
            props = json.load(r)["PropertyTable"]["Properties"][0]
        for k, v in props.items():
            if "SMILES" in k.upper() and isinstance(v, str) and v:
                return v
    except Exception:
        return None
    return None


_CHEBI = None
def _local_chebi(name: str) -> str | None:
    global _CHEBI
    p = kio.PROJECT_ROOT / "Dataset" / "chebi" / "chebi_compounds.csv"
    if not p.exists():
        return None
    if _CHEBI is None:
        _CHEBI = pd.read_csv(p, usecols=["name", "smiles"], dtype=str).dropna()
        _CHEBI["_n"] = _CHEBI["name"].str.strip().str.lower()
    hit = _CHEBI[_CHEBI["_n"] == name.strip().lower()]
    return hit["smiles"].iloc[0] if len(hit) else None


def name_to_smiles(name: str) -> tuple[str | None, str]:
    raw = name.strip()
    key = raw.lower().replace("’", "'")  # 统一弯引号
    queries = [raw]
    if key in SYNONYMS:
        queries.append(SYNONYMS[key])

    for q in queries:
        smi = _pubchem(q)
        if smi:
            c = kio.canonical_smiles(smi)
            if c:
                return c, f"pubchem:{q}"
    for q in queries:
        smi = _local_chebi(q)
        if smi:
            c = kio.canonical_smiles(smi)
            if c:
                return c, f"chebi:{q}"
    return None, "unresolved"


# ============================================================
# 反向：结构(InChIKey) → 名称
# ============================================================

_NAME_MAP = None
def load_name_map() -> dict[str, str]:
    global _NAME_MAP
    if _NAME_MAP is None:
        p = kio.FEATURES_DIR / "inchikey_name_map.parquet"
        if p.exists():
            d = pd.read_parquet(p)
            _NAME_MAP = dict(zip(d["block1"], d["name"]))
        else:
            _NAME_MAP = {}
    return _NAME_MAP


def block1_to_name(block1: str) -> str | None:
    return load_name_map().get(block1)


def ik_to_name_pubchem(inchikey: str) -> str | None:
    if not inchikey:
        return None
    url = (f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/"
           f"{urllib.parse.quote(inchikey)}/property/Title/JSON")
    try:
        with urllib.request.urlopen(url, timeout=15) as r:
            t = json.load(r)["PropertyTable"]["Properties"][0].get("Title")
            return t or None
    except Exception:
        return None


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("name")
    a = ap.parse_args()
    smi, src = name_to_smiles(a.name)
    print(f"{a.name} -> {smi}  [{src}]")
