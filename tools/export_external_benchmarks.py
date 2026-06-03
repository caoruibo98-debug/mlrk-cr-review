from __future__ import annotations

import argparse
import csv
import json
import sys
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Any

from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "src"))

import kio  # noqa: E402
from resolve import name_to_smiles  # noqa: E402


CORE_PANEL = ROOT / "data" / "real_biochemistry_panel.csv"
CHALLENGE_PANEL = ROOT / "data" / "production_challenge_panel.csv"
OUT_DIR = ROOT / "outputs" / "external_benchmarks"


EXTERNAL_TOOL_SOURCES = {
    "biotransformer": {
        "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC9252798/",
        "input_note": "BioTransformer accepts a tab-separated structure input where an identifier can precede a SMILES or InChI, and can also take SDF.",
        "comparison_role": "Product-generation comparator in Human Gut Microbial, SuperBio, or MultiBio modes.",
    },
    "microberx": {
        "url": "https://microberx.readthedocs.io/en/stable/tutorials/PredictionMetabolites.html",
        "input_note": "MicrobeRX prediction is reaction-rule based; exported query and expected-product SMILES support product recall comparison after MicrobeRX runs.",
        "comparison_role": "Reaction-rule and EC/Rhea/PubMed evidence comparator.",
    },
    "gutbug": {
        "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC11810598/",
        "input_note": "GutBug/GutBugDB focuses on EC/enzyme and strain predictions; PubChem CIDs are exported when resolvable.",
        "comparison_role": "Enzyme/EC plausibility comparator rather than strict product-ranking comparator.",
    },
}


def read_panel(path: Path, panel_name: str) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as f:
        rows = list(csv.DictReader(f))
    for row in rows:
        row["panel"] = panel_name
    return rows


def inchikey(smiles: str | None) -> str | None:
    return kio.smiles_to_inchikey(smiles) if smiles else None


def block1(smiles: str | None) -> str | None:
    ik = inchikey(smiles)
    return kio.inchikey_block1(ik) if ik else None


def formula(smiles: str | None) -> str | None:
    mol = Chem.MolFromSmiles(smiles or "")
    return rdMolDescriptors.CalcMolFormula(mol) if mol is not None else None


def exact_mass(smiles: str | None) -> float | None:
    mol = Chem.MolFromSmiles(smiles or "")
    return round(float(Descriptors.ExactMolWt(mol)), 6) if mol is not None else None


def pubchem_cid(name: str, smiles: str | None, enabled: bool) -> str:
    if not enabled:
        return ""
    queries = [name]
    if smiles:
        queries.append(f"smiles/{urllib.parse.quote(smiles, safe='')}")
    for query in queries:
        if query.startswith("smiles/"):
            url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/{query}/cids/TXT"
        else:
            url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{urllib.parse.quote(query)}/cids/TXT"
        try:
            with urllib.request.urlopen(url, timeout=20) as response:
                first = response.read().decode("utf-8", errors="replace").splitlines()[0].strip()
            if first:
                return first
        except Exception:
            continue
    return ""


def normalized_rows(fetch_cids: bool) -> list[dict[str, Any]]:
    rows = read_panel(CORE_PANEL, "core") + read_panel(CHALLENGE_PANEL, "challenge")
    out: list[dict[str, Any]] = []
    for row in rows:
        substrate_smiles, substrate_source = name_to_smiles(row["substrate_name"])
        product_smiles, product_source = name_to_smiles(row["expected_product_name"])
        if not substrate_smiles or not product_smiles:
            raise RuntimeError(f"Could not resolve benchmark case {row['case_id']}: {row}")
        out.append(
            {
                "case_id": row["case_id"],
                "panel": row["panel"],
                "substrate_name": row["substrate_name"],
                "substrate_smiles": substrate_smiles,
                "substrate_inchikey": inchikey(substrate_smiles),
                "substrate_block1": block1(substrate_smiles),
                "substrate_formula": formula(substrate_smiles),
                "substrate_exact_mass": exact_mass(substrate_smiles),
                "substrate_pubchem_cid": pubchem_cid(row["substrate_name"], substrate_smiles, fetch_cids),
                "substrate_resolution_source": substrate_source,
                "expected_product_name": row["expected_product_name"],
                "expected_product_smiles": product_smiles,
                "expected_product_inchikey": inchikey(product_smiles),
                "expected_product_block1": block1(product_smiles),
                "expected_product_formula": formula(product_smiles),
                "expected_product_exact_mass": exact_mass(product_smiles),
                "expected_product_pubchem_cid": pubchem_cid(row["expected_product_name"], product_smiles, fetch_cids),
                "expected_product_resolution_source": product_source,
                "reaction_family": row.get("reaction_family", ""),
                "reference": row.get("reference", ""),
                "challenge_role": row.get("challenge_role", ""),
                "expected_difficulty": row.get("expected_difficulty", ""),
            }
        )
    return out


def write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_biotransformer(rows: list[dict[str, Any]]) -> None:
    with (OUT_DIR / "biotransformer_input.tsv").open("w", encoding="utf-8", newline="") as f:
        for row in rows:
            f.write(f"{row['case_id']}\t{row['substrate_smiles']}\n")
    writer = Chem.SDWriter(str(OUT_DIR / "biotransformer_input.sdf"))
    for row in rows:
        mol = Chem.MolFromSmiles(row["substrate_smiles"])
        if mol is None:
            continue
        mol.SetProp("_Name", row["case_id"])
        for prop in ("panel", "substrate_name", "expected_product_name", "reaction_family", "reference"):
            mol.SetProp(prop, str(row.get(prop, "")))
        writer.write(mol)
    writer.close()


def write_tool_tables(rows: list[dict[str, Any]]) -> None:
    common = [
        "case_id",
        "panel",
        "substrate_name",
        "substrate_smiles",
        "substrate_inchikey",
        "substrate_block1",
        "expected_product_name",
        "expected_product_smiles",
        "expected_product_inchikey",
        "expected_product_block1",
        "reaction_family",
        "reference",
    ]
    write_csv(
        OUT_DIR / "benchmark_panel_unified.csv",
        rows,
        common
        + [
            "substrate_formula",
            "substrate_exact_mass",
            "substrate_pubchem_cid",
            "expected_product_formula",
            "expected_product_exact_mass",
            "expected_product_pubchem_cid",
            "challenge_role",
            "expected_difficulty",
        ],
    )
    write_csv(
        OUT_DIR / "microberx_query_table.csv",
        rows,
        common
        + [
            "substrate_formula",
            "expected_product_formula",
            "challenge_role",
            "expected_difficulty",
        ],
    )
    write_csv(
        OUT_DIR / "gutbug_query_table.csv",
        rows,
        [
            "case_id",
            "panel",
            "substrate_name",
            "substrate_pubchem_cid",
            "substrate_smiles",
            "substrate_inchikey",
            "reaction_family",
            "expected_product_name",
            "expected_product_pubchem_cid",
            "reference",
        ],
    )
    with (OUT_DIR / "gutbug_pubchem_cids.txt").open("w", encoding="utf-8", newline="") as f:
        for row in rows:
            if row["substrate_pubchem_cid"]:
                f.write(f"{row['substrate_pubchem_cid']}\t{row['case_id']}\t{row['substrate_name']}\n")


def write_manifest(rows: list[dict[str, Any]], fetch_cids: bool) -> None:
    counts: dict[str, int] = {}
    for row in rows:
        counts[row["panel"]] = counts.get(row["panel"], 0) + 1
    manifest = {
        "export_name": "generation_11_external_benchmark_inputs",
        "case_count": len(rows),
        "panel_counts": dict(sorted(counts.items())),
        "fetch_pubchem_cids": fetch_cids,
        "files": [
            "benchmark_panel_unified.csv",
            "biotransformer_input.tsv",
            "biotransformer_input.sdf",
            "microberx_query_table.csv",
            "gutbug_query_table.csv",
            "gutbug_pubchem_cids.txt",
        ],
        "external_tool_sources": EXTERNAL_TOOL_SOURCES,
        "claim_boundary": "These are external benchmark inputs and expected-product labels; they are not external benchmark results.",
        "matching_key": "case_id plus expected_product_inchikey/expected_product_block1 for product-output tools; case_id plus substrate_pubchem_cid for GutBug-style EC/enzyme tools.",
    }
    (OUT_DIR / "manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--skip-pubchem-cids", action="store_true")
    args = parser.parse_args()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    rows = normalized_rows(fetch_cids=not args.skip_pubchem_cids)
    write_tool_tables(rows)
    write_biotransformer(rows)
    write_manifest(rows, fetch_cids=not args.skip_pubchem_cids)
    print(json.dumps({"out": str(OUT_DIR), "case_count": len(rows)}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
