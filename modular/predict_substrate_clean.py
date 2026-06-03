from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit import RDLogger

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import kio  # noqa: E402
import module_router as mr  # noqa: E402
from generate_candidates import run_reactants  # noqa: E402
from ltr_build import LTR_OUT  # noqa: E402
from ltr_candidates import load_module_rules  # noqa: E402
from mlrk_prod.glycoside_rescue import aromatic_o_glycoside_rescue_candidates  # noqa: E402


RDLogger.DisableLog("rdApp.*")

MODELS = LTR_OUT / "models_clean"
DEFAULT_EVIDENCE_POOL = Path(
    r"D:\CRB\FoodGut\positive_sample_db\data\processed\foodgut_modular_positive_candidate_pool_v2.csv"
)
EVIDENCE_POOL = Path(os.environ.get("MLRK_EVIDENCE_POOL", str(DEFAULT_EVIDENCE_POOL)))
HONEST_NOTE = (
    "Ranks only rule-generated candidates; ranking score is non-wet-lab and is not a biological "
    "occurrence probability; evidence is explanatory only and is not used for ranking."
)


def _l2(vector: np.ndarray) -> np.ndarray:
    return vector / max(float(np.linalg.norm(vector)), 1e-9)


def _safe_label(label: str) -> str:
    return "".join(c if c.isalnum() or c in "-_" else "_" for c in label)


def evidence_lookup() -> dict[tuple[str, str], dict[str, str | None]]:
    if not EVIDENCE_POOL.exists():
        return {}
    data = pd.read_csv(
        EVIDENCE_POOL,
        low_memory=False,
        usecols=[
            "substrate_inchikey",
            "product_inchikey",
            "enzyme_ec",
            "enzyme_name",
            "representative_microbes",
            "pmid",
        ],
    )
    data["sb"] = data.substrate_inchikey.astype(str).str[:14]
    data["pb"] = data.product_inchikey.astype(str).str[:14]
    evidence: dict[tuple[str, str], dict[str, str | None]] = {}
    for row in data.itertuples(index=False):
        evidence.setdefault(
            (row.sb, row.pb),
            {
                "enzyme_ec": None if pd.isna(row.enzyme_ec) else str(row.enzyme_ec),
                "enzyme": None if pd.isna(row.enzyme_name) else str(row.enzyme_name)[:50],
                "microbe": None if pd.isna(row.representative_microbes) else str(row.representative_microbes)[:50],
                "pmid": None if pd.isna(row.pmid) else str(row.pmid),
            },
        )
    return evidence


def write_prediction_payload(label: str, payload: dict) -> Path:
    path = kio.safe_output_path(f"modular/predictions/clean_{_safe_label(label)}.json")
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
    return path


def no_candidate_payload(label: str, smiles: str, module: str) -> dict:
    return {
        "input": {"name": label, "smiles": smiles},
        "module": module,
        "module_name": mr.MODULE_NAMES[module],
        "n_rule_candidates": 0,
        "candidate_generation_status": "no_candidates",
        "honest_note": HONEST_NOTE,
        "top": [],
    }


def generate_candidates(smiles: str, module: str, n_rules: int) -> tuple[dict[str, tuple[str, str]], dict[str, str]]:
    rules = load_module_rules(np.random.default_rng(0), n_rules)[module]
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise SystemExit("input SMILES cannot be parsed by RDKit")
    molh = Chem.AddHs(mol)
    substrate_block = kio.inchikey_block1(kio.smiles_to_inchikey(smiles))
    candidates: dict[str, tuple[str, str]] = {}
    sources: dict[str, str] = {}
    for smarts, _ec in rules:
        for product_smiles in run_reactants(smarts, mol, molh):
            product_inchikey = kio.smiles_to_inchikey(product_smiles)
            product_block = kio.inchikey_block1(product_inchikey)
            if product_block and product_block != substrate_block and product_block not in candidates:
                candidates[product_block] = (product_smiles, product_inchikey)
                sources[product_block] = "rule"
    if module == "B":
        for product_smiles in aromatic_o_glycoside_rescue_candidates(smiles):
            product_inchikey = kio.smiles_to_inchikey(product_smiles)
            product_block = kio.inchikey_block1(product_inchikey)
            if product_block and product_block != substrate_block and product_block not in candidates:
                candidates[product_block] = (product_smiles, product_inchikey)
                sources[product_block] = "curated_aromatic_o_glycoside_rescue"
    return candidates, sources


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--name", default=None)
    parser.add_argument("--smiles", default=None)
    parser.add_argument("--topn", type=int, default=12)
    parser.add_argument("--n-rules", type=int, default=800)
    args = parser.parse_args()
    kio.setup_logging()

    import xgboost as xgb
    from encode_molecules import FrozenEncoder
    from resolve import block1_to_name, name_to_smiles

    if args.smiles:
        smiles = kio.canonical_smiles(args.smiles)
        label = args.name or "query"
    elif args.name:
        smiles, source = name_to_smiles(args.name)
        if not smiles:
            raise SystemExit(f"could not resolve name: {args.name}")
        label = args.name
        kio.log.info("resolved [%s]: %s", source, smiles)
    else:
        raise SystemExit("provide --name or --smiles")

    module, _reason = mr.classify(smiles)
    kio.log.info("== %s -> module %s (%s) ==", label, module, mr.MODULE_NAMES[module])
    model_path = MODELS / f"ltr_chem_{module}.json"
    if not model_path.exists():
        raise SystemExit(f"missing deployed model for module {module}: {model_path}")
    meta = json.loads((MODELS / "deploy_meta.json").read_text(encoding="utf-8")).get(module, {})
    dim = int(meta.get("dim", 384))
    model = xgb.XGBRanker()
    model.load_model(str(model_path))

    candidates, candidate_sources = generate_candidates(smiles, module, args.n_rules)
    if not candidates:
        payload = no_candidate_payload(label, smiles, module)
        write_prediction_payload(label, payload)
        print(json.dumps(payload, ensure_ascii=False, indent=2))
        return
    kio.log.info("generated candidates: %d", len(candidates))

    encoder = FrozenEncoder(kio.load_config()["encoder"]["primary"])
    candidate_smiles = [value[0] for value in candidates.values()]
    encoded_smiles = [smiles] + candidate_smiles
    vectors = encoder.encode(encoded_smiles)
    embeddings = {smi: np.asarray(vec, np.float32) for smi, vec in zip(encoded_smiles, vectors)}
    substrate_embedding = _l2(embeddings[smiles])

    product_blocks = list(candidates.keys())
    features = np.zeros((len(product_blocks), 4 * dim), np.float32)
    for idx, product_block in enumerate(product_blocks):
        product_embedding = _l2(embeddings[candidates[product_block][0]])
        features[idx] = np.concatenate(
            [
                substrate_embedding,
                product_embedding,
                substrate_embedding - product_embedding,
                substrate_embedding * product_embedding,
            ]
        )
    scores = model.predict(features)

    substrate_block = kio.inchikey_block1(kio.smiles_to_inchikey(smiles))
    evidence = evidence_lookup()
    ranked = pd.DataFrame(
        {
            "product_block": product_blocks,
            "smiles": [candidates[product_block][0] for product_block in product_blocks],
            "score": scores,
            "candidate_source": [candidate_sources.get(product_block, "rule") for product_block in product_blocks],
        }
    )
    ranked = ranked.sort_values("score", ascending=False).head(args.topn).reset_index(drop=True)
    low, high = ranked.score.min(), ranked.score.max()

    rows = []
    for idx, row in ranked.iterrows():
        ev = evidence.get((substrate_block, row.product_block), {})
        evidence_text = ";".join(f"{key}={value}" for key, value in ev.items() if value) or "-"
        score = round(float((row.score - low) / (high - low)) if high > low else 0.5, 3)
        rows.append(
            {
                "rank": idx + 1,
                "product_name": block1_to_name(row.product_block) or "unnamed/novel",
                "product_smiles": row.smiles,
                "score": score,
                "evidence": evidence_text,
                "candidate_source": row.candidate_source,
            }
        )

    payload = {
        "input": {"name": label, "smiles": smiles},
        "module": module,
        "module_name": mr.MODULE_NAMES[module],
        "n_rule_candidates": len(candidates),
        "candidate_generation_status": "generated_candidates",
        "honest_note": HONEST_NOTE,
        "top": rows,
    }
    write_prediction_payload(label, payload)
    print(f"\n=== {label} -> module {module} top-{args.topn} (candidates={len(candidates)}) ===")
    print(pd.DataFrame(rows)[["rank", "product_name", "score", "evidence"]].to_string(index=False, max_colwidth=40))


if __name__ == "__main__":
    main()
