#!/usr/bin/env python
"""Round9 chemical mapping and split-risk audit.

This round checks whether Round8 exact-positive candidates are genuinely new
labels or already-existing positives with candidate-generation misses.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import pandas as pd
from rdkit import Chem
from rdkit.Chem import inchi


REPO = Path(__file__).resolve().parents[1]
CURATION = REPO / "data" / "curation"
PUBCHEM_RAW = CURATION / "pubchem_round9_raw"
RAW_POOL = Path(r"D:\CRB\FoodGut\positive_sample_db\data\processed\foodgut_modular_positive_candidate_pool_v2.csv")

ROUND8_POS = CURATION / "positive_sample_expansion_round8.csv"
ROUND8_NEG = CURATION / "negative_evidence_candidates_round8.csv"
REACTIONS = REPO / "outputs" / "modular" / "ltr" / "reactions.parquet"
CLEAN_FULL = REPO / "outputs" / "modular" / "ltr" / "clean_candidates_full.parquet"
CLEAN = REPO / "outputs" / "modular" / "ltr" / "clean_candidates.parquet"
PAIRS = REPO / "outputs" / "modular" / "ltr" / "pairs.parquet"
CHALLENGE = REPO / "data" / "production_challenge_panel.csv"

OUT_MAPPING = CURATION / "chemical_mapping_round9.csv"
OUT_IMPORT = CURATION / "positive_sample_import_round9.csv"
OUT_SPLIT = CURATION / "split_and_generation_risk_round9.csv"
OUT_NEG = CURATION / "negative_label_review_round9.csv"
OUT_QUEUE = CURATION / "round9_pipeline_queue.csv"


def b1(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    return str(value).split("-", 1)[0]


def pair_key(sb: object, pb: object) -> str:
    return f"{sb}__{pb}"


def smiles_to_inchikey(smiles: object) -> str:
    if smiles is None or pd.isna(smiles):
        return ""
    mol = Chem.MolFromSmiles(str(smiles))
    if mol is None:
        return ""
    return inchi.MolToInchiKey(mol)


def smiles_to_canonical(smiles: object, isomeric: bool = True) -> str:
    if smiles is None or pd.isna(smiles):
        return ""
    mol = Chem.MolFromSmiles(str(smiles))
    if mol is None:
        return ""
    return Chem.MolToSmiles(mol, isomericSmiles=isomeric)


def safe_name(name: str) -> str:
    return re.sub(r"[^A-Za-z0-9]+", "_", name).strip("_").lower()


def load_pubchem(name: str) -> dict[str, object]:
    path = PUBCHEM_RAW / f"{safe_name(name)}.json"
    if not path.exists():
        return {}
    raw = json.loads(path.read_text(encoding="utf-8"))
    rows = raw.get("PropertyTable", {}).get("Properties", [])
    return rows[0] if rows else {}


def name_match_mask(df: pd.DataFrame, substrate: str, product: str) -> pd.Series:
    text_sub = substrate.lower()
    text_prod = product.lower()
    cols = [c for c in ["substrate_name", "expected_product_name", "product_name"] if c in df.columns]
    if not cols:
        return pd.Series([False] * len(df), index=df.index)
    joined = df[cols].astype(str).agg(" ".join, axis=1).str.lower()
    return joined.str.contains(re.escape(text_sub), regex=True) & joined.str.contains(re.escape(text_prod), regex=True)


def build_raw_pair_keys(raw: pd.DataFrame) -> pd.DataFrame:
    out = raw.copy()
    if "substrate_inchikey" in out.columns and "product_inchikey" in out.columns:
        out["sb"] = out["substrate_inchikey"].map(b1)
        out["pb"] = out["product_inchikey"].map(b1)
        out["pair_key"] = out["sb"] + "__" + out["pb"]
    return out


def table_with_pair_key(path: Path) -> pd.DataFrame:
    df = pd.read_parquet(path)
    df["pair_key"] = df["sb"].astype(str) + "__" + df["pb"].astype(str)
    return df


def main() -> None:
    pos = pd.read_csv(ROUND8_POS)
    raw = build_raw_pair_keys(pd.read_csv(RAW_POOL, low_memory=False))
    reactions = table_with_pair_key(REACTIONS)
    clean_full = table_with_pair_key(CLEAN_FULL)
    clean = table_with_pair_key(CLEAN)
    pairs = table_with_pair_key(PAIRS)
    challenge = pd.read_csv(CHALLENGE) if CHALLENGE.exists() else pd.DataFrame()

    mapping_rows: list[dict[str, object]] = []
    import_rows: list[dict[str, object]] = []
    split_rows: list[dict[str, object]] = []
    queue_rows: list[dict[str, object]] = []

    for item in pos.itertuples(index=False):
        key = item.pair_key
        rxn_hit = reactions[reactions["pair_key"] == key]
        raw_hit = raw[raw.get("pair_key", pd.Series(dtype=str)) == key] if "pair_key" in raw.columns else pd.DataFrame()
        clean_full_hit = clean_full[clean_full["pair_key"] == key]
        clean_hit = clean[clean["pair_key"] == key]
        pairs_hit = pairs[pairs["pair_key"] == key]
        challenge_hit = challenge[name_match_mask(challenge, item.substrate_name, item.product_name)] if not challenge.empty else pd.DataFrame()

        if not rxn_hit.empty:
            rxn = rxn_hit.iloc[0]
            substrate_smiles = rxn["substrate_smiles"]
            product_smiles = rxn["product_smiles"]
            module = rxn["module"]
            tier = rxn.get("tier", "")
            is_clean = bool(rxn.get("is_clean", False))
        else:
            substrate_smiles = ""
            product_smiles = ""
            module = item.module
            tier = ""
            is_clean = False

        for role, name, smiles, expected_b1 in [
            ("substrate", item.substrate_name, substrate_smiles, key.split("__")[0]),
            ("product", item.product_name, product_smiles, key.split("__")[1]),
        ]:
            local_ik = smiles_to_inchikey(smiles)
            pub = load_pubchem(name)
            pub_ik = str(pub.get("InChIKey", ""))
            mapping_rows.append(
                {
                    "round": 9,
                    "pair_key": key,
                    "compound_role": role,
                    "compound_name": name,
                    "local_smiles": smiles,
                    "local_canonical_smiles": smiles_to_canonical(smiles, isomeric=True),
                    "local_inchikey": local_ik,
                    "local_block1": b1(local_ik),
                    "expected_block1": expected_b1,
                    "pubchem_cid": pub.get("CID", ""),
                    "pubchem_inchikey": pub_ik,
                    "pubchem_block1": b1(pub_ik),
                    "pubchem_isomeric_smiles": pub.get("IsomericSMILES", pub.get("SMILES", "")),
                    "pubchem_canonical_smiles": pub.get("CanonicalSMILES", pub.get("ConnectivitySMILES", "")),
                    "pubchem_formula": pub.get("MolecularFormula", ""),
                    "pubchem_weight": pub.get("MolecularWeight", ""),
                    "block1_matches_pair_key": b1(local_ik) == expected_b1,
                    "local_pubchem_block1_match": bool(pub_ik) and b1(local_ik) == b1(pub_ik),
                    "mapping_status_round9": (
                        "confirmed_by_local_and_pubchem_block1"
                        if bool(pub_ik) and b1(local_ik) == b1(pub_ik) == expected_b1
                        else "needs_manual_structure_review"
                    ),
                }
            )

        generation_status = (
            "reaches_clean_candidates"
            if not clean_hit.empty
            else "reaches_clean_candidates_full_only"
            if not clean_full_hit.empty
            else "positive_in_reactions_but_not_generated_in_clean_candidates"
            if not rxn_hit.empty
            else "not_in_reactions"
        )
        import_decision = (
            "do_not_import_duplicate_positive; upgrade_evidence_and_fix_generator"
            if not rxn_hit.empty and clean_full_hit.empty
            else "do_not_import_duplicate_positive"
            if not rxn_hit.empty
            else "candidate_import_after_mapping_and_split_guard"
        )
        split_risk = (
            "holdout_or_challenge_overlap"
            if not challenge_hit.empty
            else "existing_gold_positive_requires_split_guard"
            if not rxn_hit.empty
            else "new_candidate_requires_split_assignment"
        )

        import_rows.append(
            {
                "round": 9,
                "pair_key": key,
                "module": module,
                "substrate_name": item.substrate_name,
                "product_name": item.product_name,
                "pmid": item.pmid,
                "doi": item.doi,
                "round8_evidence_class": item.abstract_evidence_class,
                "raw_pool_hit_count": len(raw_hit),
                "reactions_hit_count": len(rxn_hit),
                "reaction_tier": tier,
                "reaction_is_clean": is_clean,
                "clean_candidates_full_hit_count": len(clean_full_hit),
                "clean_candidates_hit_count": len(clean_hit),
                "challenge_name_hit_count": len(challenge_hit),
                "generation_status_round9": generation_status,
                "round9_import_decision": import_decision,
                "training_allowed_round9": False,
                "why_training_still_blocked": (
                    "All four Round8 exact candidates already exist as gold reactions but are absent from clean candidate outputs; "
                    "the next improvement is candidate-generation recall/rule provenance, not duplicate positive import."
                ),
            }
        )
        split_rows.append(
            {
                "round": 9,
                "pair_key": key,
                "substrate_name": item.substrate_name,
                "product_name": item.product_name,
                "split_or_holdout_risk": split_risk,
                "challenge_case_ids": ";".join(challenge_hit.get("case_id", pd.Series(dtype=str)).astype(str).tolist())
                if not challenge_hit.empty and "case_id" in challenge_hit.columns
                else "",
                "scaffold_fold_values_in_pairs": ";".join(
                    sorted(pairs_hit.get("scaffold_fold", pd.Series(dtype=str)).dropna().astype(str).unique().tolist())
                ),
                "candidate_generation_status": generation_status,
                "recommended_action": (
                    "keep out of training until challenge/holdout policy is explicit"
                    if not challenge_hit.empty
                    else "trace why clean candidate generation misses this gold reaction"
                ),
            }
        )
        queue_rows.append(
            {
                "round": 9,
                "queue_type": "candidate_generation_recall_fix",
                "priority": "high",
                "pair_key": key,
                "substrate_name": item.substrate_name,
                "product_name": item.product_name,
                "next_action": "recover or import source-traceable rule so deployment generator can enumerate the known gold product",
            }
        )

    neg = pd.read_csv(ROUND8_NEG)
    neg_out = neg.copy()
    if not neg_out.empty:
        neg_out["source_round"] = neg_out["round"]
        neg_out["round"] = 9
        neg_out["round9_decision"] = "keep_as_conditional_negative_candidate_not_global_training_negative"
        neg_out["training_allowed_round9"] = False
        neg_out["required_before_any_use"] = (
            "exact tested substrate/product-or-product-class, enzyme/organism, assay condition, detection statement, and scope boundary"
        )

    pd.DataFrame(mapping_rows).to_csv(OUT_MAPPING, index=False)
    pd.DataFrame(import_rows).to_csv(OUT_IMPORT, index=False)
    pd.DataFrame(split_rows).to_csv(OUT_SPLIT, index=False)
    neg_out.to_csv(OUT_NEG, index=False)
    pd.DataFrame(queue_rows).to_csv(OUT_QUEUE, index=False)

    print(f"wrote {OUT_MAPPING} rows={len(mapping_rows)}")
    print(pd.DataFrame(mapping_rows)["mapping_status_round9"].value_counts(dropna=False).to_string())
    print(f"wrote {OUT_IMPORT} rows={len(import_rows)}")
    print(pd.DataFrame(import_rows)["generation_status_round9"].value_counts(dropna=False).to_string())
    print(f"wrote {OUT_SPLIT} rows={len(split_rows)}")
    print(f"wrote {OUT_NEG} rows={len(neg_out)}")
    print(f"wrote {OUT_QUEUE} rows={len(queue_rows)}")


if __name__ == "__main__":
    main()
