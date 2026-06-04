#!/usr/bin/env python
"""Round15 structure mapping and generator recall audit.

Round15 starts from the source-backed Rhea candidates found in Round14. It
checks whether they are truly new training samples, whether ChEBI/PubChem/local
structures agree at InChIKey block1, and whether the deployment-like generator
can enumerate the target product.

This script is audit-only: it does not modify model training data or generator
code.
"""

from __future__ import annotations

import html
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger


REPO = Path(__file__).resolve().parents[1]
SRC = REPO / "src"
MODULAR = REPO / "modular"
sys.path.insert(0, str(SRC))
sys.path.insert(0, str(MODULAR))

import kio  # noqa: E402
from generate_candidates import run_reactants  # noqa: E402
from ltr_candidates import MACRO2MOD, POOL, ec_class, load_module_rules  # noqa: E402


RDLogger.DisableLog("rdApp.*")

CURATION = REPO / "data" / "curation"
DOCS = REPO / "docs" / "reviews"
CHEBI_RAW = CURATION / "chebi_round15_raw"
PUBCHEM_RAW = CURATION / "pubchem_round15_raw"

RAW_POSITIVE_POOL = Path(POOL)
REACTIONS = REPO / "outputs" / "modular" / "ltr" / "reactions.parquet"
CLEAN_FULL = REPO / "outputs" / "modular" / "ltr" / "clean_candidates_full.parquet"
CLEAN = REPO / "outputs" / "modular" / "ltr" / "clean_candidates.parquet"
FULLRULE = REPO / "outputs" / "gen_gap" / "fullrule_hit_miss.parquet"

OUT_COMPOUNDS = CURATION / "round15_compound_structure_mapping.csv"
OUT_PAIRS = CURATION / "round15_candidate_pair_gate.csv"
OUT_GENERATOR = CURATION / "round15_generator_recall_check.csv"
OUT_DECISIONS = CURATION / "round15_training_decisions.csv"
OUT_NEXT = CURATION / "round15_next_actions.csv"
OUT_DOC = DOCS / "structure_and_generator_round15.md"


CANDIDATE_REACTIONS = [
    {
        "source_id": "RHEA:16309",
        "priority_family": "bile_acid_deconjugation_and_lipid_context",
        "proposed_module": "D",
        "equation": "taurochenodeoxycholate + H2O = chenodeoxycholate + taurine",
        "substrate_chebi": "9407",
        "primary_product_chebi": "36234",
        "coproduct_chebis": ["507393"],
    },
    {
        "source_id": "RHEA:19353",
        "priority_family": "bile_acid_deconjugation_and_lipid_context",
        "proposed_module": "D",
        "equation": "glycocholate + H2O = cholate + glycine",
        "substrate_chebi": "29746",
        "primary_product_chebi": "29747",
        "coproduct_chebis": ["57305"],
    },
    {
        "source_id": "RHEA:20689",
        "priority_family": "hydroxycinnamate_reduction_and_hydrolysis",
        "proposed_module": "A",
        "equation": "chlorogenate + H2O = L-quinate + (E)-caffeate + H(+)",
        "substrate_chebi": "57644",
        "primary_product_chebi": "57770",
        "coproduct_chebis": ["29751"],
    },
    {
        "source_id": "RHEA:69683",
        "priority_family": "isoflavone_reductive_and_glycoside_metabolism",
        "proposed_module": "B",
        "equation": "daidzein 7-O-beta-D-glucoside + H2O = daidzein + beta-D-glucose + H(+)",
        "substrate_chebi": "42202",
        "primary_product_chebi": "77764",
        "coproduct_chebis": ["15903"],
    },
]


PUBCHEM_BY_CHEBI = {
    "9407": "taurochenodeoxycholate",
    "36234": "chenodeoxycholate",
    "29746": "glycocholate",
    "29747": "cholate",
    "57644": "chlorogenate",
    "57770": "trans_caffeate",
    "42202": "daidzein_7_o_beta_d_glucoside",
    "77764": "daidzein",
}


def b1(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    return str(value).split("-", 1)[0][:14]


def clean_html(value: object) -> str:
    text = "" if value is None else str(value)
    text = re.sub(r"<[^>]+>", "", text)
    return html.unescape(text).strip()


def read_json(path: Path) -> dict:
    if not path.exists():
        return {}
    return json.loads(path.read_text(encoding="utf-8"))


def chebi_record(chebi_id: str) -> dict[str, object]:
    raw_path = CHEBI_RAW / f"CHEBI_{chebi_id}.json"
    raw = read_json(raw_path)
    rec = raw.get("summary", raw)
    structure = rec.get("default_structure") or {}
    chem = rec.get("chemical_data") or {}
    smiles = structure.get("smiles", "")
    inchikey = structure.get("standard_inchi_key", "")
    return {
        "chebi_id": f"CHEBI:{chebi_id}",
        "chebi_name": clean_html(rec.get("name", "")),
        "chebi_ascii_name": clean_html(rec.get("ascii_name", "")),
        "chebi_formula": chem.get("formula", ""),
        "chebi_charge": chem.get("charge", ""),
        "chebi_smiles": smiles,
        "chebi_inchikey": inchikey,
        "chebi_block1": b1(inchikey),
        "chebi_structure_status": "mapped" if smiles and inchikey else "missing_structure",
        "chebi_raw_file": str(raw_path.relative_to(REPO)) if raw_path.exists() else "",
    }


def pubchem_record(chebi_id: str) -> dict[str, object]:
    stem = PUBCHEM_BY_CHEBI.get(chebi_id)
    if not stem:
        return {
            "pubchem_query_id": "",
            "pubchem_cid": "",
            "pubchem_formula": "",
            "pubchem_smiles": "",
            "pubchem_inchikey": "",
            "pubchem_block1": "",
            "pubchem_structure_status": "not_queried",
            "pubchem_raw_file": "",
        }
    raw_path = PUBCHEM_RAW / f"{stem}_pubchem.json"
    raw = read_json(raw_path)
    records = raw.get("PropertyTable", {}).get("Properties", [])
    rec = records[0] if records else {}
    smiles = rec.get("SMILES") or rec.get("IsomericSMILES") or rec.get("CanonicalSMILES") or ""
    inchikey = rec.get("InChIKey", "")
    return {
        "pubchem_query_id": stem,
        "pubchem_cid": rec.get("CID", ""),
        "pubchem_formula": rec.get("MolecularFormula", ""),
        "pubchem_smiles": smiles,
        "pubchem_inchikey": inchikey,
        "pubchem_block1": b1(inchikey),
        "pubchem_structure_status": "mapped" if smiles and inchikey else "missing_or_not_queried",
        "pubchem_raw_file": str(raw_path.relative_to(REPO)) if raw_path.exists() else "",
    }


def build_compound_mapping() -> pd.DataFrame:
    rows = []
    role_rows = []
    for rxn in CANDIDATE_REACTIONS:
        role_rows.append((rxn["source_id"], rxn["substrate_chebi"], "substrate"))
        role_rows.append((rxn["source_id"], rxn["primary_product_chebi"], "primary_product"))
        for chebi_id in rxn["coproduct_chebis"]:
            role_rows.append((rxn["source_id"], chebi_id, "coproduct_context"))

    for source_id, chebi_id, role in role_rows:
        c = chebi_record(chebi_id)
        p = pubchem_record(chebi_id)
        if p["pubchem_structure_status"] == "mapped":
            status = "consistent_block1" if c["chebi_block1"] == p["pubchem_block1"] else "block1_conflict"
        else:
            status = "chebi_only"
        rows.append(
            {
                "round": 15,
                "source_id": source_id,
                "compound_role": role,
                **c,
                **p,
                "cross_source_structure_status": status,
                "training_mapping_note": "Use ChEBI as primary source; PubChem is a block1 cross-check when available.",
            }
        )
    return pd.DataFrame(rows)


def build_pair_seed(compounds: pd.DataFrame) -> pd.DataFrame:
    by_source_role = {
        (r.source_id, r.compound_role): r
        for r in compounds.itertuples(index=False)
        if r.compound_role in {"substrate", "primary_product"}
    }
    rows = []
    for spec in CANDIDATE_REACTIONS:
        sub = by_source_role[(spec["source_id"], "substrate")]
        prod = by_source_role[(spec["source_id"], "primary_product")]
        rows.append(
            {
                **spec,
                "substrate_name": sub.chebi_ascii_name or sub.chebi_name,
                "product_name": prod.chebi_ascii_name or prod.chebi_name,
                "sb": sub.chebi_block1,
                "pb": prod.chebi_block1,
                "pair_key": f"{sub.chebi_block1}__{prod.chebi_block1}",
                "substrate_smiles_chebi": sub.chebi_smiles,
                "product_smiles_chebi": prod.chebi_smiles,
                "substrate_structure_status": sub.cross_source_structure_status,
                "product_structure_status": prod.cross_source_structure_status,
            }
        )
    return pd.DataFrame(rows)


def load_raw_pool_pair_counts() -> pd.DataFrame:
    cols = ["substrate_inchikey", "product_inchikey", "source_dataset", "source_origin_type", "training_use_recommendation"]
    pool = pd.read_csv(RAW_POSITIVE_POOL, usecols=cols, low_memory=False)
    pool["sb"] = pool["substrate_inchikey"].map(b1)
    pool["pb"] = pool["product_inchikey"].map(b1)
    pool = pool[pool["sb"].ne("") & pool["pb"].ne("")]
    return (
        pool.groupby(["sb", "pb"], dropna=False)
        .agg(
            raw_positive_pool_pair_rows=("sb", "size"),
            raw_source_datasets=("source_dataset", lambda x: ";".join(sorted(set(map(str, x))))),
            raw_source_origin_types=("source_origin_type", lambda x: ";".join(sorted(set(map(str, x))))),
            raw_training_recommendations=("training_use_recommendation", lambda x: ";".join(sorted(set(map(str, x))))),
        )
        .reset_index()
    )


def table_pair_counts(path: Path, prefix: str) -> pd.DataFrame:
    df = pd.read_parquet(path)
    aggs: dict[str, tuple[str, object]] = {f"{prefix}_pair_rows": ("sb", "size")}
    if "module" in df.columns:
        aggs[f"{prefix}_modules"] = ("module", lambda x: ";".join(sorted(set(map(str, x)))))
    if "tier" in df.columns:
        aggs[f"{prefix}_tiers"] = ("tier", lambda x: ";".join(sorted(set(map(str, x)))))
    if "is_clean" in df.columns:
        aggs[f"{prefix}_any_clean"] = ("is_clean", "max")
    if "y" in df.columns:
        aggs[f"{prefix}_positive_rows"] = ("y", "sum")
    return df.groupby(["sb", "pb"], dropna=False).agg(**aggs).reset_index()


def first_reaction_row(reactions: pd.DataFrame, sb: str, pb: str) -> pd.Series | None:
    hit = reactions[reactions["sb"].eq(sb) & reactions["pb"].eq(pb)]
    if hit.empty:
        return None
    return hit.sort_values(["is_clean", "tier"], ascending=[False, True]).iloc[0]


def build_pair_gate(compounds: pd.DataFrame) -> pd.DataFrame:
    pairs = build_pair_seed(compounds)
    for extra in [
        load_raw_pool_pair_counts(),
        table_pair_counts(REACTIONS, "reactions"),
        table_pair_counts(CLEAN_FULL, "clean_full"),
        table_pair_counts(CLEAN, "clean"),
    ]:
        pairs = pairs.merge(extra, on=["sb", "pb"], how="left")

    reactions = pd.read_parquet(REACTIONS)
    modules = []
    generator_smiles = []
    generator_smiles_source = []
    for row in pairs.itertuples(index=False):
        rxn_row = first_reaction_row(reactions, row.sb, row.pb)
        if rxn_row is not None:
            modules.append(rxn_row["module"])
            generator_smiles.append(rxn_row["substrate_smiles"])
            generator_smiles_source.append("local_reactions")
        else:
            modules.append(row.proposed_module)
            generator_smiles.append(row.substrate_smiles_chebi)
            generator_smiles_source.append("chebi_source")
    pairs["module"] = modules
    pairs["generator_input_substrate_smiles"] = generator_smiles
    pairs["generator_input_smiles_source"] = generator_smiles_source

    for col in [
        "raw_positive_pool_pair_rows",
        "reactions_pair_rows",
        "clean_full_pair_rows",
        "clean_pair_rows",
        "clean_full_positive_rows",
        "clean_positive_rows",
    ]:
        if col in pairs.columns:
            pairs[col] = pd.to_numeric(pairs[col], errors="coerce").fillna(0).astype(int)
    for col in pairs.columns:
        if pairs[col].dtype == object:
            pairs[col] = pairs[col].fillna("")

    pairs["already_in_raw_positive_pool"] = pairs["raw_positive_pool_pair_rows"].gt(0)
    pairs["already_in_reactions"] = pairs["reactions_pair_rows"].gt(0)
    pairs["already_in_clean_candidates_full"] = pairs["clean_full_pair_rows"].gt(0)
    pairs["already_in_clean_candidates"] = pairs["clean_pair_rows"].gt(0)
    pairs["structure_mapping_pass"] = (
        pairs["substrate_structure_status"].isin(["consistent_block1", "chebi_only"])
        & pairs["product_structure_status"].isin(["consistent_block1", "chebi_only"])
        & pairs["sb"].ne("")
        & pairs["pb"].ne("")
    )
    return pairs


def load_all_positive_pool_rules() -> dict[str, list[tuple[str, int]]]:
    d = pd.read_csv(POOL, low_memory=False, usecols=["macro_module", "reaction_rule_smarts", "enzyme_ec"])
    d["module"] = d["macro_module"].map(MACRO2MOD)
    d = d.dropna(subset=["module", "reaction_rule_smarts"])
    d = d[d["reaction_rule_smarts"].astype(str).str.contains(">>", regex=False)]
    d["ecc"] = d["enzyme_ec"].map(ec_class)
    out: dict[str, list[tuple[str, int]]] = {}
    for module, group in d.groupby("module"):
        unique = group.drop_duplicates("reaction_rule_smarts")[["reaction_rule_smarts", "ecc"]]
        out[module] = list(unique.itertuples(index=False, name=None))
    return out


def generator_hit_for_rules(substrate_smiles: str, substrate_block: str, target_block: str, rules: list[tuple[str, int]]) -> dict[str, object]:
    mol = Chem.MolFromSmiles(substrate_smiles)
    if mol is None:
        return {"target_generated": False, "generated_product_count": 0, "target_rule_count": 0, "target_product_smiles_examples": ""}
    molh = Chem.AddHs(mol)
    generated = set()
    target_rule_count = 0
    target_examples = []
    for smarts, _ecc in rules:
        for product_smiles in run_reactants(str(smarts), mol, molh):
            generated_block = b1(kio.smiles_to_inchikey(product_smiles))
            if not generated_block or generated_block == substrate_block:
                continue
            generated.add(generated_block)
            if generated_block == target_block:
                target_rule_count += 1
                target_examples.append(product_smiles)
    return {
        "target_generated": target_block in generated,
        "generated_product_count": len(generated),
        "target_rule_count": target_rule_count,
        "target_product_smiles_examples": ";".join(target_examples[:3]),
    }


def fullrule_lookup() -> pd.DataFrame:
    df = pd.read_parquet(FULLRULE)
    return (
        df.groupby(["sb", "pb"], dropna=False)
        .agg(
            fullrule_existing_rows=("hit", "size"),
            fullrule_existing_hit_rows=("hit", "sum"),
            fullrule_existing_any_hit=("hit", "max"),
        )
        .reset_index()
    )


def build_generator_recall(pairs: pd.DataFrame) -> pd.DataFrame:
    sampled_rules = load_module_rules(np.random.default_rng(0), 800)
    full_pool_rules = load_all_positive_pool_rules()
    fullrule = fullrule_lookup()
    rows = []
    for item in pairs.itertuples(index=False):
        module = item.module
        sampled = generator_hit_for_rules(
            item.generator_input_substrate_smiles,
            item.sb,
            item.pb,
            sampled_rules.get(module, []),
        )
        full_pool = generator_hit_for_rules(
            item.generator_input_substrate_smiles,
            item.sb,
            item.pb,
            full_pool_rules.get(module, []),
        )
        rows.append(
            {
                "round": 15,
                "source_id": item.source_id,
                "pair_key": item.pair_key,
                "sb": item.sb,
                "pb": item.pb,
                "module": module,
                "substrate_name": item.substrate_name,
                "product_name": item.product_name,
                "generator_input_smiles_source": item.generator_input_smiles_source,
                "sampled_clean_rules_n": len(sampled_rules.get(module, [])),
                "sampled_clean_rules_target_generated": sampled["target_generated"],
                "sampled_clean_rules_generated_product_count": sampled["generated_product_count"],
                "sampled_clean_rules_target_rule_count": sampled["target_rule_count"],
                "sampled_clean_rules_target_smiles_examples": sampled["target_product_smiles_examples"],
                "full_positive_pool_rules_n": len(full_pool_rules.get(module, [])),
                "full_positive_pool_rules_target_generated": full_pool["target_generated"],
                "full_positive_pool_rules_generated_product_count": full_pool["generated_product_count"],
                "full_positive_pool_rules_target_rule_count": full_pool["target_rule_count"],
                "full_positive_pool_rules_target_smiles_examples": full_pool["target_product_smiles_examples"],
            }
        )
    out = pd.DataFrame(rows)
    out = out.merge(fullrule, on=["sb", "pb"], how="left")
    for col in ["fullrule_existing_rows", "fullrule_existing_hit_rows", "fullrule_existing_any_hit"]:
        out[col] = pd.to_numeric(out[col], errors="coerce").fillna(0).astype(int)
    return out


def decision_for(pair: pd.Series, gen: pd.Series) -> tuple[str, bool, str, str]:
    if pair["already_in_reactions"] and pair["already_in_clean_candidates"]:
        return (
            "do_not_import_existing_positive_already_in_clean_candidates",
            False,
            "Existing positive is already available to the ranker; use Rhea as provenance enrichment only.",
            "provenance_enrichment_or_source_manifest_only",
        )
    if pair["already_in_reactions"] and not pair["already_in_clean_candidates"]:
        if pair["already_in_clean_candidates_full"]:
            return (
                "do_not_import_existing_positive_dropped_after_clean_full",
                False,
                "Existing positive reaches clean_candidates_full but not final clean_candidates; inspect filtering/sampling/evaluation path before retraining.",
                "trace_clean_full_to_clean_candidates_loss",
            )
        if bool(gen["fullrule_existing_any_hit"]) or bool(gen["full_positive_pool_rules_target_generated"]):
            return (
                "do_not_import_existing_positive_fix_generator_path",
                False,
                "Existing positive has source/database support and full-rule evidence, but clean candidates miss it; repair generator/rule path, not duplicate label.",
                "promote_or_wire_source_traceable_rule_then_rerun_clean_candidates",
            )
        return (
            "do_not_import_existing_positive_requires_rule_recovery",
            False,
            "Existing positive is absent from clean candidates and checked generator paths did not recover it.",
            "manual_rule_recovery_or_source_rule_import",
        )
    if not pair["structure_mapping_pass"]:
        return (
            "blocked_structure_mapping",
            False,
            "ChEBI/PubChem/local structure mapping is not stable enough for training.",
            "manual_structure_review",
        )
    if bool(gen["sampled_clean_rules_target_generated"]):
        return (
            "candidate_new_positive_after_full_curation",
            False,
            "The deployment-like generator can enumerate the target, but import still needs manual context/license/split review.",
            "build_import_manifest_after_manual_curation",
        )
    return (
        "blocked_generator_cannot_enumerate_target",
        False,
        "Source-backed pair is not useful for retraining until generator recall is repaired.",
        "recover_source_traceable_rule_before_import",
    )


def build_decisions(pairs: pd.DataFrame, generator: pd.DataFrame) -> pd.DataFrame:
    gen_by_key = {r.pair_key: r for r in generator.itertuples(index=False)}
    rows = []
    for pair in pairs.itertuples(index=False):
        gen = gen_by_key[pair.pair_key]
        status, allowed, reason, next_action = decision_for(pd.Series(pair._asdict()), pd.Series(gen._asdict()))
        rows.append(
            {
                "round": 15,
                "source_id": pair.source_id,
                "pair_key": pair.pair_key,
                "module": pair.module,
                "substrate_name": pair.substrate_name,
                "product_name": pair.product_name,
                "training_decision_round15": status,
                "training_allowed_round15": allowed,
                "reason": reason,
                "next_action": next_action,
            }
        )
    return pd.DataFrame(rows)


def build_next_actions(decisions: pd.DataFrame) -> pd.DataFrame:
    priority = {
        "trace_clean_full_to_clean_candidates_loss": 1,
        "promote_or_wire_source_traceable_rule_then_rerun_clean_candidates": 2,
        "manual_rule_recovery_or_source_rule_import": 3,
        "provenance_enrichment_or_source_manifest_only": 4,
        "build_import_manifest_after_manual_curation": 5,
        "recover_source_traceable_rule_before_import": 6,
        "manual_structure_review": 7,
    }
    rows = []
    for action, group in decisions.groupby("next_action", dropna=False):
        rows.append(
            {
                "round": 15,
                "priority": priority.get(str(action), 99),
                "next_action": action,
                "n_pairs": len(group),
                "pair_keys": ";".join(group["pair_key"].astype(str)),
                "why": "; ".join(group["reason"].astype(str).head(2)),
            }
        )
    return pd.DataFrame(rows).sort_values(["priority", "next_action"])


def write_doc(compounds: pd.DataFrame, pairs: pd.DataFrame, generator: pd.DataFrame, decisions: pd.DataFrame, next_actions: pd.DataFrame) -> None:
    decision_counts = decisions["training_decision_round15"].value_counts().to_dict()
    clean_counts = {
        "already_in_reactions": int(pairs["already_in_reactions"].sum()),
        "already_in_clean_candidates_full": int(pairs["already_in_clean_candidates_full"].sum()),
        "already_in_clean_candidates": int(pairs["already_in_clean_candidates"].sum()),
    }
    gen_rows = "\n".join(
        f"- `{r.source_id}` `{r.substrate_name} -> {r.product_name}`: sampled={r.sampled_clean_rules_target_generated}, full_pool={r.full_positive_pool_rules_target_generated}, fullrule_hit={bool(r.fullrule_existing_any_hit)}"
        for r in generator.itertuples(index=False)
    )
    decision_rows = "\n".join(
        f"- `{r.source_id}` `{r.substrate_name} -> {r.product_name}`: `{r.training_decision_round15}`"
        for r in decisions.itertuples(index=False)
    )
    action_rows = "\n".join(
        f"- P{r.priority} `{r.next_action}` ({r.n_pairs} pair): {r.pair_keys}"
        for r in next_actions.itertuples(index=False)
    )
    text = f"""# Round15 Structure And Generator Audit

Date: 2026-06-04

## Direct Answer

Round15 makes the diagnosis sharper:

> The four Rhea source-backed candidates from Round14 are not simply new samples to add. All four already exist in `reactions.parquet`. The production question is whether they survive the clean generator path.

Local availability counts: `{clean_counts}`

Decision counts: `{decision_counts}`

No Round15 row is allowed directly into training.

## Structure Mapping

ChEBI was used as the primary structure source and PubChem as an InChIKey block cross-check. Primary substrate/product mappings are consistent at block1 level where PubChem was queried.

Structure table:

- `data/curation/round15_compound_structure_mapping.csv`

## Generator Recall

{gen_rows}

Interpretation:

- If a pair is in `reactions.parquet` but absent from `clean_candidates`, duplicate-importing it as a new positive is wrong.
- If full or positive-pool rules can generate it but clean candidates do not contain it, the next action is generator/rule wiring or clean-path trace.
- If the pair reaches `clean_candidates_full` but not `clean_candidates`, the loss is a filtering/sampling/evaluation-path issue.

## Training Decisions

{decision_rows}

## Next Actions

{action_rows}

## Files Created

- `data/curation/round15_compound_structure_mapping.csv`
- `data/curation/round15_candidate_pair_gate.csv`
- `data/curation/round15_generator_recall_check.csv`
- `data/curation/round15_training_decisions.csv`
- `data/curation/round15_next_actions.csv`
- `data/curation/chebi_round15_raw/*.json`
- `data/curation/pubchem_round15_raw/*.json`

## Claim Boundary

Allowed after Round15:

> Round15 confirmed that the strongest Rhea candidates are already represented in the normalized reaction table, and that production progress depends on generator recall/path retention rather than duplicate positive import.

Forbidden:

> Round15 expanded the training set or proved production readiness.
"""
    OUT_DOC.write_text(text.strip() + "\n", encoding="utf-8")


def main() -> None:
    CURATION.mkdir(parents=True, exist_ok=True)
    DOCS.mkdir(parents=True, exist_ok=True)
    compounds = build_compound_mapping()
    pairs = build_pair_gate(compounds)
    generator = build_generator_recall(pairs)
    decisions = build_decisions(pairs, generator)
    next_actions = build_next_actions(decisions)
    for path, table in [
        (OUT_COMPOUNDS, compounds),
        (OUT_PAIRS, pairs),
        (OUT_GENERATOR, generator),
        (OUT_DECISIONS, decisions),
        (OUT_NEXT, next_actions),
    ]:
        table.to_csv(path, index=False)
        print(f"wrote {path} rows={len(table)}")
    write_doc(compounds, pairs, generator, decisions, next_actions)
    print(f"wrote {OUT_DOC}")
    print("round15_decisions=", decisions["training_decision_round15"].value_counts().to_dict())
    print("round15_training_allowed=", int(decisions["training_allowed_round15"].sum()))


if __name__ == "__main__":
    main()
