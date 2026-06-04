#!/usr/bin/env python
"""Round13 dry-run for source-traceable generator overlay candidates."""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
from rdkit import Chem
from rdkit import RDLogger


REPO = Path(__file__).resolve().parents[1]
SRC = REPO / "src"
sys.path.insert(0, str(SRC))

import kio  # noqa: E402
from generate_candidates import run_reactants  # noqa: E402


RDLogger.DisableLog("rdApp.*")

CURATION = REPO / "data" / "curation"
ROUND12_OVERLAY = CURATION / "rule_promotion_overlay_round12.csv"
ROUND5_RECOVERY = CURATION / "master_rule_recovery_candidates_round5.csv"
REACTIONS = REPO / "outputs" / "modular" / "ltr" / "reactions.parquet"

OUT_DRYRUN = CURATION / "rule_overlay_dryrun_round13.csv"
OUT_DECISION = CURATION / "round13_generator_repair_decision.csv"


def b1(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    return str(value).split("-", 1)[0]


def main() -> None:
    overlay = pd.read_csv(ROUND12_OVERLAY)
    recovery = pd.read_csv(ROUND5_RECOVERY)
    recovery["pair_key"] = recovery["sb"].astype(str) + "__" + recovery["pb"].astype(str)
    reactions = pd.read_parquet(REACTIONS)
    reactions["pair_key"] = reactions["sb"].astype(str) + "__" + reactions["pb"].astype(str)

    allowed = overlay[overlay["overlay_dryrun_allowed_round12"].astype(bool)].copy()
    rows = []
    decisions = []
    for item in allowed.itertuples(index=False):
        rule_rows = recovery[
            recovery["pair_key"].eq(item.pair_key)
            & recovery["rule_id"].eq(item.rule_id)
        ]
        reaction_rows = reactions[reactions["pair_key"].eq(item.pair_key)]
        if rule_rows.empty or reaction_rows.empty:
            rows.append(
                {
                    "round": 13,
                    "pair_key": item.pair_key,
                    "rule_id": item.rule_id,
                    "dryrun_status": "blocked_missing_rule_or_reaction_row",
                    "target_generated": False,
                    "generated_product_count": 0,
                }
            )
            continue
        rule = rule_rows.iloc[0]
        rxn = reaction_rows.iloc[0]
        mol = Chem.MolFromSmiles(rxn["substrate_smiles"])
        molh = Chem.AddHs(mol) if mol is not None else None
        generated = []
        for product_smiles in run_reactants(rule["rule_smarts"], mol, molh):
            product_inchikey = kio.smiles_to_inchikey(product_smiles)
            product_block = b1(product_inchikey)
            if product_block and product_block != str(rxn["sb"]):
                generated.append(
                    {
                        "product_smiles": product_smiles,
                        "product_inchikey": product_inchikey,
                        "product_block": product_block,
                    }
                )
        target_block = str(rxn["pb"])
        target_hits = [g for g in generated if g["product_block"] == target_block]
        rows.append(
            {
                "round": 13,
                "pair_key": item.pair_key,
                "module": item.module,
                "substrate_name": item.substrate_name,
                "product_name": item.product_name,
                "target_product_block": target_block,
                "rule_id": item.rule_id,
                "rule_source": item.rule_source,
                "rule_legacy_id": item.rule_legacy_id,
                "source_reaction_check_status": item.source_reaction_check_status,
                "dryrun_status": "target_generated_by_overlay_rule" if target_hits else "target_not_generated_by_overlay_rule",
                "target_generated": bool(target_hits),
                "generated_product_count": len({g["product_block"] for g in generated}),
                "target_product_smiles_from_rule": target_hits[0]["product_smiles"] if target_hits else "",
                "target_product_inchikey_from_rule": target_hits[0]["product_inchikey"] if target_hits else "",
                "training_label_action": "no_label_change_existing_gold_positive",
                "deployment_action": "dryrun_only_not_deployment_promoted",
            }
        )
        decisions.append(
            {
                "round": 13,
                "pair_key": item.pair_key,
                "rule_id": item.rule_id,
                "generator_repair_result": "overlay_rule_recovers_target_candidate" if target_hits else "overlay_rule_failed_to_recover_target",
                "can_modify_generator_next": bool(target_hits),
                "training_allowed_round13": False,
                "deployment_promotion_allowed_round13": False,
                "next_required_gate": "wire overlay into clean generator for target-family test, then rerun clean candidate output and leakage checks",
            }
        )

    dryrun = pd.DataFrame(rows)
    decision = pd.DataFrame(decisions)
    dryrun.to_csv(OUT_DRYRUN, index=False)
    decision.to_csv(OUT_DECISION, index=False)
    print(f"wrote {OUT_DRYRUN} rows={len(dryrun)}")
    print(f"wrote {OUT_DECISION} rows={len(decision)}")
    if not dryrun.empty:
        print("round13_dryrun_status=", dryrun["dryrun_status"].value_counts().to_dict())


if __name__ == "__main__":
    main()

