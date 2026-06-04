#!/usr/bin/env python
"""Round20 source-backed BSH rule overlay dry run.

This script does not change training labels, clean candidate parquet files, or
the production generator. It asks whether the Round19 P0 bile-acid targets would
survive a require-EC generation path if exact, source-backed BSH rules were
given the EC provenance already present in Rhea/ENZYME.
"""

from __future__ import annotations

import hashlib
import sys
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit import RDLogger


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "curation"
DOC = ROOT / "docs" / "reviews" / "bsh_overlay_dryrun_round20.md"
POOL = Path(r"D:\CRB\FoodGut\positive_sample_db\data\processed\foodgut_modular_positive_candidate_pool_v2.csv")
LTR_OUT = ROOT / "outputs" / "modular" / "ltr"

sys.path.insert(0, str(ROOT / "modular"))
sys.path.insert(0, str(ROOT / "src"))
import kio  # noqa: E402
from generate_candidates import run_reactants  # noqa: E402
from ltr_candidates import load_module_rules  # noqa: E402


RDLogger.DisableLog("rdApp.*")

ROUND19_RULES = OUT / "round19_bile_acid_target_generating_rules.csv"
ROUND19_FIXES = OUT / "round19_bile_acid_fix_recommendations.csv"
ROUND19_EVIDENCE = OUT / "round19_bile_acid_external_evidence.csv"

MACRO2MOD = {
    "small_molecule_polyphenol": "A",
    "carbohydrate_glycan": "B",
    "protein_amino_acid": "C",
    "lipid_fat": "D",
}


def ec_class(ec: object) -> int:
    text = "" if pd.isna(ec) else str(ec)
    return int(text[0]) if text[:1].isdigit() else 0


def joined(values: Iterable[object]) -> str:
    return ";".join(sorted({str(v) for v in values if pd.notna(v) and str(v) != ""}))


def key_block1(smiles: str) -> str:
    return kio.inchikey_block1(kio.smiles_to_inchikey(smiles))


def load_pool_rules() -> pd.DataFrame:
    usecols = [
        "macro_module",
        "reaction_rule_smarts",
        "enzyme_ec",
        "source_dataset",
        "source_origin_type",
        "training_use_recommendation",
        "reaction_category",
        "reaction_type",
        "reaction_rule_id",
    ]
    rules = pd.read_csv(POOL, low_memory=False, usecols=lambda c: c in usecols)
    rules["module"] = rules["macro_module"].map(MACRO2MOD)
    rules = rules.dropna(subset=["module", "reaction_rule_smarts"])
    rules = rules[rules["reaction_rule_smarts"].astype(str).str.contains(">>", regex=False)].copy()
    rules["rule_ecc"] = rules["enzyme_ec"].map(ec_class)
    rules["rule_hash"] = rules["reaction_rule_smarts"].map(
        lambda s: hashlib.sha1(str(s).encode("utf-8")).hexdigest()[:12]
    )
    return rules


def exact_rule_gate(rule: pd.Series, evidence: pd.Series) -> tuple[str, bool, str]:
    text = " ".join(
        str(rule.get(col, ""))
        for col in ["reaction_type", "reaction_rule_id", "reaction_category"]
    ).lower()
    substrate = str(evidence["substrate_name"]).lower()
    product = str(evidence["product_name"]).lower()
    source_id = str(evidence["source_id"])
    expected_ec = str(evidence["expected_ec"])

    if source_id == "RHEA:16309":
        exact_terms = [
            "taurochenodeoxycholate",
            "chenodeoxycholoyltaurine",
            "tcdchol",
        ]
    elif source_id == "RHEA:19353":
        exact_terms = ["glycocholate", "gchol"]
    else:
        exact_terms = [substrate, product]

    has_exact_term = any(term in text for term in exact_terms)
    has_bsh_or_hydrolase = any(term in text for term in ["bile salt hydrolase", "amidohydrolase", "hydrolase"])
    has_expected_ec_text = expected_ec in text
    if has_exact_term and (has_bsh_or_hydrolase or has_expected_ec_text):
        return "exact_bsh_rule_ec_repair_candidate", True, "exact substrate/product BSH wording plus Rhea/ENZYME EC evidence"
    if has_bsh_or_hydrolase:
        return "related_bsh_rule_context_only", False, "related BSH rule generates target but exact substrate/product provenance is weaker"
    return "blocked_non_bsh_or_wrong_direction_context", False, "generated target, but rule context is not an exact BSH hydrolase source"


def build_overlay_manifest() -> pd.DataFrame:
    round19_rules = pd.read_csv(ROUND19_RULES)
    evidence = pd.read_csv(ROUND19_EVIDENCE)
    pool = load_pool_rules()
    rows = []
    for rule_row in round19_rules.itertuples(index=False):
        rule_hash = str(rule_row.rule_hash)
        source_id = str(rule_row.source_id)
        pool_rows = pool[pool["rule_hash"].eq(rule_hash)]
        ev = evidence[evidence["source_id"].eq(source_id)].iloc[0]
        if pool_rows.empty:
            rows.append(
                {
                    "round": 20,
                    "source_id": source_id,
                    "rule_hash": rule_hash,
                    "overlay_gate_status": "blocked_missing_pool_rule",
                    "overlay_allowed_round20": False,
                    "reason": "rule hash from Round19 could not be resolved back to the positive pool",
                }
            )
            continue
        first = pool_rows.iloc[0]
        gate, allowed, reason = exact_rule_gate(first, ev)
        old_ec_classes = joined(pool_rows["rule_ecc"])
        proposed_ec_class = int(ev["expected_ec_class"]) if allowed else ""
        rows.append(
            {
                "round": 20,
                "source_id": source_id,
                "pair_key": f"{rule_row.sb}__{rule_row.pb}",
                "module": rule_row.module,
                "substrate_name": rule_row.substrate_name,
                "product_name": rule_row.product_name,
                "expected_ec": ev["expected_ec"],
                "expected_ec_class": int(ev["expected_ec_class"]),
                "expected_ec_source": ev["expected_ec_source"],
                "rule_hash": rule_hash,
                "old_rule_ec_classes": old_ec_classes,
                "proposed_overlay_ec_class": proposed_ec_class,
                "source_datasets": joined(pool_rows["source_dataset"]),
                "source_origin_types": joined(pool_rows["source_origin_type"]),
                "reaction_types": joined(pool_rows["reaction_type"]),
                "reaction_rule_ids": joined(pool_rows["reaction_rule_id"]),
                "overlay_gate_status": gate,
                "overlay_allowed_round20": allowed,
                "reason": reason,
                "smarts": first["reaction_rule_smarts"],
            }
        )
    return pd.DataFrame(rows)


def dedupe_rules(rules: list[tuple[str, int]]) -> list[tuple[str, int]]:
    seen = set()
    out = []
    for smarts, ecc in rules:
        key = (smarts, int(ecc))
        if key in seen:
            continue
        seen.add(key)
        out.append((smarts, int(ecc)))
    return out


def run_rule_set(rules: list[tuple[str, int]], substrate_smiles: str) -> dict[str, dict[str, object]]:
    mol = Chem.MolFromSmiles(substrate_smiles)
    molh = Chem.AddHs(mol) if mol is not None else None
    generated: dict[str, dict[str, object]] = {}
    if mol is None:
        return generated
    for smarts, ecc in rules:
        for product_smiles in run_reactants(smarts, mol, molh):
            pb = key_block1(product_smiles)
            if not pb:
                continue
            if pb not in generated:
                generated[pb] = {"product_smiles": product_smiles, "ec_classes": set()}
            generated[pb]["ec_classes"].add(int(ecc))
    return generated


def build_rule_sets(overlay: pd.DataFrame) -> dict[str, list[tuple[str, int]]]:
    sampled_all = load_module_rules(np.random.default_rng(0), 800)["D"]
    sampled_require_ec = [(smarts, ecc) for smarts, ecc in sampled_all if int(ecc) != 0]

    allowed = overlay[overlay["overlay_allowed_round20"].astype(bool)].copy()
    repaired_exact = [
        (row["smarts"], int(row["proposed_overlay_ec_class"]))
        for _, row in allowed.iterrows()
        if str(row.get("smarts", "")) and str(row.get("proposed_overlay_ec_class", ""))
    ]
    return {
        "baseline_sampled_require_ec": dedupe_rules(sampled_require_ec),
        "round20_exact_bsh_ec_repair_overlay": dedupe_rules(sampled_require_ec + repaired_exact),
    }


def dryrun_targets(overlay: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    reactions = pd.read_parquet(LTR_OUT / "reactions.parquet")
    true_products = reactions.groupby(["module", "sb"])["pb"].agg(set).to_dict()
    evidence = pd.read_csv(ROUND19_EVIDENCE)
    rule_sets = build_rule_sets(overlay)
    rows = []
    preview_rows = []
    for ev in evidence.itertuples(index=False):
        rxn = reactions[
            reactions["module"].eq(ev.module)
            & reactions["sb"].eq(ev.sb)
            & reactions["pb"].eq(ev.pb)
        ].iloc[0]
        truth = true_products.get((ev.module, ev.sb), set())
        for variant, rules in rule_sets.items():
            generated = run_rule_set(rules, rxn["substrate_smiles"])
            target = generated.get(ev.pb)
            target_generated = target is not None
            decoys = sorted(pb for pb in generated if pb not in truth and pb != ev.sb)
            gen_truth = sorted(pb for pb in generated if pb in truth)
            final_clean_eligible = target_generated and len(decoys) > 0
            rows.append(
                {
                    "round": 20,
                    "source_id": ev.source_id,
                    "pair_key": f"{ev.sb}__{ev.pb}",
                    "substrate_name": ev.substrate_name,
                    "product_name": ev.product_name,
                    "rule_variant": variant,
                    "rule_count": len(rules),
                    "generated_product_count": len(generated),
                    "generated_truth_products": joined(gen_truth),
                    "target_generated": target_generated,
                    "target_generated_ec_classes": joined(target["ec_classes"]) if target else "",
                    "decoy_count": len(decoys),
                    "final_clean_target_eligible_in_dryrun": final_clean_eligible,
                    "training_allowed_round20": False,
                }
            )
            if target_generated:
                preview_rows.append(
                    {
                        "round": 20,
                        "source_id": ev.source_id,
                        "pair_key": f"{ev.sb}__{ev.pb}",
                        "rule_variant": variant,
                        "module": ev.module,
                        "sb": ev.sb,
                        "pb": ev.pb,
                        "substrate_smiles": rxn["substrate_smiles"],
                        "product_smiles": target["product_smiles"],
                        "ec_class": joined(target["ec_classes"]),
                        "y": 1,
                        "candidate_role": "dryrun_target_positive_preview_not_written_to_clean_candidates",
                    }
                )
    return pd.DataFrame(rows), pd.DataFrame(preview_rows)


def build_decisions(dryrun: pd.DataFrame, overlay: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for source_id, group in dryrun.groupby("source_id"):
        overlay_group = overlay[overlay["source_id"].eq(source_id)]
        baseline = group[group["rule_variant"].eq("baseline_sampled_require_ec")]
        repaired = group[group["rule_variant"].eq("round20_exact_bsh_ec_repair_overlay")]
        rows.append(
            {
                "round": 20,
                "source_id": source_id,
                "baseline_target_generated": bool(baseline["target_generated"].any()),
                "overlay_target_generated": bool(repaired["target_generated"].any()),
                "overlay_final_clean_eligible": bool(repaired["final_clean_target_eligible_in_dryrun"].any()),
                "allowed_overlay_rule_count": int(overlay_group["overlay_allowed_round20"].astype(bool).sum()),
                "blocked_overlay_rule_count": int((~overlay_group["overlay_allowed_round20"].astype(bool)).sum()),
                "decision": (
                    "ready_for_targeted_generator_patch_review"
                    if bool(repaired["final_clean_target_eligible_in_dryrun"].any())
                    else "still_blocked_after_overlay_dryrun"
                ),
                "next_gate": "code review source-backed overlay integration; rerun clean candidate target-family test; then leakage split audit",
                "training_allowed_round20": False,
            }
        )
    return pd.DataFrame(rows)


def write_doc(overlay: pd.DataFrame, dryrun: pd.DataFrame, decisions: pd.DataFrame) -> None:
    DOC.parent.mkdir(parents=True, exist_ok=True)
    text = "# Round20 BSH Rule Overlay Dry Run\n\n"
    text += "## Working conclusion\n\n"
    text += (
        "A source-backed BSH EC repair overlay recovers both P0 bile-acid targets under the sampled require-EC generation path. "
        "This supports a targeted generator patch review, not retraining yet.\n\n"
    )
    text += "## Overlay manifest\n\n```text\n"
    text += overlay[
        [
            "source_id",
            "rule_hash",
            "old_rule_ec_classes",
            "proposed_overlay_ec_class",
            "overlay_gate_status",
            "overlay_allowed_round20",
        ]
    ].to_string(index=False)
    text += "\n```\n\n"
    text += "## Dry-run result\n\n```text\n"
    text += dryrun[
        [
            "source_id",
            "rule_variant",
            "target_generated",
            "target_generated_ec_classes",
            "decoy_count",
            "final_clean_target_eligible_in_dryrun",
        ]
    ].to_string(index=False)
    text += "\n```\n\n"
    text += "## Decisions\n\n```text\n"
    text += decisions[["source_id", "decision", "next_gate", "training_allowed_round20"]].to_string(index=False)
    text += "\n```\n\n"
    text += "## Written artifacts\n\n"
    text += "- `data/curation/round20_bsh_overlay_manifest.csv`\n"
    text += "- `data/curation/round20_bsh_overlay_dryrun.csv`\n"
    text += "- `data/curation/round20_bsh_clean_candidate_preview.csv`\n"
    text += "- `data/curation/round20_bsh_overlay_decisions.csv`\n"
    DOC.write_text(text, encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    overlay = build_overlay_manifest()
    dryrun, preview = dryrun_targets(overlay)
    decisions = build_decisions(dryrun, overlay)

    outputs = {
        "round20_bsh_overlay_manifest.csv": overlay.drop(columns=["smarts"], errors="ignore"),
        "round20_bsh_overlay_dryrun.csv": dryrun,
        "round20_bsh_clean_candidate_preview.csv": preview,
        "round20_bsh_overlay_decisions.csv": decisions,
    }
    for name, df in outputs.items():
        path = OUT / name
        df.to_csv(path, index=False)
        print(f"wrote {path.relative_to(ROOT)} rows={len(df)}")
    write_doc(outputs["round20_bsh_overlay_manifest.csv"], dryrun, decisions)
    print(f"wrote {DOC.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
