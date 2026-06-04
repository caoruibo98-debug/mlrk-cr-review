from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger


REPO = Path(__file__).resolve().parents[1]
ROUND3 = REPO / "data/curation/sample_route_trace_round3.csv"
OUT = REPO / "data/curation/targeted_rule_generation_trace_round4.csv"

sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(REPO / "modular"))
import kio  # noqa: E402
from generate_candidates import run_reactants  # noqa: E402
from ltr_candidates import load_module_rules  # noqa: E402


RDLogger.DisableLog("rdApp.*")

ACTIONABLE_STAGES = {
    "selected_substrate_missing_from_clean_output",
    "generated_in_audit_but_missing_from_clean",
    "target_product_not_generated",
    "generation_miss_drops_substrate",
    "no_true_product_generated_for_substrate",
    "not_covered_by_fullrule_audit",
}


def _as_bool(value: object) -> bool:
    if isinstance(value, bool):
        return value
    if pd.isna(value):
        return False
    if isinstance(value, (int, float)):
        return value > 0
    return str(value).strip().lower() in {"true", "1", "yes", "y"}


def _generated_products(smiles: str, module: str, rules: dict[str, list[tuple[str, int]]]) -> set[str]:
    cmol = Chem.MolFromSmiles(smiles)
    if cmol is None:
        return set()
    cmolh = Chem.AddHs(cmol)
    generated: set[str] = set()
    for smarts, _ecc in rules.get(module, []):
        for psmi in run_reactants(smarts, cmol, cmolh):
            pb = kio.inchikey_block1(kio.smiles_to_inchikey(psmi))
            if pb:
                generated.add(pb)
    return generated


def _classify(row: pd.Series) -> str:
    if _as_bool(row.get("exact_holdout_pair")):
        return "holdout_do_not_train"
    if _as_bool(row.get("target_in_clean_positive")):
        return "already_in_clean_candidates"
    if not _as_bool(row.get("in_reactions_parquet")):
        return "upstream_positive_pool_filter"
    if _as_bool(row.get("clean_default_hit")):
        return "post_generation_filter_or_duplicate_drop"
    if _as_bool(row.get("all_pool_rule_hit")):
        return "rule_sampling_gap_in_clean_builder"
    if _as_bool(row.get("fullrule_pair_hit")):
        return "master_rule_available_but_not_in_positive_pool_rules"
    return "rule_family_absent_or_not_validated_locally"


def _recommend(issue: str) -> str:
    return {
        "holdout_do_not_train": "Keep as evaluation only; do not train exact pair.",
        "already_in_clean_candidates": "Use as positive control, not a blocker.",
        "upstream_positive_pool_filter": "Audit single-step, structure, module, and canonicalization filters before rule work.",
        "post_generation_filter_or_duplicate_drop": "Inspect clean_candidates_full output filtering, duplicate handling, and y labeling.",
        "rule_sampling_gap_in_clean_builder": "Replace random rule sampling with deterministic evidence-tiered rule inclusion.",
        "master_rule_available_but_not_in_positive_pool_rules": "Import evidence-backed RetroRules/Rhea/MetaNetX rule provenance into the module rule pool.",
        "rule_family_absent_or_not_validated_locally": "Curate exact literature/database reaction and add family-specific rule only after provenance verification.",
    }[issue]


def build_trace(round3_path: Path) -> pd.DataFrame:
    trace = pd.read_csv(round3_path)
    target = trace[trace["route_loss_stage"].isin(ACTIONABLE_STAGES)].copy()
    reactions = pd.read_parquet(REPO / "outputs/modular/ltr/reactions.parquet")
    clean_rules = load_module_rules(np.random.default_rng(0), 800)
    broad_rules = load_module_rules(np.random.default_rng(0), 6000)
    all_pool_rules = load_module_rules(np.random.default_rng(0), 1_000_000)

    generation_cache: dict[tuple[str, str, str], dict[str, set[str]]] = {}
    rows: list[dict[str, object]] = []
    for source in target.itertuples(index=False):
        module = str(source.module)
        sb = str(source.sb)
        pb = str(source.pb)
        # Prefer the canonical substrate SMILES from reactions.parquet when available.
        match = reactions[reactions["module"].eq(module) & reactions["sb"].eq(sb)]
        substrate_smiles = str(match.iloc[0]["substrate_smiles"]) if not match.empty else ""
        key = (module, sb, substrate_smiles)
        if key not in generation_cache:
            generation_cache[key] = {
                "clean": _generated_products(substrate_smiles, module, clean_rules),
                "broad": _generated_products(substrate_smiles, module, broad_rules),
                "all_pool": _generated_products(substrate_smiles, module, all_pool_rules),
            }
        generated = generation_cache[key]
        row = source._asdict()
        row.update(
            {
                "substrate_smiles_used": substrate_smiles,
                "clean_default_rule_count": len(clean_rules.get(module, [])),
                "broad_rule_count": len(broad_rules.get(module, [])),
                "all_pool_rule_count": len(all_pool_rules.get(module, [])),
                "clean_default_generated_count": len(generated["clean"]),
                "broad_generated_count": len(generated["broad"]),
                "all_pool_generated_count": len(generated["all_pool"]),
                "clean_default_hit": pb in generated["clean"],
                "broad_rule_hit": pb in generated["broad"],
                "all_pool_rule_hit": pb in generated["all_pool"],
            }
        )
        issue = _classify(pd.Series(row))
        row["round4_rule_issue"] = issue
        row["round4_recommended_action"] = _recommend(issue)
        rows.append(row)
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description="Trace Round 3 losses against clean/broad/all module rule pools.")
    parser.add_argument("--round3", type=Path, default=ROUND3)
    parser.add_argument("--out", type=Path, default=OUT)
    args = parser.parse_args()
    out = build_trace(args.round3)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out, index=False)
    print(
        {
            "rows": len(out),
            "issue_counts": out["round4_rule_issue"].value_counts(dropna=False).to_dict(),
            "hit_counts": {
                "clean_default_hit": int(out["clean_default_hit"].sum()),
                "broad_rule_hit": int(out["broad_rule_hit"].sum()),
                "all_pool_rule_hit": int(out["all_pool_rule_hit"].sum()),
            },
        }
    )
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
