from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd
from rdkit import Chem, RDLogger


REPO = Path(__file__).resolve().parents[1]
ROUND4_TRACE = REPO / "data/curation/targeted_rule_generation_trace_round4.csv"
RULE_MASTER = REPO.parents[2] / "Dataset/reaction_rule_master.csv"
OUT = REPO / "data/curation/master_rule_recovery_candidates_round5.csv"

sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(REPO / "modular"))
import kio  # noqa: E402
from generate_candidates import run_reactants  # noqa: E402


RDLogger.DisableLog("rdApp.*")


def _read_targets(round4_trace: Path) -> pd.DataFrame:
    trace = pd.read_csv(round4_trace)
    target = trace[trace["round4_rule_issue"].eq("master_rule_available_but_not_in_positive_pool_rules")].copy()
    if target.empty:
        return target
    grouped = (
        target.groupby(["module", "sb", "pb", "substrate_name", "product_name"], dropna=False)
        .agg(
            round4_row_ids=("row_id", lambda x: ";".join(sorted(set(map(str, x))))),
            gap_families=("gap_family", lambda x: ";".join(sorted(set(map(str, x))))),
            route_loss_stages=("route_loss_stage", lambda x: ";".join(sorted(set(map(str, x))))),
            fullrule_pair_status=("fullrule_pair_status", "first"),
        )
        .reset_index()
    )
    reactions = pd.read_parquet(REPO / "outputs/modular/ltr/reactions.parquet")
    smiles_map = (
        reactions.drop_duplicates(["module", "sb"])
        .set_index(["module", "sb"])["substrate_smiles"]
        .to_dict()
    )
    grouped["substrate_smiles"] = [smiles_map.get((r.module, r.sb), "") for r in grouped.itertuples(index=False)]
    return grouped


def _load_master(rule_master: Path, max_rules: int | None = None) -> pd.DataFrame:
    cols = [
        "rule_id",
        "rule_smarts",
        "ec_number",
        "diameter",
        "source",
        "score",
        "applicable_structure_type",
        "substrate_name",
        "product_name",
        "legacy_id",
        "ec_class",
    ]
    rules = pd.read_csv(rule_master, usecols=cols, low_memory=False).dropna(subset=["rule_smarts"])
    rules = rules[rules["rule_smarts"].astype(str).str.contains(">>", regex=False)].copy()
    rules = rules.drop_duplicates("rule_smarts").reset_index(drop=True)
    if max_rules is not None:
        rules = rules.head(max_rules).copy()
    return rules


def _lhs_pattern(smarts: str):
    lhs = smarts.split(">>", 1)[0].strip().lstrip("(").rstrip(")")
    try:
        return Chem.MolFromSmarts(lhs)
    except Exception:
        return None


def _find_hits_for_target(target: pd.Series, rules: pd.DataFrame, max_hits: int) -> list[dict[str, object]]:
    substrate_smiles = str(target["substrate_smiles"])
    module = str(target["module"])
    sb = str(target["sb"])
    pb = str(target["pb"])
    cmol = Chem.MolFromSmiles(substrate_smiles)
    if cmol is None:
        return []
    cmolh = Chem.AddHs(cmol)
    rows: list[dict[str, object]] = []
    lhs_checked = 0
    lhs_matched = 0
    rule_fired = 0
    for idx, rule in rules.iterrows():
        lhs_checked += 1
        patt = rule["_lhs_mol"]
        if patt is None:
            continue
        try:
            if not (cmol.HasSubstructMatch(patt) or cmolh.HasSubstructMatch(patt)):
                continue
        except Exception:
            continue
        lhs_matched += 1
        products = set()
        for psmi in run_reactants(str(rule["rule_smarts"]), cmol, cmolh):
            product_b1 = kio.inchikey_block1(kio.smiles_to_inchikey(psmi))
            if product_b1:
                products.add(product_b1)
        if products:
            rule_fired += 1
        if pb not in products:
            continue
        rows.append(
            {
                "module": module,
                "sb": sb,
                "pb": pb,
                "substrate_name": target["substrate_name"],
                "product_name": target["product_name"],
                "substrate_smiles": substrate_smiles,
                "gap_families": target["gap_families"],
                "round4_row_ids": target["round4_row_ids"],
                "route_loss_stages": target["route_loss_stages"],
                "rule_id": rule["rule_id"],
                "rule_smarts": rule["rule_smarts"],
                "rule_source": rule["source"],
                "rule_legacy_id": rule["legacy_id"],
                "rule_ec_number": rule["ec_number"],
                "rule_ec_class": rule["ec_class"],
                "rule_diameter": rule["diameter"],
                "rule_score": rule["score"],
                "rule_applicable_structure_type": rule["applicable_structure_type"],
                "master_rule_substrate_name": rule["substrate_name"],
                "master_rule_product_name": rule["product_name"],
                "generated_target_product": True,
                "lhs_checked_count": lhs_checked,
                "lhs_matched_count_at_hit": lhs_matched,
                "rule_fired_count_at_hit": rule_fired,
                "candidate_training_role": "rule_recovery_candidate_not_positive_label",
                "positive_label_status": "requires_exact_literature_or_database_pair_verification",
                "round5_recommended_action": "import_rule_with_provenance_after_license_and_source_reaction_check",
            }
        )
        if len(rows) >= max_hits:
            break
    if not rows:
        rows.append(
            {
                "module": module,
                "sb": sb,
                "pb": pb,
                "substrate_name": target["substrate_name"],
                "product_name": target["product_name"],
                "substrate_smiles": substrate_smiles,
                "gap_families": target["gap_families"],
                "round4_row_ids": target["round4_row_ids"],
                "route_loss_stages": target["route_loss_stages"],
                "rule_id": "",
                "rule_smarts": "",
                "rule_source": "",
                "rule_legacy_id": "",
                "rule_ec_number": "",
                "rule_ec_class": "",
                "rule_diameter": "",
                "rule_score": "",
                "rule_applicable_structure_type": "",
                "master_rule_substrate_name": "",
                "master_rule_product_name": "",
                "generated_target_product": False,
                "lhs_checked_count": lhs_checked,
                "lhs_matched_count_at_hit": lhs_matched,
                "rule_fired_count_at_hit": rule_fired,
                "candidate_training_role": "needs_manual_rule_recovery",
                "positive_label_status": "requires_exact_literature_or_database_pair_verification",
                "round5_recommended_action": "rerun_with_full_master_or_curate_family_specific_rule",
            }
        )
    return rows


def build_recovery(args: argparse.Namespace) -> pd.DataFrame:
    targets = _read_targets(args.round4_trace)
    if targets.empty:
        return pd.DataFrame()
    rules = _load_master(args.rule_master, args.max_rules)
    rules["_lhs_mol"] = rules["rule_smarts"].map(_lhs_pattern)
    rows: list[dict[str, object]] = []
    for i, target in targets.iterrows():
        rows.extend(_find_hits_for_target(target, rules, args.max_hits_per_pair))
        print(f"[{i + 1}/{len(targets)}] {target['substrate_name']} -> {target['product_name']} hits={len(rows)}")
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description="Extract master-rule recovery candidates for Round 4 route losses.")
    parser.add_argument("--round4-trace", type=Path, default=ROUND4_TRACE)
    parser.add_argument("--rule-master", type=Path, default=RULE_MASTER)
    parser.add_argument("--out", type=Path, default=OUT)
    parser.add_argument("--max-rules", type=int, default=None)
    parser.add_argument("--max-hits-per-pair", type=int, default=5)
    args = parser.parse_args()

    out = build_recovery(args)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out, index=False)
    print(
        {
            "rows": len(out),
            "unique_pairs": out[["module", "sb", "pb"]].drop_duplicates().shape[0] if not out.empty else 0,
            "hit_rows": int(out.get("generated_target_product", pd.Series(dtype=bool)).fillna(False).sum()) if not out.empty else 0,
            "sources": out.get("rule_source", pd.Series(dtype=str)).value_counts(dropna=False).to_dict() if not out.empty else {},
        }
    )
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
