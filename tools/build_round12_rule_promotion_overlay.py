#!/usr/bin/env python
"""Build a gated generator-rule promotion overlay manifest.

Round12 does not alter training labels. It decides which recovered rules are
safe enough for a deployment-generator dry run and which remain blocked by
source provenance, exact-pair, direction, or license gates.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd


REPO = Path(__file__).resolve().parents[1]
CURATION = REPO / "data" / "curation"

ROUND11_GENERATOR = CURATION / "known_gold_generator_trace_round11.csv"
ROUND11_DECISIONS = CURATION / "round11_training_gate_decisions.csv"
ROUND5_RECOVERY = CURATION / "master_rule_recovery_candidates_round5.csv"
ROUND7_SOURCE = CURATION / "source_reaction_check_round7.csv"

OUT_OVERLAY = CURATION / "rule_promotion_overlay_round12.csv"
OUT_DRYRUN = CURATION / "rule_promotion_dryrun_plan_round12.csv"
OUT_SUMMARY = CURATION / "round12_rule_promotion_summary.csv"


def as_bool(value: object) -> bool:
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"true", "1", "yes"}


def choose_candidate_rules(recovery: pd.DataFrame) -> pd.DataFrame:
    if recovery.empty:
        return recovery
    df = recovery.copy()
    df["generated_target_product_bool"] = df.get("generated_target_product", False).map(as_bool)
    df["source_rank"] = df["rule_source"].map({"retrorules": 0, "rhea": 1, "metanetx": 1, "microberx": 5}).fillna(9)
    df["ec_rank"] = df["rule_ec_number"].notna().map({True: 0, False: 1})
    df["score_sort"] = pd.to_numeric(df.get("rule_score", 0), errors="coerce").fillna(0)
    df = df[df["generated_target_product_bool"]].copy()
    df = df.sort_values(["pair_key", "source_rank", "ec_rank", "score_sort", "rule_id"], ascending=[True, True, True, False, True])
    return df.groupby("pair_key", as_index=False).head(3)


def status_for(row: pd.Series) -> tuple[str, bool, bool, str]:
    exact_count = int(row.get("exact_source_candidate_count", 0) or 0)
    source_status = str(row.get("source_reaction_check_status", ""))
    rule_source = str(row.get("rule_source", ""))
    has_ec = bool(str(row.get("rule_ec_number", "")).strip() and str(row.get("rule_ec_number", "")).strip().lower() != "nan")

    if exact_count <= 0:
        return (
            "blocked_no_exact_pair_source",
            False,
            False,
            "Exact substrate-product evidence is not strong enough for a generator overlay tied to this pair.",
        )
    if "verified_not_imported" in source_status and rule_source == "retrorules" and has_ec:
        return (
            "dryrun_candidate_source_traceable_rule",
            True,
            False,
            "RetroRules source reaction is verified and exact-pair evidence exists; dry-run overlay is allowed, but deployment promotion still needs license/direction/stereochemistry review.",
        )
    if rule_source == "microberx":
        return (
            "blocked_microberx_source_recovery_required",
            False,
            False,
            "Rule can generate target, but source reaction provenance is MicrobeRX-only or unrecovered.",
        )
    return (
        "blocked_manual_source_direction_license_review",
        False,
        False,
        "Rule needs manual source, direction, stereochemistry, and license review before any overlay dry run.",
    )


def build_overlay() -> pd.DataFrame:
    generator = pd.read_csv(ROUND11_GENERATOR)
    decisions = pd.read_csv(ROUND11_DECISIONS)
    recovery = pd.read_csv(ROUND5_RECOVERY)
    source = pd.read_csv(ROUND7_SOURCE)
    recovery["pair_key"] = recovery["sb"].astype(str) + "__" + recovery["pb"].astype(str)
    source["pair_key"] = source["sb"].astype(str) + "__" + source["pb"].astype(str)
    candidate_rules = choose_candidate_rules(recovery)

    rows = []
    for gen in generator.itertuples(index=False):
        key = gen.pair_key
        decision = decisions[decisions["pair_key"].eq(key)].head(1)
        source_rows = source[source["pair_key"].eq(key)]
        rules = candidate_rules[candidate_rules["pair_key"].eq(key)]
        if rules.empty:
            rows.append(
                {
                    "round": 12,
                    "pair_key": key,
                    "module": gen.module,
                    "substrate_name": gen.substrate_name,
                    "product_name": gen.product_name,
                    "rule_id": "",
                    "rule_source": "",
                    "rule_ec_number": "",
                    "rule_legacy_id": "",
                    "source_reaction_check_status": "",
                    "exact_source_candidate_count": int(decision["exact_source_candidate_count"].iloc[0]) if not decision.empty else 0,
                    "rule_generated_target": False,
                    "overlay_status_round12": "blocked_no_recovered_rule_row",
                    "overlay_dryrun_allowed_round12": False,
                    "deployment_promotion_allowed_round12": False,
                    "blocker_or_note": "No recovered rule row was found for this known gold miss.",
                    "candidate_generator_role": "none",
                }
            )
            continue
        for _, rule in rules.iterrows():
            src = source_rows[source_rows["rule_id"].eq(rule["rule_id"])].head(1)
            if src.empty and "rule_legacy_id" in source_rows.columns:
                src = source_rows[source_rows["rule_legacy_id"].astype(str).eq(str(rule.get("rule_legacy_id", "")))].head(1)
            merged = rule.to_dict()
            merged.update(
                {
                    "source_reaction_check_status": src["source_reaction_check_status"].iloc[0] if not src.empty else "",
                    "exact_source_candidate_count": int(decision["exact_source_candidate_count"].iloc[0]) if not decision.empty else 0,
                }
            )
            status, dryrun, deployment, note = status_for(pd.Series(merged))
            rows.append(
                {
                    "round": 12,
                    "pair_key": key,
                    "module": gen.module,
                    "substrate_name": gen.substrate_name,
                    "product_name": gen.product_name,
                    "rule_id": rule.get("rule_id", ""),
                    "rule_source": rule.get("rule_source", ""),
                    "rule_ec_number": rule.get("rule_ec_number", ""),
                    "rule_ec_class": rule.get("rule_ec_class", ""),
                    "rule_legacy_id": rule.get("rule_legacy_id", ""),
                    "rule_diameter": rule.get("rule_diameter", ""),
                    "rule_score": rule.get("rule_score", ""),
                    "source_reaction_check_status": merged["source_reaction_check_status"],
                    "exact_source_candidate_count": merged["exact_source_candidate_count"],
                    "rule_generated_target": as_bool(rule.get("generated_target_product", False)),
                    "overlay_status_round12": status,
                    "overlay_dryrun_allowed_round12": dryrun,
                    "deployment_promotion_allowed_round12": deployment,
                    "blocker_or_note": note,
                    "candidate_generator_role": "source_traceable_generator_overlay" if dryrun else "blocked_generator_overlay_candidate",
                }
            )
    return pd.DataFrame(rows)


def build_dryrun_plan(overlay: pd.DataFrame) -> pd.DataFrame:
    rows = []
    allowed = overlay[overlay["overlay_dryrun_allowed_round12"].map(as_bool)].copy()
    if allowed.empty:
        rows.append(
            {
                "round": 12,
                "plan_status": "no_dryrun_rules_allowed",
                "pair_key": "",
                "rule_id": "",
                "dryrun_step": "continue source recovery before generator modification",
                "success_criterion": "at least one source-traceable rule passes exact/source/license gates",
            }
        )
    else:
        for item in allowed.itertuples(index=False):
            rows.append(
                {
                    "round": 12,
                    "plan_status": "dryrun_rule_overlay_candidate",
                    "pair_key": item.pair_key,
                    "rule_id": item.rule_id,
                    "dryrun_step": "add this rule to a temporary generator overlay and rerun clean candidate generation for the target substrate only",
                    "success_criterion": "target product appears in clean_candidates_full without duplicate positive import",
                }
            )
    return pd.DataFrame(rows)


def build_summary(overlay: pd.DataFrame) -> pd.DataFrame:
    return (
        overlay.groupby("overlay_status_round12", dropna=False)
        .agg(
            n_rules=("rule_id", "size"),
            n_pairs=("pair_key", "nunique"),
            dryrun_allowed=("overlay_dryrun_allowed_round12", lambda x: int(pd.Series(x).map(as_bool).sum())),
            deployment_allowed=("deployment_promotion_allowed_round12", lambda x: int(pd.Series(x).map(as_bool).sum())),
        )
        .reset_index()
        .assign(round=12)
    )


def main() -> None:
    overlay = build_overlay()
    dryrun = build_dryrun_plan(overlay)
    summary = build_summary(overlay)
    for path, table in [
        (OUT_OVERLAY, overlay),
        (OUT_DRYRUN, dryrun),
        (OUT_SUMMARY, summary),
    ]:
        table.to_csv(path, index=False)
        print(f"wrote {path} rows={len(table)}")
    print("round12_overlay_status=", overlay["overlay_status_round12"].value_counts().to_dict())
    print("round12_dryrun_allowed=", int(overlay["overlay_dryrun_allowed_round12"].map(as_bool).sum()))


if __name__ == "__main__":
    main()
