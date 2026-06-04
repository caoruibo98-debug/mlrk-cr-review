from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


REPO = Path(__file__).resolve().parents[1]
LTR_OUT = REPO / "outputs/modular/ltr"


def _block1(value: object) -> str:
    if pd.isna(value):
        return ""
    return str(value).split("-")[0].strip()


def _as_bool(value: object) -> bool:
    if isinstance(value, bool):
        return value
    if pd.isna(value):
        return False
    text = str(value).strip().lower()
    return text in {"true", "1", "yes", "y"}


def _truth_set_map(reactions: pd.DataFrame) -> dict[tuple[str, str], set[str]]:
    return reactions.groupby(["module", "sb"])["pb"].agg(set).to_dict()


def _selected_substrates(reactions: pd.DataFrame, weak_cap: int) -> set[tuple[str, str]]:
    rng = np.random.default_rng(0)
    selected: set[tuple[str, str]] = set()
    for module in ["A", "B", "C", "D"]:
        sub = reactions[reactions["module"].eq(module)].drop_duplicates("sb").copy()
        if "tier" not in sub.columns:
            sub["tier"] = np.where(sub["is_clean"], "gold", "weak")
        prio = sub[sub["tier"].isin(["gold", "silver"])]
        weak = sub[sub["tier"].eq("weak")]
        if len(weak) > weak_cap:
            weak = weak.sample(n=weak_cap, random_state=0)
        chosen = pd.concat([prio, weak]).drop_duplicates("sb")
        selected.update((module, sb) for sb in chosen["sb"].astype(str))
    return selected


def _fullrule_status(fullrule: pd.DataFrame, module: str, sb: str, pb: str) -> tuple[str, int | str, int]:
    if fullrule.empty or not {"module", "sb", "pb", "hit"}.issubset(fullrule.columns):
        return "not_available", "", 0
    same_pair = fullrule[fullrule["module"].eq(module) & fullrule["sb"].eq(sb) & fullrule["pb"].eq(pb)]
    same_sub = fullrule[fullrule["module"].eq(module) & fullrule["sb"].eq(sb)]
    same_sub_hits = int(same_sub["hit"].sum()) if not same_sub.empty else 0
    if same_pair.empty:
        return "pair_not_in_fullrule_audit", "", same_sub_hits
    hit = int(same_pair["hit"].max())
    return "pair_hit" if hit else "pair_miss", hit, same_sub_hits


def _classify(row: pd.Series) -> tuple[str, str]:
    status = row["candidate_status"]
    if status == "do_not_train_exact_holdout" or _as_bool(row["exact_holdout_pair"]):
        return (
            "holdout_leakage_guard",
            "Keep as evaluation/diagnosis only; replace with independent non-holdout cases before training.",
        )
    if status == "already_in_clean_candidate_pool" or _as_bool(row["target_in_clean_positive"]):
        return (
            "already_reaches_clean_candidates",
            "Use as a coverage control; this row is not the current production blocker.",
        )
    if status == "absent_from_local_positive_pool":
        return (
            "external_curation_needed",
            "Curate exact substrate/product pairs from cited literature or databases before training.",
        )
    if not _as_bool(row["in_reactions_parquet"]):
        return (
            "positive_pool_to_reactions_filter",
            "Inspect ltr_build.py filters: single-step, structure_parseable, module mapping, canonicalization, and tier assignment.",
        )
    if not _as_bool(row["selected_by_clean_builder"]):
        return (
            "substrate_not_selected_by_clean_builder",
            "Inspect weak sampling or substrate-level deduplication before candidate generation.",
        )
    if not _as_bool(row["same_substrate_in_clean_candidates_full"]):
        if row["fullrule_pair_status"] == "pair_miss":
            return (
                "generation_miss_drops_substrate",
                "Candidate-generation rules do not generate this true product; add or import evidence-backed reaction rules before retraining.",
            )
        if int(row["fullrule_same_substrate_hits"]) == 0:
            return (
                "no_true_product_generated_for_substrate",
                "No true products are generated for the selected substrate, so clean candidate generation drops the substrate.",
            )
        return (
            "selected_substrate_missing_from_clean_output",
            "Run a targeted generation trace; selected substrate has evidence but no clean output row.",
        )
    if row["fullrule_pair_status"] == "pair_miss":
        return (
            "target_product_not_generated",
            "Other products for the same substrate may be generated, but this target product is not; add a family-specific rule.",
        )
    if row["fullrule_pair_status"] == "pair_not_in_fullrule_audit":
        return (
            "not_covered_by_fullrule_audit",
            "Add this reaction to the generation audit or run targeted rule generation before deciding training status.",
        )
    if row["fullrule_pair_status"] == "pair_hit":
        return (
            "generated_in_audit_but_missing_from_clean",
            "Inspect canonicalization, clean_candidates_full generation parameters, duplicate handling, and post-generation filtering.",
        )
    return (
        "unclassified_route_gap",
        "Needs manual route review before training promotion.",
    )


def build_trace(args: argparse.Namespace) -> pd.DataFrame:
    round2 = pd.read_csv(args.round2_manifest)
    reactions = pd.read_parquet(args.ltr_out / "reactions.parquet")
    clean = pd.read_parquet(args.ltr_out / "clean_candidates_full.parquet")
    fullrule_path = args.kernel_root / "outputs/gen_gap/fullrule_hit_miss.parquet"
    fullrule = pd.read_parquet(fullrule_path) if fullrule_path.exists() else pd.DataFrame()

    selected = _selected_substrates(reactions, args.weak_cap)
    truth = _truth_set_map(reactions)
    clean_pos = clean[clean["y"].eq(1)].copy()
    clean_pos_pairs = set(zip(clean_pos["module"], clean_pos["sb"], clean_pos["pb"]))
    clean_substrate_pairs = set(zip(clean["module"], clean["sb"]))

    rows: list[dict[str, object]] = []
    for _, source in round2.iterrows():
        module = str(source.get("macro_module", ""))
        # Round 2 stores A/B/C/D as macro_module only when inherited from reactions; fall back from reaction row.
        sb = _block1(source.get("substrate_inchikey"))
        pb = _block1(source.get("product_inchikey"))
        reaction_rows = reactions[reactions["sb"].eq(sb) & reactions["pb"].eq(pb)]
        if not reaction_rows.empty:
            module = str(reaction_rows.iloc[0]["module"])
        elif module in {"small_molecule_polyphenol", "carbohydrate_glycan", "protein_amino_acid", "lipid_fat"}:
            module = {
                "small_molecule_polyphenol": "A",
                "carbohydrate_glycan": "B",
                "protein_amino_acid": "C",
                "lipid_fat": "D",
            }[module]
        elif module not in {"A", "B", "C", "D"}:
            module = ""

        fullrule_pair_status, fullrule_pair_hit, same_sub_hits = _fullrule_status(fullrule, module, sb, pb)
        clean_same_sub = clean[clean["module"].eq(module) & clean["sb"].eq(sb)] if module else pd.DataFrame()
        row = {
            "row_id": source.get("row_id", ""),
            "gap_family": source.get("gap_family", ""),
            "candidate_status": source.get("candidate_status", ""),
            "module": module,
            "substrate_name": source.get("substrate_name", ""),
            "product_name": source.get("product_name", ""),
            "sb": sb,
            "pb": pb,
            "in_reactions_parquet": _as_bool(source.get("in_reactions_parquet", False)),
            "reaction_tiers": ";".join(sorted(reaction_rows["tier"].dropna().astype(str).unique())) if not reaction_rows.empty else "",
            "selected_by_clean_builder": (module, sb) in selected if module and sb else False,
            "target_in_clean_positive": (module, sb, pb) in clean_pos_pairs if module and sb and pb else False,
            "same_substrate_in_clean_candidates_full": (module, sb) in clean_substrate_pairs if module and sb else False,
            "clean_positive_products_for_substrate": int(clean_same_sub[clean_same_sub["y"].eq(1)]["pb"].nunique()) if not clean_same_sub.empty else 0,
            "truth_products_for_substrate": len(truth.get((module, sb), set())) if module and sb else 0,
            "fullrule_pair_status": fullrule_pair_status,
            "fullrule_pair_hit": fullrule_pair_hit,
            "fullrule_same_substrate_hits": same_sub_hits,
            "exact_holdout_pair": _as_bool(source.get("exact_holdout_pair", False)),
            "source_id": source.get("source_id", ""),
            "validation_status_round2": source.get("validation_status", ""),
        }
        loss_stage, recommended_fix = _classify(pd.Series(row))
        row["route_loss_stage"] = loss_stage
        row["recommended_fix"] = recommended_fix
        rows.append(row)
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description="Trace Round 2 sample candidates through the LTR route.")
    parser.add_argument("--round2-manifest", type=Path, default=REPO / "data/curation/reaction_gap_positive_manifest_round2.csv")
    parser.add_argument("--ltr-out", type=Path, default=LTR_OUT)
    parser.add_argument("--kernel-root", type=Path, default=REPO)
    parser.add_argument("--weak-cap", type=int, default=600)
    parser.add_argument("--out", type=Path, default=REPO / "data/curation/sample_route_trace_round3.csv")
    args = parser.parse_args()

    trace = build_trace(args)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    trace.to_csv(args.out, index=False)
    print(
        {
            "rows": len(trace),
            "families": trace["gap_family"].nunique(),
            "route_loss_stage_counts": trace["route_loss_stage"].value_counts(dropna=False).to_dict(),
        }
    )
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
