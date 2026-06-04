#!/usr/bin/env python
"""Convert Round7 Rhea checks into an auditable CSV."""

from __future__ import annotations

import json
import re
from pathlib import Path

import pandas as pd


REPO = Path(__file__).resolve().parents[1]
CURATION = REPO / "data" / "curation"
SOURCE_CHECK = CURATION / "source_reaction_check_round7.csv"
RAW = CURATION / "rhea_round7_raw.json"
OUT = CURATION / "rhea_source_reaction_summary_round7.csv"


def extract_rhea_ids(value: object) -> list[str]:
    if value is None or pd.isna(value):
        return []
    ids = []
    for match in re.finditer(r"rhea:(\d+)", str(value), flags=re.I):
        rid = match.group(1)
        if rid not in ids:
            ids.append(rid)
    return ids


def pair_exactness(equation: str, substrate: str, product: str) -> str:
    eq = equation.lower()
    sub = str(substrate).lower()
    prod = str(product).lower()
    has_sub = sub in eq
    has_prod = prod in eq
    if has_sub and has_prod:
        return "rhea_equation_names_both_pair_compounds"
    if has_sub or has_prod:
        return "rhea_equation_names_one_pair_compound"
    return "rhea_equation_not_exact_pair_by_name"


def missing_rhea_status(rid: str, returned_ids: set[str]) -> str:
    """Rhea often exposes directional companion ids around a master reaction id."""
    try:
        numeric = int(rid)
    except ValueError:
        return "rhea_xref_not_returned_live"
    for offset in (-3, -2, -1, 1, 2, 3):
        if str(numeric + offset) in returned_ids:
            return "rhea_directional_companion_not_returned_live"
    return "rhea_xref_not_returned_live"


def main() -> None:
    source = pd.read_csv(SOURCE_CHECK)
    raw = json.loads(RAW.read_text(encoding="utf-8"))
    results = raw.get("results") or raw.get("records") or []
    by_id = {str(r.get("id")): r for r in results if isinstance(r, dict)}
    returned_ids = set(by_id)

    rows: list[dict[str, object]] = []
    for item in source.itertuples(index=False):
        rhea_ids = extract_rhea_ids(item.rhea_xrefs)
        if not rhea_ids:
            rows.append(
                {
                    "round": 7,
                    "substrate_name": item.substrate_name,
                    "product_name": item.product_name,
                    "rule_id": item.rule_id,
                    "rhea_id": "",
                    "rhea_found_live_round7": False,
                    "rhea_status": "",
                    "rhea_equation": "",
                    "rhea_balanced": "",
                    "rhea_pair_exactness_by_name": "no_rhea_xref",
                    "round7_training_implication": "no_direct_rhea_support_for_pair",
                }
            )
            continue
        for rid in rhea_ids:
            rec = by_id.get(rid, {})
            equation = str(rec.get("equation", ""))
            exactness = (
                pair_exactness(equation, item.substrate_name, item.product_name)
                if equation
                else missing_rhea_status(rid, returned_ids)
            )
            implication = "do_not_train_positive_from_rhea_alone"
            if exactness == "rhea_equation_names_both_pair_compounds":
                implication = "candidate_exact_pair_check_still_needs_direction_stereo_context"
            rows.append(
                {
                    "round": 7,
                    "substrate_name": item.substrate_name,
                    "product_name": item.product_name,
                    "rule_id": item.rule_id,
                    "rhea_id": rid,
                    "rhea_found_live_round7": bool(rec),
                    "rhea_status": rec.get("status", ""),
                    "rhea_equation": equation,
                    "rhea_balanced": rec.get("balanced", ""),
                    "rhea_pair_exactness_by_name": exactness,
                    "round7_training_implication": implication,
                }
            )

    out = pd.DataFrame(rows)
    out.to_csv(OUT, index=False)
    print(f"wrote {OUT} rows={len(out)}")
    print(out["rhea_pair_exactness_by_name"].value_counts(dropna=False).to_string())


if __name__ == "__main__":
    main()
