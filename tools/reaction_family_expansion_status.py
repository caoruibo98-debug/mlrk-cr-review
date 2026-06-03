from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import Counter
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from mlrk_prod.io_utils import atomic_write_json, atomic_write_text  # noqa: E402

SEED_MANIFEST = ROOT / "data" / "reaction_family_gap_seed_manifest.tsv"
REACTION_KPIS = ROOT / "outputs" / "appraisal" / "reaction_family_kpis.json"
DEFAULT_JSON_OUT = ROOT / "outputs" / "appraisal" / "reaction_family_expansion_status.json"
DEFAULT_MD_OUT = ROOT / "docs" / "reviews" / "reaction_family_expansion_status.md"

REQUIRED_SEED_FIELDS = [
    "case_id",
    "split_role",
    "training_allowed",
    "holdout_policy",
    "reaction_family",
    "substrate_name_en",
    "substrate_name_zh",
    "substrate_smiles",
    "substrate_inchikey",
    "product_name_en",
    "product_name_zh",
    "product_smiles",
    "product_inchikey",
    "transformation_class",
    "evidence_type",
    "reference",
    "source_url",
    "evidence_strength",
    "current_gap",
    "required_model_change",
    "representation_ready",
]

FAMILY_TRAINING_TARGETS = {
    "phenolic_ester_hydrolysis": 20,
    "ellagitannin_hydrolysis": 20,
    "ellagitannin_urolithin_multistep": 30,
    "isoflavone_reduction": 25,
    "isoflavone_reductive_metabolism": 30,
    "isoflavone_c_glycoside_conversion": 20,
    "lignan_deglucosylation": 20,
    "lignan_multistep_demethylation_dehydroxylation": 30,
    "lignan_oxidation": 15,
    "flavanol_ring_fission": 30,
    "flavanone_ring_fission": 25,
    "flavonol_ring_fission": 25,
}


def read_seed(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        missing = [field for field in REQUIRED_SEED_FIELDS if field not in (reader.fieldnames or [])]
        if missing:
            raise ValueError(f"Seed manifest is missing required fields: {missing}")
        return list(reader)


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def as_bool(value: str) -> bool:
    return str(value).strip().lower() in {"true", "1", "yes"}


def family_kpi_map(report: dict[str, Any]) -> dict[str, dict[str, Any]]:
    rows = report.get("combined_family_kpis", [])
    return {str(row.get("reaction_family")): row for row in rows}


def build_report(seed_rows: list[dict[str, str]], kpis: dict[str, Any]) -> dict[str, Any]:
    kpi_by_family = family_kpi_map(kpis)
    seed_by_family: dict[str, list[dict[str, str]]] = {}
    for row in seed_rows:
        seed_by_family.setdefault(row["reaction_family"], []).append(row)

    blocked_families = sorted(
        family
        for family, row in kpi_by_family.items()
        if row.get("production_gap") == "candidate_generation_blocked"
    )
    training_allowed = [row for row in seed_rows if as_bool(row["training_allowed"])]
    holdout = [row for row in seed_rows if not as_bool(row["training_allowed"])]
    seed_family_counts = Counter(row["reaction_family"] for row in seed_rows)
    transformation_counts = Counter(row["transformation_class"] for row in seed_rows)

    family_rows = []
    for family in blocked_families:
        kpi = kpi_by_family.get(family, {})
        seeds = seed_by_family.get(family, [])
        train_seeds = [row for row in seeds if as_bool(row["training_allowed"])]
        holdout_seeds = [row for row in seeds if not as_bool(row["training_allowed"])]
        target = FAMILY_TRAINING_TARGETS.get(family, 15)
        family_rows.append(
            {
                "reaction_family": family,
                "current_case_count": kpi.get("case_count", 0),
                "current_candidate_pool_recall": kpi.get("candidate_pool_recall", 0),
                "current_top5_hit_rate": kpi.get("top5_hit_rate", 0),
                "seed_rows": len(seeds),
                "holdout_seed_rows": len(holdout_seeds),
                "training_allowed_seed_rows": len(train_seeds),
                "target_independent_training_pairs": target,
                "additional_training_pairs_needed": max(0, target - len(train_seeds)),
                "current_status": "seeded_holdout_only" if seeds and not train_seeds else "needs_seed",
                "next_action": kpi.get("recommended_next_action", "Add independently sourced training cases."),
                "required_model_changes": sorted({row["required_model_change"] for row in seeds}),
            }
        )

    blocked_without_seed = [row["reaction_family"] for row in family_rows if row["seed_rows"] == 0]
    ready_for_retraining = [
        row["reaction_family"]
        for row in family_rows
        if row["training_allowed_seed_rows"] >= row["target_independent_training_pairs"]
    ]

    return {
        "status": "ready",
        "branch_intent": "reaction-family-expansion",
        "source_seed_manifest": str(SEED_MANIFEST.relative_to(ROOT)).replace("\\", "/"),
        "source_reaction_kpis": str(REACTION_KPIS.relative_to(ROOT)).replace("\\", "/"),
        "policy": {
            "challenge_holdouts_are_not_training_data": True,
            "reason": "Seed rows resolve known blocked benchmark cases into machine-readable structures, but must be replaced by independent cases before training to avoid leakage.",
        },
        "seed_summary": {
            "rows": len(seed_rows),
            "families": len(seed_family_counts),
            "training_allowed_rows": len(training_allowed),
            "holdout_rows": len(holdout),
            "family_counts": dict(sorted(seed_family_counts.items())),
            "transformation_counts": dict(sorted(transformation_counts.items())),
        },
        "blocked_family_count": len(blocked_families),
        "blocked_families": blocked_families,
        "blocked_without_seed": blocked_without_seed,
        "ready_for_retraining_families": ready_for_retraining,
        "family_status": family_rows,
        "overall_next_steps": [
            "Collect independent non-holdout positive pairs for each seeded family.",
            "Attach source IDs, English names, Chinese names, SMILES, InChIKey, PMIDs/PMCIDs/DOIs, enzymes, microbes, and evidence tiers.",
            "Add reaction-family candidate-generation rules only after at least one independent training seed and one holdout case exist for that family.",
            "Benchmark ChemBERTa-2 against MoLFormer and MolT5/ChemT5-style encoders after the new cases are versioned.",
            "Retrain only after data audit confirms no challenge-holdout leakage.",
        ],
    }


def render_markdown(report: dict[str, Any]) -> str:
    rows = report["family_status"]
    table_lines = [
        "| Reaction family | Holdout seeds | Training seeds | Target independent training pairs | Additional needed | Current status |",
        "| --- | ---: | ---: | ---: | ---: | --- |",
    ]
    for row in rows:
        table_lines.append(
            "| {reaction_family} | {holdout_seed_rows} | {training_allowed_seed_rows} | {target_independent_training_pairs} | {additional_training_pairs_needed} | {current_status} |".format(
                **row
            )
        )
    next_steps = "\n".join(f"- {item}" for item in report["overall_next_steps"])
    return f"""# Reaction Family Expansion Status

## Version Intent

- Branch intent: `{report['branch_intent']}`
- Seed manifest: `{report['source_seed_manifest']}`
- KPI source: `{report['source_reaction_kpis']}`
- Challenge holdouts are training data: `false`

## Seed Summary

- Seed rows: `{report['seed_summary']['rows']}`
- Seeded families: `{report['seed_summary']['families']}`
- Training-allowed rows: `{report['seed_summary']['training_allowed_rows']}`
- Holdout rows: `{report['seed_summary']['holdout_rows']}`
- Blocked families with no seed: `{report['blocked_without_seed']}`
- Ready for retraining families: `{report['ready_for_retraining_families']}`

## Family Status

{chr(10).join(table_lines)}

## Policy

{report['policy']['reason']}

## Next Steps

{next_steps}
"""


def write_outputs(report: dict[str, Any], json_out: Path, md_out: Path) -> None:
    atomic_write_json(json_out, report)
    atomic_write_text(md_out, render_markdown(report))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--json-out", default=str(DEFAULT_JSON_OUT.relative_to(ROOT)).replace("\\", "/"))
    parser.add_argument("--md-out", default=str(DEFAULT_MD_OUT.relative_to(ROOT)).replace("\\", "/"))
    args = parser.parse_args()

    report = build_report(read_seed(SEED_MANIFEST), read_json(REACTION_KPIS))
    json_out = ROOT / args.json_out
    md_out = ROOT / args.md_out
    write_outputs(report, json_out, md_out)
    print(
        json.dumps(
            {
                "status": report["status"],
                "seed_rows": report["seed_summary"]["rows"],
                "blocked_family_count": report["blocked_family_count"],
                "blocked_without_seed": report["blocked_without_seed"],
                "ready_for_retraining_families": report["ready_for_retraining_families"],
                "json_out": str(json_out.relative_to(ROOT)).replace("\\", "/"),
                "md_out": str(md_out.relative_to(ROOT)).replace("\\", "/"),
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
