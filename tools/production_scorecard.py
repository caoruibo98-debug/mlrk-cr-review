from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from collections import Counter
from pathlib import Path
from statistics import mean
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from mlrk_prod.manifest import readiness_status, validate_readiness  # noqa: E402


METRICS = ROOT / "outputs" / "modular" / "ltr" / "clean2_metrics.csv"
CORE_REPORT = ROOT / "outputs" / "appraisal" / "life_science_appraisal.json"
CHALLENGE_REPORT = ROOT / "outputs" / "appraisal" / "challenge_appraisal.json"


def parse_mean(value: str) -> float:
    match = re.search(r"-?\d+(?:\.\d+)?", str(value))
    if not match:
        raise ValueError(f"cannot parse metric value {value!r}")
    return float(match.group(0))


def read_json(path: Path) -> dict[str, Any] | None:
    if not path.exists():
        return None
    return json.loads(path.read_text(encoding="utf-8"))


def internal_metric_summary() -> dict[str, Any]:
    rows = list(csv.DictReader(METRICS.open("r", encoding="utf-8-sig", newline="")))
    by_method: dict[str, list[float]] = {}
    by_module: dict[str, dict[str, float]] = {}
    for row in rows:
        method = row["method"]
        module = row["module"]
        r5 = parse_mean(row["r@5"])
        by_method.setdefault(method, []).append(r5)
        by_module.setdefault(module, {})[method] = r5

    means = {method: round(mean(values), 3) for method, values in sorted(by_method.items())}
    ltr = means.get("LTR_chem")
    deltas = {}
    if ltr is not None:
        for baseline in ("random", "ec_only", "tanimoto"):
            if baseline in means:
                deltas[f"LTR_chem_minus_{baseline}_r5"] = round(ltr - means[baseline], 3)
    if "LTR_full" in means and ltr is not None:
        deltas["LTR_full_minus_LTR_chem_r5"] = round(means["LTR_full"] - ltr, 3)

    weak_modules = []
    for module, scores in sorted(by_module.items()):
        module_ltr = scores.get("LTR_chem")
        if module_ltr is None:
            continue
        strongest_baseline = max(scores.get("random", 0.0), scores.get("ec_only", 0.0), scores.get("tanimoto", 0.0))
        if module_ltr <= strongest_baseline:
            weak_modules.append(module)

    return {
        "metric_file": str(METRICS.relative_to(ROOT)).replace("\\", "/"),
        "r5_mean_by_method": means,
        "ablation_deltas": deltas,
        "anti_cheat_pass": not weak_modules and deltas.get("LTR_chem_minus_ec_only_r5", 0.0) > 0,
        "modules_without_ltr_gain": weak_modules,
    }


def panel_summary(report: dict[str, Any] | None, label: str) -> dict[str, Any]:
    if report is None:
        return {"label": label, "status": "missing_report"}
    cases = report.get("cases", [])
    top5 = [c for c in cases if c.get("expected_rank") is not None and int(c["expected_rank"]) <= 5]
    generated = [c for c in cases if c.get("expected_rank") is not None]
    expected_in_pool = [c for c in cases if c.get("expected_in_candidate_pool") is True]
    evidence = Counter(str(c.get("expected_evidence_source", "none")) for c in cases)
    families = Counter(str(c.get("reaction_family", "unknown")) for c in cases)
    failure_types = Counter(str(c.get("failure_type", "not_recorded")) for c in cases)
    failures = [
        {
            "case_id": c.get("case_id"),
            "substrate": c.get("substrate_name"),
            "expected_product": c.get("expected_product_name"),
            "reaction_family": c.get("reaction_family"),
            "failure_type": c.get("failure_type"),
            "expected_in_candidate_pool": c.get("expected_in_candidate_pool"),
            "expected_rank": c.get("expected_rank"),
            "match_type": c.get("expected_match_type"),
            "top_product": c.get("top_product"),
        }
        for c in cases
        if c.get("expected_rank") is None or int(c["expected_rank"]) > 5
    ]
    return {
        "label": label,
        "status": report.get("status"),
        "panel_file": report.get("panel_file"),
        "case_count": len(cases),
        "score": report.get("score", {}),
        "expected_generated_count": len(generated),
        "expected_in_candidate_pool_count": len(expected_in_pool),
        "top5_hit_count": len(top5),
        "top5_hit_rate": round(len(top5) / max(1, len(cases)), 3),
        "evidence_source_counts": dict(sorted(evidence.items())),
        "failure_type_counts": dict(sorted(failure_types.items())),
        "reaction_family_counts": dict(sorted(families.items())),
        "failure_cases": failures,
    }


def external_comparison_matrix() -> list[dict[str, str]]:
    return [
        {
            "tool_or_model": "BioTransformer 3/4",
            "comparison_role": "small-molecule biotransformation product generation",
            "current_status": "not_executed",
            "required_next_step": "Export core and challenge panels as SMILES inputs; compare generation recall and enzyme annotations in gut microbial/SuperBio modes.",
            "source": "Nucleic Acids Research 2022 BioTransformer 3.0; BioTransformer 4.0 reported as successor in 2025 records.",
        },
        {
            "tool_or_model": "MicrobeRX",
            "comparison_role": "enzyme-reaction-based human/gut metabolite prediction from GEM-derived reaction rules",
            "current_status": "not_executed",
            "required_next_step": "Run the same substrate SMILES through MicrobeRX and compare product InChIKey block-1 recall plus reaction evidence fields.",
            "source": "MicrobeRX paper and documentation describe metabolite SMILES, reaction identifiers, EC, Rhea, KEGG, and PubMed fields.",
        },
        {
            "tool_or_model": "GutBug",
            "comparison_role": "gut bacterial enzyme prediction for biotic/xenobiotic molecules",
            "current_status": "not_executed",
            "required_next_step": "Use as an enzyme/gene plausibility comparator rather than a strict product-ranking comparator if product output is not directly aligned.",
            "source": "GutBug Journal of Molecular Biology 2023 description.",
        },
        {
            "tool_or_model": "MIMOSA2 / AGORA / AGREDA",
            "comparison_role": "community metabolic potential and diet-compound network coverage",
            "current_status": "scope_reference_only",
            "required_next_step": "Do not compare rank@k directly; use these to audit pathway coverage, gene/reaction evidence, and sample-level microbiome context requirements.",
            "source": "MIMOSA2 Bioinformatics 2022; AGREDA Nature Communications 2021.",
        },
    ]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", default="outputs/appraisal/production_scorecard.json")
    args = parser.parse_args()

    issues = validate_readiness()
    report = {
        "status": readiness_status(issues),
        "readiness_issues": [issue.__dict__ for issue in issues],
        "internal_comparison": internal_metric_summary(),
        "core_panel": panel_summary(read_json(CORE_REPORT), "core_food_glycoside_panel"),
        "challenge_panel": panel_summary(read_json(CHALLENGE_REPORT), "production_challenge_panel"),
        "external_comparison_matrix": external_comparison_matrix(),
        "current_claim": "Internal research MVP for food-polyphenol candidate generation/ranking, strongest on glycoside aglycone release.",
        "remaining_gap_to_full_food_microbiome_metabolite_prediction": [
            "Challenge-panel failures are dominated by candidate-generation recall, not ranking.",
            "Reaction-family coverage is still narrow for reductions, dehydroxylations, ring fission, decarboxylation, and multi-step gut microbial pathways.",
            "External tools have not yet been run on the exact same substrate panel.",
            "No gene/genome or strain-level abundance context is connected to predictions.",
            "Candidate generation still constrains the ceiling; ranking cannot recover products absent from the candidate set.",
            "Model-output evidence coverage remains incomplete even when benchmark traceability is present.",
        ],
    }

    out = ROOT / args.out
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps({"status": report["status"], "out": str(out), "core_top5": report["core_panel"].get("top5_hit_rate"), "challenge_top5": report["challenge_panel"].get("top5_hit_rate")}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
