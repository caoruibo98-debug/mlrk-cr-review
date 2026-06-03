from __future__ import annotations

import argparse
import csv
import json
import re
import subprocess
import sys
from pathlib import Path
from statistics import mean

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from mlrk_prod.biosanity import annotate_prediction_payload  # noqa: E402
from mlrk_prod.manifest import readiness_status, validate_readiness  # noqa: E402


PANEL = ROOT / "data" / "real_biochemistry_panel.csv"
PREDICTIONS = ROOT / "outputs" / "modular" / "predictions"
METRICS = ROOT / "outputs" / "modular" / "ltr" / "clean2_metrics.csv"


def safe_label(label: str) -> str:
    return "".join(c if c.isalnum() or c in "-_" else "_" for c in label)


def run_prediction(name: str, topn: int) -> None:
    cmd = [sys.executable, "-m", "mlrk_prod.cli", "predict", "--name", name, "--topn", str(topn)]
    proc = subprocess.run(
        cmd,
        cwd=str(ROOT),
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
        check=False,
        timeout=240,
    )
    if proc.returncode != 0:
        raise RuntimeError(f"prediction failed for {name}: {proc.stderr or proc.stdout}")


def read_panel() -> list[dict[str, str]]:
    with PANEL.open("r", encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f))


def read_prediction(name: str) -> dict:
    path = PREDICTIONS / f"clean_{safe_label(name)}.json"
    with path.open("r", encoding="utf-8") as f:
        return json.load(f)


def parse_mean(value: str) -> float:
    match = re.search(r"-?\d+(?:\.\d+)?", str(value))
    if not match:
        raise ValueError(f"cannot parse numeric mean from {value!r}")
    return float(match.group(0))


def metric_summary() -> dict:
    rows = list(csv.DictReader(METRICS.open("r", encoding="utf-8-sig", newline="")))
    ltr = [parse_mean(r["r@5"]) for r in rows if r["method"] == "LTR_chem"]
    ec = [parse_mean(r["r@5"]) for r in rows if r["method"] == "ec_only"]
    tan = [parse_mean(r["r@5"]) for r in rows if r["method"] == "tanimoto"]
    return {
        "ltr_r5_mean": round(mean(ltr), 3),
        "ec_only_r5_mean": round(mean(ec), 3),
        "tanimoto_r5_mean": round(mean(tan), 3),
        "anti_cheat_pass": mean(ltr) > mean(ec) and mean(ltr) > mean(tan),
    }


def find_expected_rank(payload: dict, expected: str) -> int | None:
    target = expected.strip().lower()
    for row in payload.get("top", []):
        if str(row.get("product_name", "")).strip().lower() == target:
            return int(row["rank"])
    return None


def score_appraisal(cases: list[dict], metrics: dict, readiness: str) -> dict:
    rank_hits = [c for c in cases if c["expected_rank"] is not None and c["expected_rank"] <= 5]
    high_quality_top = [c for c in cases if c["top_quality_tier"] in {"high", "medium"}]
    evidence_hits = [c for c in cases if c["expected_has_evidence"]]
    flagged_bad_top = [c for c in cases if c["top_quality_tier"] == "reject"]

    real_case_score = len(rank_hits) / max(1, len(cases))
    quality_score = max(0.0, (len(high_quality_top) - len(flagged_bad_top)) / max(1, len(cases)))
    evidence_score = len(evidence_hits) / max(1, len(cases))
    metric_score = 1.0 if metrics["anti_cheat_pass"] and metrics["ltr_r5_mean"] >= 0.80 else 0.55
    deployment_score = 0.65 if readiness == "internal_mvp_only" else 0.40

    dimensions = {
        "internal_ranker_and_anti_cheat": round(metric_score, 3),
        "real_biochemistry_panel": round(real_case_score, 3),
        "candidate_biochemical_sanity": round(quality_score, 3),
        "mechanistic_evidence_coverage": round(evidence_score, 3),
        "deployment_and_claim_boundaries": round(deployment_score, 3),
    }
    final = round(5.0 * mean(dimensions.values()), 2)
    return {"score_0_to_5": final, "dimensions": dimensions}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-panel", action="store_true")
    parser.add_argument("--topn", type=int, default=10)
    parser.add_argument("--out", default="outputs/appraisal/life_science_appraisal.json")
    args = parser.parse_args()

    panel = read_panel()
    if args.run_panel:
        for case in panel:
            run_prediction(case["substrate_name"], args.topn)

    cases = []
    for case in panel:
        payload = annotate_prediction_payload(read_prediction(case["substrate_name"]))
        expected_rank = find_expected_rank(payload, case["expected_product_name"])
        top = payload.get("top", [{}])[0] if payload.get("top") else {}
        expected = next(
            (row for row in payload.get("top", []) if str(row.get("product_name", "")).lower() == case["expected_product_name"].lower()),
            {},
        )
        cases.append(
            {
                **case,
                "module": payload.get("module_name"),
                "n_rule_candidates": payload.get("n_rule_candidates"),
                "expected_rank": expected_rank,
                "expected_has_evidence": "pmid=" in str(expected.get("evidence", "")).lower() or "enzyme=" in str(expected.get("evidence", "")).lower(),
                "top_product": top.get("product_name"),
                "top_quality_tier": top.get("biochem_quality", {}).get("tier"),
                "top_quality_flags": [f["code"] for f in top.get("biochem_quality", {}).get("flags", [])],
            }
        )

    issues = validate_readiness()
    readiness = readiness_status(issues)
    metrics = metric_summary()
    score = score_appraisal(cases, metrics, readiness)
    report = {
        "status": readiness,
        "score": score,
        "metrics": metrics,
        "cases": cases,
        "known_readiness_issues": [issue.__dict__ for issue in issues],
        "claim": "Internal research score only; not a wet-lab or clinical validation score.",
    }
    out = ROOT / args.out
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(report["score"], indent=2))
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
