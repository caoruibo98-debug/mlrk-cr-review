# ML Ranking Kernel Production Candidate

This clone preserves the original `ml_ranking_kernel` code and adds a production-candidate shell around it.

## Validate readiness

From this directory:

```bash
python tools/validate_production_readiness.py
```

or:

```bash
python -m mlrk_prod.cli validate-readiness
```

Expected current status:

```text
internal_mvp_only
```

That status is intentional. External benchmark and wet-lab gates are still incomplete.

## Run contract tests

If pytest is installed:

```bash
python -m pytest -q
```

If pytest is not installed:

```bash
python tools/run_contract_tests.py
```

## Run one prediction

```bash
python -m mlrk_prod.cli predict --name rutin --topn 8
```

Prediction payloads include `biochem_quality`, `quality_summary`, and `interpretation_ready_top`. The ranking score is still not a wet-lab probability.

## Life-science appraisal

Run the internal, non-clinical life-science appraisal:

```bash
python tools/life_science_appraisal.py --run-panel --topn 10
```

The appraisal uses a real curated food-glycoside panel, structure-aware expected-product matching, full-InChIKey stereo-aware review fields, anti-cheat internal metrics, biochemical quality flags, and claim-boundary checks. The current strict core-panel score is `3.88 / 5`; this is an internal research-product score, not wet-lab validation.

Run the production challenge panel:

```bash
python scripts/life_science_appraisal.py --panel data/production_challenge_panel.csv --run-panel --topn 10 --out outputs/appraisal/challenge_appraisal.json
```

The current challenge-panel score is `2.81 / 5` with `8 / 22` strict top-5 hits. The failure taxonomy shows `14` expected products not generated and `8` benchmark-only top-5 hits. That result is used to track the remaining gap to broad food microbiome metabolite prediction.

Build the combined production scorecard:

```bash
python scripts/production_scorecard.py
```

Build the reaction-family KPI report:

```bash
python scripts/reaction_family_kpis.py
```

This writes `outputs/appraisal/reaction_family_kpis.json` and `.csv`, splitting each family into candidate-pool recall, strict top-5 recall, benchmark traceability, model-output evidence coverage, production gap, and recommended next action.

Check repository layout and reviewer entrypoints:

```bash
python scripts/repo_doctor.py
```

This writes `outputs/appraisal/repo_doctor.json`.

Generate the scorecard-backed model card:

```bash
python scripts/model_card.py
```

This writes `docs/MODEL_CARD.md` and `outputs/appraisal/model_card_summary.json`.

Generate the final release summary:

```bash
python scripts/release_summary.py
```

This writes `docs/RELEASE_SUMMARY.md` and `outputs/appraisal/release_summary.json`.

Record external AI/code-review status:

```bash
python scripts/external_review_status.py --tool CodeRabbit --status timed_out --command "coderabbit review --agent --base origin/master" --duration-ms 604046
```

This writes `outputs/appraisal/external_review_status.json`. A timed-out review is not a passed review.

Record fresh-clone release-candidate status:

```bash
python scripts/fresh_clone_report.py --source-url https://github.com/caoruibo98-debug/mlrk-cr-review.git --source-ref generation_18 --commit fda7f6c8eff8ccd99eb3a6b7275037c2cfcacd2c --clone-path D:\CRB\FoodGut\fresh_clone_generation_19 --contract-tests "PASSED 49 contract tests" --scorecard-status internal_mvp_only --core-top5 1.0 --challenge-top5 0.364 --readiness-status internal_mvp_only --repo-doctor-status ready --repo-doctor-checks 38 --model-card-status ready
```

This writes `outputs/appraisal/fresh_clone_report.json`.

Export external benchmark inputs:

```bash
python scripts/export_external_benchmarks.py
```

This writes BioTransformer, MicrobeRX, and GutBug-style benchmark input files under `outputs/external_benchmarks/`. These files are not external validation results; they are the case-level handoff for external tool execution.

Score imported external outputs:

```bash
python scripts/score_external_results.py
```

Without real external tool outputs, the scorecard remains `awaiting_external_outputs`. Product-output tools are scored by expected-product InChIKey recall; EC/enzyme tools are imported as evidence coverage only unless expected EC labels are later curated.

Frozen iterations are recorded under:

```text
freezes/generation_01/manifest.json
freezes/generation_02/manifest.json
freezes/generation_03/manifest.json
freezes/generation_04/manifest.json
freezes/generation_05/manifest.json
freezes/generation_06/manifest.json
freezes/generation_07/manifest.json
freezes/generation_08/manifest.json
freezes/generation_09/manifest.json
freezes/generation_10/manifest.json
freezes/generation_11/manifest.json
freezes/generation_12/manifest.json
freezes/generation_13/manifest.json
freezes/generation_14/manifest.json
freezes/generation_15/manifest.json
freezes/generation_16/manifest.json
freezes/generation_17/manifest.json
freezes/generation_18/manifest.json
freezes/generation_19/manifest.json
freezes/generation_20/manifest.json
```

## Retrain deployment models

```bash
python tools/train_deploy_models.py
```

This retrains `outputs/modular/ltr/models_clean/ltr_chem_{A,B,C,D}.json` from `clean_candidates.parquet`, then runs readiness validation.

## Start internal API

```bash
python -m mlrk_prod.cli serve --host 127.0.0.1 --port 8765
```

API endpoints:

```text
GET  /health
GET  /readiness
POST /api/jobs
GET  /api/jobs/{job_id}
GET  /api/jobs/{job_id}/result
```

The API returns stable error codes for invalid input, schema validation, missing jobs, pending jobs, worker failures, and timeouts. See:

```text
docs/API_CONTRACT.md
```

Submit a job with Python:

```python
import json
import urllib.request

body = json.dumps({"name": "rutin", "topn": 5}).encode("utf-8")
req = urllib.request.Request(
    "http://127.0.0.1:8765/api/jobs",
    data=body,
    headers={"Content-Type": "application/json"},
)
job = json.loads(urllib.request.urlopen(req).read().decode("utf-8"))
print(job)
```

If the server was started in the background, check:

```powershell
Get-Content runtime\api_server.json
```

Stop it with:

```powershell
Stop-Process -Id <pid>
```

## Added production-candidate assets

- `production_artifact_manifest.json`
- `mlrk_prod/`
- `tools/validate_production_readiness.py`
- `tools/run_contract_tests.py`
- `tools/train_deploy_models.py`
- `tests/`
- `docs/PRODUCTION_READINESS_REVIEW.md`
- `docs/WEB_APP_MVP_SPEC.md`
- `docs/SECURITY_AND_DEPLOYMENT_BOUNDARIES.md`
- `docs/EVALUATION_PROTOCOL.md`

## Current safe product shape

Internal research MVP for ranking and reviewing rule-generated metabolite candidates.

## Current unsafe product shape

Public consumer or clinical product.
