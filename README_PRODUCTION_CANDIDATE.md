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
