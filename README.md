# ML Ranking Kernel Production Candidate

Production-candidate shell for an internal food and gut-microbiome metabolite ranking model.

This repository does not claim to be a broad, externally validated predictor. It is currently an internal research MVP for prioritizing rule-generated food-polyphenol metabolite candidates.

## Reviewer Summary

Current status: `internal_mvp_only`

Current safe claim:

> Internal research MVP for prioritizing rule-generated food-polyphenol metabolite candidates, strongest on glycoside-to-aglycone transformations.

Current unsafe claim:

> Broad prediction of all food-derived gut microbial metabolites, strain-aware metabolism, clinical effects, consumer health recommendations, or wet-lab occurrence probabilities.

As of `generation_19`:

- Core curated food-glycoside panel: `6 / 6` strict top-5 hits, score `3.88 / 5`.
- Production challenge panel: `8 / 22` strict top-5 hits, score `2.81 / 5`.
- Challenge failure taxonomy: `14` expected products not generated, `8` benchmark-only top-5 hits.
- Reaction-family KPI report: `19` families, separating candidate-generation-blocked families from evidence-integration-blocked families.
- API contract: stable error codes for invalid input, schema validation, pending jobs, failed jobs, and timeout states.
- Model card: generated from `outputs/appraisal/production_scorecard.json`.
- External review gate: CodeRabbit CLI/auth passed, review command timed out after `604046 ms`; no external review pass is claimed.
- Fresh-clone release-candidate check: `passed` on GitHub tag `generation_18` at commit `fda7f6c`.
- Repository doctor: `ready`, with `40` layout and entrypoint checks passing.
- Contract tests: `51`.
- Internal ranking comparison: LTR_chem mean recall@5 `0.94`, above random `0.532`, EC-only `0.525`, and Tanimoto `0.716`.

The challenge-panel result is intentional and important: it shows the model is not ready for broad food microbiome metabolite prediction.

## Fast Review Path

Install dependencies:

```bash
python -m pip install -r requirements.txt
```

Run the main reviewer checks:

```bash
python scripts/run_contract_tests.py
python scripts/validate_production_readiness.py
python scripts/repo_doctor.py
python scripts/model_card.py
python scripts/external_review_status.py --tool CodeRabbit --status timed_out --command "coderabbit review --agent --base origin/master" --duration-ms 604046
python scripts/production_scorecard.py
```

Expected high-level results:

```text
PASSED 51 contract tests
readiness status: internal_mvp_only
repo doctor: ready
scorecard status: internal_mvp_only
```

Run one prediction:

```bash
python -m mlrk_prod.cli predict --name rutin --topn 10
```

Optional evidence overlays can be enabled by setting `MLRK_EVIDENCE_POOL` to a CSV with substrate/product InChIKey and enzyme, microbe, PMID fields. Reproducible scorecard runs should not depend on hidden local evidence files.

## Current Evidence

Run the core panel:

```bash
python scripts/life_science_appraisal.py --run-panel --topn 10
```

Run the challenge panel:

```bash
python scripts/life_science_appraisal.py --panel data/production_challenge_panel.csv --run-panel --topn 10 --out outputs/appraisal/challenge_appraisal.json
```

Build reaction-family KPI reporting:

```bash
python scripts/reaction_family_kpis.py
```

Export external benchmark inputs:

```bash
python scripts/export_external_benchmarks.py
```

Score imported external benchmark outputs:

```bash
python scripts/score_external_results.py
```

The external benchmark input files under `outputs/external_benchmarks/` are not external validation results. They are the reproducible handoff for running BioTransformer, MicrobeRX, GutBug-style EC/enzyme comparison, and related external checks.

## Important Outputs

| Path | Meaning |
| --- | --- |
| `outputs/appraisal/production_scorecard.json` | Combined internal readiness and evaluation scorecard. |
| `outputs/appraisal/life_science_appraisal.json` | Core food-glycoside appraisal. |
| `outputs/appraisal/challenge_appraisal.json` | Harder production challenge panel appraisal. |
| `outputs/appraisal/reaction_family_kpis.json` | Reaction-family coverage and next-action report. |
| `outputs/appraisal/model_card_summary.json` | Scorecard-backed model-card summary. |
| `outputs/appraisal/external_review_status.json` | External AI/code-review attempt status. |
| `outputs/appraisal/fresh_clone_report.json` | Fresh-clone release-candidate verification report. |
| `outputs/appraisal/repo_doctor.json` | Repository layout and reviewer-entrypoint check report. |
| `outputs/external_benchmarks/manifest.json` | External benchmark export manifest. |
| `outputs/external_benchmarks/external_result_scorecard.json` | External result import status and scores, currently awaiting real external outputs. |
| `freezes/generation_*/manifest.json` | Frozen generation file hashes and review notes. |
| `docs/MODEL_CARD.md` | Scorecard-backed model card and claim boundary. |

## API

Start the internal API:

```bash
python -m mlrk_prod.cli serve --host 127.0.0.1 --port 8765
```

Endpoints:

```text
GET  /health
GET  /readiness
POST /api/jobs
GET  /api/jobs/{job_id}
GET  /api/jobs/{job_id}/result
```

The internal API exposes a stable error contract for web-app clients. See `docs/API_CONTRACT.md`.

## Repository Layout

```text
mlrk_prod/      Internal API, CLI, schema, readiness, and quality wrappers.
scripts/        Reviewer-facing commands that wrap matching tools.
tools/          Appraisal, scorecards, exports, freezing, and repository checks.
tests/          Contract tests runnable through scripts/run_contract_tests.py.
docs/           Production boundary, API contract, evaluation, and iteration docs.
data/           Curated core and challenge evaluation panels.
outputs/        Versioned metrics, model artifacts, benchmark exports, and appraisal reports.
freezes/        Generation manifests with file hashes and score snapshots.
modular/        Legacy candidate-generation, prediction, and LTR scripts.
src/            Legacy research/data-building utilities.
```

The layout policy is documented in `docs/REPOSITORY_GUIDE.md`: keep legacy scientific code in place, and add production wrappers, tests, contracts, and reports around it.

## Known Limitations

- External head-to-head outputs from BioTransformer, MicrobeRX, and GutBug-style enzyme tools have not been imported yet.
- No wet-lab validation is included.
- No strain-level abundance or genome context is connected to predictions.
- Candidate generation is still the ceiling: ranking cannot recover metabolites absent from the candidate pool.
- Model-output evidence coverage remains incomplete even when benchmark traceability exists.
- Public consumer, clinical, and health recommendation claims are out of scope.

## Reviewer Documents

- `docs/PRODUCTION_READINESS_REVIEW.md`
- `docs/EVALUATION_PROTOCOL.md`
- `docs/API_CONTRACT.md`
- `docs/MODEL_CARD.md`
- `docs/REPOSITORY_GUIDE.md`
- `docs/WEB_APP_MVP_SPEC.md`
- `docs/SECURITY_AND_DEPLOYMENT_BOUNDARIES.md`
- `docs/production/SCIENTIFIC_POSITIONING.md`
- `docs/production/ITERATION_LEDGER.md`
- `docs/reviews/`

## Iteration Process

Each production round must include:

1. one material improvement,
2. contract tests,
3. relevant appraisal or scorecard output,
4. a review note under `docs/reviews/`,
5. a git commit,
6. a freeze manifest,
7. a tag.
