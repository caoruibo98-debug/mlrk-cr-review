# ML Ranking Kernel Production Candidate

This repository is a production-candidate shell for a food and gut-microbiome metabolite ranking model.

The model combines:

1. rule-based candidate generation,
2. chemistry-only learned-to-rank scoring,
3. evidence overlays from enzyme, microbe, and literature records,
4. biochemical quality filters,
5. explicit claim boundaries.

## Current Position

Current safe claim:

> Internal research MVP for prioritizing rule-generated food-polyphenol metabolite candidates, strongest on glycoside-to-aglycone transformations.

Current unsafe claim:

> Broad prediction of all food-derived gut microbial metabolites, strain-aware metabolism, clinical effects, consumer health recommendations, or wet-lab occurrence probabilities.

## Current Performance Snapshot

As of `generation_10`:

- Core curated food-glycoside panel: `6 / 6` strict top-5 hits, score `3.88 / 5`.
- Production challenge panel: `8 / 22` strict top-5 hits, score `2.81 / 5`.
- Challenge failure taxonomy: `14` expected products not generated, `8` benchmark-only top-5 hits.
- Internal ranking comparison: LTR_chem mean recall@5 `0.94`, above random `0.532`, EC-only `0.525`, and Tanimoto `0.716`.
- Readiness status: `internal_mvp_only`.

The core score dropped from `4.43 / 5` after generation 09 because the default evidence pool is no longer a local Windows-only path. External evidence can still be supplied with `MLRK_EVIDENCE_POOL`, but reproducible scorecard runs should not depend on hidden local files.

The challenge-panel result is intentional and important: it shows that the current system is not yet ready for broad food microbiome metabolite prediction.

## Repository Layout

```text
mlrk_prod/      Production-candidate API, CLI, schema, readiness, and quality wrappers.
modular/        Legacy modular candidate generation, deployment prediction, and LTR scripts.
src/            Original research/data-building utilities.
scripts/        GitHub-friendly command entrypoints that wrap production tools.
tools/          Internal production, appraisal, freezing, and test utilities.
tests/          Contract tests for requests, prediction payloads, appraisal matching, and quality filters.
data/           Curated core and challenge evaluation panels.
docs/           Model cards, readiness reviews, production ledger, and scientific positioning.
outputs/        Versioned deployment metrics/models plus small appraisal reports.
freezes/        Generation manifests with file hashes and score snapshots.
```

## Quickstart

Run contract tests:

```bash
python scripts/run_contract_tests.py
```

Validate readiness:

```bash
python scripts/validate_production_readiness.py
```

Run one prediction:

```bash
python -m mlrk_prod.cli predict --name rutin --topn 10
```

Optional evidence overlays can be enabled by setting `MLRK_EVIDENCE_POOL` to a CSV with substrate/product InChIKey and enzyme, microbe, PMID fields. Without it, predictions still rank candidates but evidence fields may be `-`.

Run the core panel:

```bash
python scripts/life_science_appraisal.py --run-panel --topn 10
```

Run the challenge panel:

```bash
python scripts/life_science_appraisal.py --panel data/production_challenge_panel.csv --run-panel --topn 10 --out outputs/appraisal/challenge_appraisal.json
```

Build the production scorecard:

```bash
python scripts/production_scorecard.py
```

Export external benchmark inputs:

```bash
python scripts/export_external_benchmarks.py
```

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

## Evaluation Model

The production scorecard separates:

1. internal comparison against random, EC-only, Tanimoto, and full-feature baselines,
2. core curated food-glycoside appraisal,
3. challenge panel appraisal across harder food-metabolism classes,
4. candidate-generation failure, ranking failure, evidence failure, and no-candidate states,
5. external comparison readiness and exported inputs for BioTransformer, MicrobeRX, GutBug, MIMOSA2, and AGREDA,
6. remaining production blockers.

The external benchmark input files are under `outputs/external_benchmarks/`. They are not external validation results; they are the reproducible handoff for running those tools and importing their outputs later.

## Iteration Process

The repository uses explicit generation freezes. Each production round must include:

1. a material improvement,
2. contract tests,
3. relevant appraisal or scorecard output,
4. a review note,
5. a git commit,
6. a freeze manifest,
7. a tag.

See:

- `docs/production/ITERATION_LEDGER.md`
- `docs/production/SCIENTIFIC_POSITIONING.md`
- `docs/reviews/`
- `freezes/`
