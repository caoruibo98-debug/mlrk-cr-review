# Repository Guide

This repository is a production-candidate shell around the legacy food and gut-microbiome metabolite ranking kernel.

The standardization rule is conservative:

> Keep legacy scientific code in place. Add production wrappers, contracts, tests, and reports around it.

## Main directories

| Path | Role |
| --- | --- |
| `mlrk_prod/` | Internal API, CLI, schema, readiness, and quality wrappers. |
| `scripts/` | Reviewer-facing commands. Each script wraps a matching implementation under `tools/`. |
| `tools/` | Implementation code for appraisal, scorecards, exports, freezing, and repository checks. |
| `tests/` | Contract tests runnable without relying on pytest fixtures. |
| `docs/` | Production boundary, API contract, evaluation, web-app, and iteration documentation. |
| `data/` | Curated core and challenge evaluation panels. |
| `outputs/appraisal/` | Small versioned appraisal and scorecard reports. |
| `outputs/external_benchmarks/` | External-tool input exports and external result import templates. |
| `outputs/modular/ltr/` | Versioned deployment metrics and clean model artifacts. |
| `freezes/` | Generation manifests with file hashes. |
| `modular/` | Legacy candidate-generation, prediction, and training scripts. Do not move casually. |
| `src/` | Legacy research/data-building utilities. |

## Canonical commands

```bash
python scripts/run_contract_tests.py
python scripts/validate_production_readiness.py
python scripts/production_scorecard.py
python scripts/reaction_family_kpis.py
python scripts/model_card.py
python scripts/release_summary.py
python scripts/export_external_benchmarks.py
python scripts/score_external_results.py
python scripts/external_review_status.py
python scripts/fresh_clone_report.py
python scripts/repo_doctor.py
python -m mlrk_prod.cli predict --name rutin --topn 10
python -m mlrk_prod.cli serve --host 127.0.0.1 --port 8765
```

## Entry point policy

Reviewer-facing commands should live under `scripts/`.

Implementation logic should live under `tools/` or `mlrk_prod/`.

Legacy scientific logic should remain under `modular/` and `src/` unless a dedicated migration round moves it with tests and compatibility wrappers.

## Repository doctor

Run:

```bash
python scripts/repo_doctor.py
```

The doctor writes:

```text
outputs/appraisal/repo_doctor.json
```

It checks required docs, entrypoint pairs, legacy boundaries, and versioned outputs. It is a maintainability check, not a biological validation score.
