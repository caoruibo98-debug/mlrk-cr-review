# Generation 06 Review

## Frozen state

- Target tag: `generation_06`
- evaluated_commit: `12f5f2e1f2daa0c4f0c10ecb3b926e6511f4bbb3`
- Strict life-science application score: `4.43 / 5`
- Contract tests: `PASSED 15`
- Real food-glycoside panel: `6 / 6` expected products ranked in top 5.

## Trigger for this generation

After `generation_05` was pushed to the GitHub review mirror, contract tests failed in the mirror because a fresh checkout did not contain deployment artifacts under `outputs/modular/ltr/...` and did not have `outputs/modular/predictions/clean_rutin.json`.

This exposed a production-readiness gap: local tests could pass because old generated outputs existed on disk, while a new clone could not reproduce the same readiness state.

## Fixes in this generation

1. The deployment model JSON files and clean2 metric CSVs are now versionable despite the broader `outputs/` ignore rule.
2. `tools/freeze_generation.py` now includes `.gitignore` and required deployment artifacts in freeze manifests.
3. The prediction contract test now generates `clean_rutin.json` through the production CLI if the file is absent.

## Current position

This generation improves deployability and reviewer reproducibility. It still remains `internal_mvp_only` because external head-to-head benchmarks, wet-lab validation, and broader model-output evidence coverage remain incomplete.
