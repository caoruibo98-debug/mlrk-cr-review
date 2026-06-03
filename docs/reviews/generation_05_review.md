# Generation 05 Review

## Frozen state

- Target tag: `generation_05`
- evaluated_commit: `dd38d9977932151b807463225a0f41917441b163`
- Strict life-science application score: `4.43 / 5`
- Contract tests: `PASSED 15`
- Real food-glycoside panel: `6 / 6` expected products ranked in top 5.

## CodeRabbit PR review result

CodeRabbit completed the PR review for `generation_01..generation_03` on GitHub PR #1:

- Repository: `caoruibo98-debug/mlrk-cr-review`
- PR: `Review production candidate generations 01 to 03`
- Review state: `COMMENTED`
- Check status: `SUCCESS`

The actionable findings were treated as production-readiness blockers for the appraisal path.

## Fixes in this generation

1. `find_expected_row` now requires InChIKey block-1 structural agreement whenever an expected block is available. This prevents a product with the right display name but wrong structure from being counted as a hit.
2. Invalid or unparsable product SMILES now return `None` from `block1_from_smiles` instead of crashing the appraisal.
3. Rejected candidate ranks are parsed defensively, so malformed rank values do not become fake rank `0` entries or raise during payload annotation.
4. `generation_02` traceability now separates `evaluated_commit` from `freeze_commit` in both the review note and manifest.

## Regression coverage

New contract tests cover:

1. Structure-first matching rejects name-only false positives.
2. Invalid SMILES are tolerated by the evaluator.
3. Malformed rejected ranks are skipped in `quality_summary`.

## Current position

This generation is stronger as an internal MVP and reviewable product candidate. It is still not production-ready for public or clinical-facing claims because model-output evidence coverage is `0.667`, external head-to-head benchmark coverage is still incomplete, and wet-lab validation is absent.
