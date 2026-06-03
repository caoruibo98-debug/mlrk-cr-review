# Branch Strategy

## Stable Line

- Stable branch/tag baseline: `master` at `generation_20`.
- `master` should stay frozen unless a feature branch has passed data audit, contract tests, scorecard checks, and human review.
- New model-development work should start from the latest reviewed branch or from `generation_20`, then merge intentionally.

## Active Branches

| Branch | Base | Purpose | Merge condition |
| --- | --- | --- | --- |
| `biochem-logic-review` | `generation_20` | Tighten biochemical evidence tiers so structural-only candidates cannot be called high-evidence outputs. | Passed contract tests; should be reviewed before merging to stable. |
| `reaction-family-expansion` | `biochem-logic-review` | Resolve candidate-generation-blocked reaction families into machine-readable holdout seeds and retraining targets. | Do not merge as a retraining branch until independent non-holdout samples are collected. |

## Merge Policy

1. Keep challenge benchmark rows out of training.
2. Add independent training rows in a data branch before changing candidate-generation rules.
3. Rebuild evaluation artifacts only after leakage checks pass.
4. Run `python scripts/run_contract_tests.py`, `python scripts/repo_doctor.py`, and `python scripts/validate_production_readiness.py`.
5. Commit and freeze the version before merge.

## Current Reaction-Family Expansion Status

The `reaction-family-expansion` branch records `14` blocked benchmark seed rows across `12` reaction families. They are machine-readable but marked `training_allowed=false`. The branch therefore improves traceability and planning, not model training yet.
