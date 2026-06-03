# Generation 17 Review

## Production question

Can model cards and claim boundaries be tied to scorecard outputs?

## Changes

1. Added `tools/model_card.py` and `scripts/model_card.py`.
2. Generated `docs/MODEL_CARD.md` from `outputs/appraisal/production_scorecard.json`.
3. Generated `outputs/appraisal/model_card_summary.json`.
4. Added model-card contract tests.
5. Added model card checks to `repo_doctor`.
6. Updated README, repository guide, production manifest, and iteration ledger.

## Results

- Model card generation: `ready`
- Model-card source: `outputs/appraisal/production_scorecard.json`
- Repository doctor: `ready`
- Repository doctor checks: `36`
- Contract tests: `PASSED 47`
- Production scorecard status: `internal_mvp_only`
- Core strict top-5 hit rate: `6 / 6`
- Challenge strict top-5 hit rate: `8 / 22`
- Readiness: `internal_mvp_only`

## Audit interpretation

This generation reduces claim drift. The model card now reads the scorecard rather than relying only on manually copied prose.

The generated model card keeps these boundaries visible:

1. internal benchmark evidence only,
2. no imported external tool outputs yet,
3. no wet-lab validation,
4. scores are ranking signals over generated candidates, not occurrence probabilities.

## Next round

Generation 18 should run a fresh AI/code review pass on the production branch and address actionable findings, while avoiding repeated CodeRabbit triggers if rate limits are active.
