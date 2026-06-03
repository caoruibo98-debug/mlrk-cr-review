# Generation 12 Review

## Production question

Can external tool outputs be imported and scored case-by-case without pretending external validation is already complete?

## Changes

1. Added `tools/score_external_results.py` and `scripts/score_external_results.py`.
2. Added result templates:
   - `outputs/external_benchmarks/result_templates/external_product_results_template.csv`
   - `outputs/external_benchmarks/result_templates/external_enzyme_results_template.csv`
3. Added `outputs/external_benchmarks/external_result_scorecard.json`.
4. Added test fixtures for product-output and EC/enzyme-output external result examples.
5. Added contract tests for:
   - product result import,
   - case_id and InChIKey-based recall scoring,
   - enzyme result import without mislabeling it as product recall,
   - default `awaiting_external_outputs` status.
6. Added external-result scorecard integration to the production scorecard.
7. Updated readiness messaging to distinguish exported inputs plus import harness from completed external head-to-head validation.

## Scoring behavior

For product-output tools such as BioTransformer and MicrobeRX-style exports, the harness computes:

- top-1 block-1 recall,
- top-5 block-1 recall,
- any-rank block-1 recall,
- any-rank full-InChIKey recall,
- per-case hit rank and prediction count.

For EC/enzyme tools such as GutBug-style exports, the harness imports EC/enzyme/microbe rows and reports coverage only. It does not score product recall because the current benchmark does not have expected EC labels.

## Results

- Contract tests: `PASSED 29`
- External export case count: `28`
- External result scorecard status: `awaiting_external_outputs`
- Product tools scored from real external output: none yet
- Readiness: `internal_mvp_only`
- Readiness warning: external benchmark inputs and import harness are ready, but real external outputs are not scored yet.
- Core strict top-5 hit rate: `6 / 6`
- Challenge strict top-5 hit rate: `8 / 22`

## Interpretation

This generation adds the missing bridge between exported benchmark inputs and real external comparison metrics. It still does not claim external validation. The correct claim is now:

> The repository can export benchmark inputs and import/score external tool outputs, but no real BioTransformer, MicrobeRX, or GutBug output has been imported yet.

This is a necessary production step because it prevents future external comparisons from becoming prose-only or cherry-picked.

## Next round

Generation 13 should add substrate/reaction-family KPI reporting so the model's coverage can be read as a product metric: which reaction families are solved, which are candidate-generation failures, and which are ranking/evidence failures.
