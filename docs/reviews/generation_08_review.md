# Generation 08 Review

## Production question

Can the challenge panel be expanded beyond the narrow glycoside niche, and can the evaluator distinguish candidate-generation failures from ranking failures?

## Changes

1. Expanded `data/production_challenge_panel.csv` from 6 to 22 literature-backed food and gut-microbiome metabolism cases.
2. Added quoted CSV fields for comma-containing metabolite names so the resolver evaluates the intended structures.
3. Added `candidate_summary` to prediction payloads, including the full generated candidate InChIKey block-1 pool and candidate source counts.
4. Added failure taxonomy to `tools/life_science_appraisal.py`.
5. Added failure-type aggregation to `tools/production_scorecard.py`.
6. Added regression tests for candidate-pool matching and generation-vs-ranking failure classification.

## Literature and external position

The challenge panel targets food-polyphenol gut microbial metabolism classes that appear repeatedly in the literature: hydroxycinnamate reduction, ellagitannin-to-urolithin metabolism, isoflavone reduction to equol-pathway intermediates, lignan enterolignan formation, flavanol/flavonol ring fission, and gallate decarboxylation.

External tools set a higher production bar:

- BioTransformer 3.0 combines machine learning and rule-based metabolism prediction, supports human gut microbial transformations, and reports predicted products with reaction and enzyme provenance.
- MicrobeRX is explicitly reaction-rule based and states that prediction depends on its reaction-rule and evidence databases.
- GutBug predicts gut bacterial enzyme EC numbers and reaction centers using ML/chemoinformatics and validates against known gut-bacterial biotransformations.

Against that landscape, this repository is still an internal research MVP. Its current value is a transparent, testable ranking shell around rule-generated candidates, not a broad gut-microbial metabolism predictor.

## Results

- Contract tests: `PASSED 18`
- Challenge panel resolver check: `22` rows, `0` unresolved substrate or product names.
- Readiness: `internal_mvp_only`
- Core panel score: `4.43 / 5`
- Core panel top-5 hit rate: `6 / 6`
- Core panel failure taxonomy: `4` hit_top5_model_evidence, `2` hit_top5_benchmark_only
- Challenge panel score: `1.79 / 5`
- Challenge panel top-5 hit rate: `1 / 22`
- Challenge panel expected product in candidate pool: `1 / 22`
- Challenge panel failure taxonomy:
  - `18` expected_product_not_generated
  - `3` no_generated_candidates
  - `1` hit_top5_benchmark_only

The only challenge-panel top-5 hit was `glycitin -> glycitein` at rank 1. No challenge-panel expected product had model-output evidence attached.

## Interpretation

This generation makes the main production blocker much sharper. The ranker is not the primary bottleneck on broad food gut-microbial metabolism. In 21 of 22 challenge cases, the expected product is either absent from the generated candidate pool or no candidates are generated at all. Ranking cannot recover products that candidate generation never emits.

The current model is therefore positioned as:

- usable for internal triage of rule-generated food-polyphenol candidates,
- strongest on glycoside-to-aglycone transformations,
- weak on multi-step reductive, dehydroxylation, ring-fission, decarboxylation, hydrolysis, and enterolignan/urolithin pathways,
- not production-ready for broad food-derived gut microbial metabolite prediction.

## Next round

Generation 09 should prioritize generation recall before ranking optimization. The most valuable next improvement is to add a controlled challenge-family candidate-generation layer for a small number of well-scoped transformations, then rerun the same 22-case panel without changing the scoring rubric.

Candidate families to consider first:

1. hydroxycinnamate double-bond reduction,
2. isoflavone first-step reduction,
3. simple glycoside/aglycone release beyond current flavonoid cases,
4. gallate decarboxylation,
5. conservative ellagitannin hydrolysis to ellagic acid.

The review gate for generation 09 should be candidate-pool recall, not final top-5 score alone.
