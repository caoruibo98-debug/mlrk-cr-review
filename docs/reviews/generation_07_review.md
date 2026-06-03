# Generation 07 Review

## Production question

Can the model expose the next real production gap beyond the already-stable glycoside hydrolysis panel?

## Changes

1. Added `data/production_challenge_panel.csv` with six literature-backed food microbiome metabolism cases outside the narrow core glycoside panel.
2. Added `--panel` support to `tools/life_science_appraisal.py` so core and challenge panels can be run through the same evaluator.
3. Added `tools/production_scorecard.py` and `scripts/production_scorecard.py` to summarize internal comparison, ablation deltas, core panel, challenge panel, external comparison readiness, and remaining blockers.
4. Fixed prediction behavior when candidate generation returns no products. The CLI now writes a schema-compliant payload with `candidate_generation_status=no_candidates` instead of letting the production runner fail.
5. Added a regression test for the no-candidate prediction path.
6. Added GitHub-friendly `scripts/` entrypoints while preserving legacy `tools/` paths.
7. Replaced the garbled root `README.md` with a reviewer-first production-candidate README.
8. Expanded dependency declarations in `pyproject.toml` and added `requirements.txt`.

## Results

- Contract tests: `PASSED 16`
- Readiness: `internal_mvp_only`
- Core panel score: `4.43 / 5`
- Core panel top-5 hit rate: `6 / 6`
- Challenge panel score: `1.79 / 5`
- Challenge panel top-5 hit rate: `0 / 6`
- Internal LTR_chem recall@5 mean: `0.94`
- LTR_chem delta over EC-only recall@5: `+0.415`
- LTR_chem delta over Tanimoto recall@5: `+0.224`

## Interpretation

This generation makes the model more honest and more production-reviewable. The core model is stable for the current food-glycoside niche, but the challenge panel shows that the system does not yet solve broad food microbiome metabolite prediction.

The strongest blocker is candidate generation coverage: ranking cannot recover expected metabolites that are absent from the candidate set. The challenge failures include phenolic ester hydrolysis, ellagic-acid-to-urolithin formation, daidzein-to-equol conversion, and lignan enterolignan metabolism.

## Next round

Generation 08 should expand the challenge panel toward at least 20 literature-backed cases and add a failure taxonomy that separates:

1. no generated candidates,
2. expected product absent from generated candidates,
3. expected product generated but ranked too low,
4. expected product generated but rejected by biochemical quality filters,
5. expected product generated with missing evidence.
