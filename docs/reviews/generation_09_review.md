# Generation 09 Review

## Production question

Can CodeRabbit's generation 08 review findings be fixed so the production-candidate evaluation is more reproducible and less likely to overstate biochemical identity?

## CodeRabbit findings addressed

1. Dependency declarations were unpinned, risking non-reproducible installs and appraisal drift.
2. InChIKey block-1 matching was connectivity-only and could treat stereochemistry mismatches as full expected-product hits.
3. Freeze manifests silently skipped missing tracked paths.
4. `glycoside_rescue.py` disabled RDKit logging globally at import time and rebuilt the same reaction on every call.
5. `modular/predict_substrate_clean.py` used a local Windows absolute path as the default evidence pool.
6. A no-candidate prediction contract test could fail with an unclear `FileNotFoundError`.
7. `production_scorecard.py` produced a raw missing-file traceback if internal metrics were absent.

## Changes

1. Pinned runtime, API, research, and test dependencies in both `pyproject.toml` and `requirements.txt`.
2. Added full-InChIKey metadata to the appraisal output:
   - `expected_inchikey`
   - `expected_matched_product_inchikey`
   - `expected_full_inchikey_match`
3. Added strict top-5 scoring that excludes explicit full-InChIKey mismatches while still reporting connectivity-level top-5 hits.
4. Added `hit_top5_connectivity_only_stereo_mismatch` and `hit_below_top5_connectivity_only_stereo_mismatch` failure classes.
5. Added scorecard fields for connectivity top-5 hits and stereo-mismatch top-5 hits.
6. Made freeze generation warn to stderr when a configured `TRACKED` path is missing.
7. Made the default evidence pool repository-relative. External evidence must now be provided explicitly through `MLRK_EVIDENCE_POOL`.
8. Scoped RDKit log suppression to the rescue reaction build/run path and precompiled the rescue reaction once.
9. Added regression coverage for stereo mismatch and clearer no-candidate output-file failure.

## Results

- Contract tests: `PASSED 19`
- Core panel score: `3.88 / 5`
- Core panel strict top-5 hit rate: `6 / 6`
- Core panel connectivity top-5 hit rate: `6 / 6`
- Core panel stereo-mismatch top-5 count: `0`
- Core panel model-output evidence coverage: `0 / 6`
- Challenge panel score: `1.79 / 5`
- Challenge panel strict top-5 hit rate: `1 / 22`
- Challenge panel connectivity top-5 hit rate: `1 / 22`
- Challenge panel stereo-mismatch top-5 count: `0`
- Challenge panel failure taxonomy:
  - `18` expected_product_not_generated
  - `3` no_generated_candidates
  - `1` hit_top5_benchmark_only

## Interpretation

This generation intentionally lowers the core appraisal score from `4.43 / 5` to `3.88 / 5` because the previous evidence overlay depended on a local unversioned evidence pool. The prediction/ranking top-5 behavior did not degrade; the score became more reproducible and less dependent on Ray's machine.

The model's current honest position is:

- structure/ranking behavior remains stable on the core glycoside panel,
- evidence provenance is not reproducible unless the evidence pool is versioned or explicitly supplied,
- strict identity scoring now has room to catch future stereo mismatches,
- challenge failures are still dominated by candidate-generation recall.

## Next round

Generation 10 should return to model capability. The immediate target should be candidate-generation recall on the 22-case challenge panel, with a strict rule that the score rubric and challenge panel remain fixed while generation logic changes.
