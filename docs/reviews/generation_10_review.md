# Generation 10 Review

## Production question

Can candidate-generation recall improve on the fixed 22-case challenge panel without changing the scoring rubric or using case-specific lookup?

## Anti-cheat rule

This generation did not add substrate-name or case-id lookup. The added rescue layer is structure-driven: it applies generic RDKit reaction templates to any compatible molecule and allows at most two transformation steps.

## Changes

1. Added `mlrk_prod/microbial_rescue.py` with generic food/gut microbial transformation templates:
   - hydroxycinnamate side-chain reduction,
   - stilbene double-bond reduction,
   - aromatic carboxyl decarboxylation,
   - phenol dehydroxylation.
2. Allowed up to two generic rescue steps so downstream products such as reduced-and-dehydroxylated stilbenes can be enumerated.
3. Integrated microbial rescue candidates into `modular/predict_substrate_clean.py` after legacy rules and aromatic O-glycoside rescue.
4. Added contract tests for:
   - caffeic acid to dihydrocaffeic acid,
   - gallic acid to pyrogallol,
   - resveratrol to lunularin through a two-step rescue.
5. Updated the no-candidate contract test from `ellagic acid` to `methane`, because ellagic acid now generates rescue candidates.

## Results

- Contract tests: `PASSED 22`
- Core panel score: `3.88 / 5`
- Core panel strict top-5 hit rate: `6 / 6`
- Challenge panel score: `2.81 / 5`
- Challenge panel strict top-5 hit rate: `8 / 22`
- Challenge expected product in candidate pool: `8 / 22`
- Challenge failure taxonomy:
  - `8` hit_top5_benchmark_only
  - `14` expected_product_not_generated

Generation 09 challenge state was `1 / 22` top-5, score `1.79 / 5`, with `18` expected products not generated and `3` no-candidate substrates. Generation 10 therefore improves the fixed challenge panel by `+7` top-5 hits and removes the no-candidate category for the current 22-case panel.

## Recovered challenge cases

1. `caffeic acid -> dihydrocaffeic acid`
2. `ferulic acid -> dihydroferulic acid`
3. `p-coumaric acid -> phloretic acid`
4. `resveratrol -> dihydroresveratrol`
5. `resveratrol -> lunularin`
6. `glycitin -> glycitein`
7. `gallic acid -> pyrogallol`
8. `urolithin C -> urolithin A`

## Interpretation

This is the first generation where the challenge panel improves materially through model logic rather than documentation or scoring changes. The improvement is still limited to reaction families with simple structural templates. It does not yet solve:

- isoflavone reductive pathways to dihydrodaidzein, dihydrogenistein, or equol,
- enterolignan formation from secoisolariciresinol,
- flavanol/flavonol/flavanone ring-fission chemistry,
- ellagitannin hydrolysis and ellagic-acid-to-urolithin multistep formation,
- gene, strain, abundance, or enzyme-context prediction.

## Next round

Generation 11 should add external benchmark export files for BioTransformer, MicrobeRX, and GutBug-style comparison, while preserving the fixed challenge panel. Generation 12 can then return to deeper reaction-family expansion with a better external comparison harness.
