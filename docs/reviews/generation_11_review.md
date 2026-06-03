# Generation 11 Review

## Production question

Can external benchmark comparison move from a text-only plan to reproducible input artifacts for BioTransformer, MicrobeRX, and GutBug-style evaluation?

## Changes

1. Added `tools/export_external_benchmarks.py` and `scripts/export_external_benchmarks.py`.
2. Exported all `28` benchmark cases:
   - `6` core food-glycoside cases,
   - `22` production challenge cases.
3. Added `outputs/external_benchmarks/` with:
   - `benchmark_panel_unified.csv`,
   - `biotransformer_input.tsv`,
   - `biotransformer_input.sdf`,
   - `microberx_query_table.csv`,
   - `gutbug_query_table.csv`,
   - `gutbug_pubchem_cids.txt`,
   - `manifest.json`.
4. Added contract tests for external benchmark export presence and schema.
5. Added external benchmark export status to the production scorecard.
6. Updated `.gitignore` and freeze tracking so external benchmark inputs are versioned.

## External comparison framing

- BioTransformer is treated as a product-generation comparator. The exported TSV/SDF preserve `case_id` and substrate structure.
- MicrobeRX is treated as a reaction-rule and EC/Rhea/PubMed evidence comparator. The exported query table includes substrate and expected-product SMILES/InChIKey labels for recall scoring after MicrobeRX runs.
- GutBug is treated as an enzyme/EC plausibility comparator rather than a strict product-ranking comparator. PubChem CIDs are exported when resolvable.

## Results

- Export case count: `28`
- Panel counts: `6` core, `22` challenge
- Contract tests: `PASSED 26`
- Readiness: `internal_mvp_only`
- Core strict top-5 hit rate: `6 / 6`
- Challenge strict top-5 hit rate: `8 / 22`
- External comparison matrix:
  - BioTransformer 3/4: `input_exported`
  - MicrobeRX: `input_exported`
  - GutBug: `input_exported`
  - MIMOSA2 / AGORA / AGREDA: `scope_reference_only`

## Interpretation

This generation does not claim external validation. It creates the reproducible handoff needed to run external tools and import their results later. That matters because external comparison should be case-level and structure-level, not a narrative comparison in a README.

The model remains at generalization evidence level 3 at best for internal tests and rule coverage. It has not reached level 4 external validation because the external tools have not yet been executed and scored on the exported inputs.

## Next round

Generation 12 should add an external-result import and scoring harness. It should accept tool outputs, map them back to `case_id`, compare product InChIKey block-1/full InChIKey where available, and report external recall side by side with the current model.
