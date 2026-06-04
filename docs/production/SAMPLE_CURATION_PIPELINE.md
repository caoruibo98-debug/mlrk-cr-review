# Sample Curation Pipeline

## Purpose

This pipeline converts the sample-gap audit into training-ready data without inventing unsupported reactions.

The current model cannot move to production mainly because candidate generation and sample routing do not cover important gut microbial reaction families. The problem is not only raw sample shortage: many high-evidence reactions are already present in the local positive pool but are absent from the final `clean_candidates_full.parquet` table.

## Current Round

Round 2 artifacts:

- Positive candidate manifest: `data/curation/reaction_gap_positive_manifest_round2.csv`
- Negative strategy manifest: `data/curation/reaction_gap_negative_strategy_round2.csv`
- Build script: `tools/build_sample_gap_manifests.py`

Run:

```powershell
python tools\build_sample_gap_manifests.py
```

Round 2 output summary:

- Positive candidate rows: 75
- Families covered: 12
- `pipeline_recovery_candidate`: 36
- `positive_pool_to_reactions_candidate`: 8
- `already_in_clean_candidate_pool`: 8
- `do_not_train_exact_holdout`: 22
- `absent_from_local_positive_pool`: 1
- Negative strategy rows: 48

All Round 2 rows keep `training_allowed_round2=false`. This is intentional. A row becomes trainable only after source-level verification, holdout leakage review, structure verification, and candidate-generation recall review.

## Evidence Rules

Round 2 positive candidates require:

1. local positive-pool match to a target reaction family;
2. `evidence_level >= 3`;
3. PMID or DOI present;
4. exact substrate/product SMILES and InChIKey present.

This keeps weak database/rule-only rows out of training-positive promotion. MicrobeRX/RetroRules-style database rules can support candidate generation, but they are not automatically gold positives.

## Validation Layers

Use the layers in this order:

1. Local source alignment: positive pool, `reactions.parquet`, `clean_candidates_full.parquet`, and holdout manifest.
2. PubMed/DOI verification: confirm the cited paper exists and the title matches the claimed pathway.
3. Structure verification: spot-check InChIKey and SMILES against PubChem/ChEBI/HMDB.
4. Biological context verification: enzyme, EC, gene, microbe, strain, or fecal/community context.
5. Leakage review: exact benchmark/challenge holdout pairs must remain `do_not_train`.
6. Candidate-generation recall: verify the rule/candidate generator can produce the product before retraining the ranker.

## Positive Sample Promotion

Promotion statuses:

- `do_not_train_exact_holdout`: keep only for evaluation/diagnosis.
- `pipeline_recovery_candidate`: evidence exists in `reactions.parquet`, but the positive is missing from `clean_candidates_full`; fix routing or candidate generation before adding new data.
- `positive_pool_to_reactions_candidate`: evidence exists upstream but not in `reactions.parquet`; inspect `ltr_build.py` filters.
- `already_in_clean_candidate_pool`: use as a control, not as new evidence.
- `absent_from_local_positive_pool`: requires external curation before any model training.

## Negative Sample Policy

Do not treat absence as a true negative. Allowed negative types:

- `unlabeled_generated_nonmatch`: keep unlabeled; do not use as `y=0`.
- `class_contrast_negative`: usable only for reaction-family classification.
- `context_negative_nonproducer`: usable only when the model includes strain/community context.
- `hard_negative_after_source_review`: usable only after a reviewed source supports a different product under the same context.

## Next Iteration

The next iteration should prioritize:

1. Trace the 36 `pipeline_recovery_candidate` rows through candidate generation and clean-candidate filtering.
2. Fix or document why those positives are removed.
3. Externally curate `isoflavone_c_glycoside_conversion`, because it is the only Round 2 family with no local high-evidence positive-pool row after stricter filtering.
4. Add a separate `training_allowed=true` manifest only after the row passes all validation layers.
5. Recompute generation recall before any LTR retraining.
