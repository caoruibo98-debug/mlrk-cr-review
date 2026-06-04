# Round23 Source-backed Sample Expansion Pipeline

## Working conclusion

The current production blocker is primarily upstream of the ranker: important reaction families are missing from the candidate-generation/rule layer. The strongest current evidence is Round22: 14 P0 reactions were checked under 56 route variants and none generated the target product.

This does not mean sample size is irrelevant. It means adding positives alone will not help unless the source-backed reaction family can be generated as a deployment candidate.

## What changed in this round

- Added source-backed query raw outputs for Rhea and NCBI PubMed under `data/curation/rhea_round23_raw/` and `data/curation/ncbi_round23_raw/`.
- Screened ECReact in `runtime/external_round23/` and froze the P0 exact-pair summary in a tracked CSV.
- Created explicit gates separating positive samples, ranking decoys, condition-specific negatives, and unknowns.

## Key findings

- `dopamine -> m-tyramine` is exact in Rhea as RHEA:61520 and has gut catechol-dehydroxylase literature support, but the current generator still cannot produce the target.
- Urolithin dehydroxylation has recent gut Enterocloster literature support, but Rhea and ECReact did not return exact coverage in this round. It needs paper-table extraction and atom mapping before training.
- `protocatechuic acid -> catechol` is exact in Rhea and ECReact, but gut-specific evidence and prior PMID mismatch still need re-curation before promotion.
- Absence from a database is never a negative label. It is an unknown until a real assay-negative or condition-specific non-production record is found.

## Generated artifacts

* `data/curation/round23_current_blocker_audit.csv`
* `data/curation/round23_p0_source_verification.csv`
* `data/curation/round23_source_query_manifest.csv`
* `data/curation/round23_external_training_source_candidates.csv`
* `data/curation/round23_negative_sample_strategy.csv`
* `data/curation/round23_sample_expansion_pipeline.csv`
* `data/curation/round23_ecreact_p0_exact_pair_screen.csv`

## Production gate

No Round23 candidate is allowed directly into training. The next production-aligned step is to derive or import source-backed, atom-mapped reaction rules for the missing P0 families, run generation dry-runs, then only promote reactions whose targets are generated and whose labels survive false-negative screening.
