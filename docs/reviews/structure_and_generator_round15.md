# Round15 Structure And Generator Audit

Date: 2026-06-04

## Direct Answer

Round15 makes the diagnosis sharper:

> The four Rhea source-backed candidates from Round14 are not simply new samples to add. All four already exist in `reactions.parquet`. The production question is whether they survive the clean generator path.

Local availability counts: `{'already_in_reactions': 4, 'already_in_clean_candidates_full': 3, 'already_in_clean_candidates': 1}`

Decision counts: `{'do_not_import_existing_positive_dropped_after_clean_full': 2, 'do_not_import_existing_positive_fix_generator_path': 1, 'do_not_import_existing_positive_already_in_clean_candidates': 1}`

No Round15 row is allowed directly into training.

## Structure Mapping

ChEBI was used as the primary structure source and PubChem as an InChIKey block cross-check. Primary substrate/product mappings are consistent at block1 level where PubChem was queried.

Structure table:

- `data/curation/round15_compound_structure_mapping.csv`

## Generator Recall

- `RHEA:16309` `taurochenodeoxycholate -> chenodeoxycholate`: sampled=True, full_pool=True, fullrule_hit=False
- `RHEA:19353` `glycocholate -> cholate`: sampled=True, full_pool=True, fullrule_hit=True
- `RHEA:20689` `chlorogenate -> trans-caffeate`: sampled=False, full_pool=False, fullrule_hit=True
- `RHEA:69683` `daidzein 7-O-beta-D-glucoside -> daidzein(1-)`: sampled=True, full_pool=True, fullrule_hit=True

Interpretation:

- If a pair is in `reactions.parquet` but absent from `clean_candidates`, duplicate-importing it as a new positive is wrong.
- If full or positive-pool rules can generate it but clean candidates do not contain it, the next action is generator/rule wiring or clean-path trace.
- If the pair reaches `clean_candidates_full` but not `clean_candidates`, the loss is a filtering/sampling/evaluation-path issue.

## Training Decisions

- `RHEA:16309` `taurochenodeoxycholate -> chenodeoxycholate`: `do_not_import_existing_positive_dropped_after_clean_full`
- `RHEA:19353` `glycocholate -> cholate`: `do_not_import_existing_positive_dropped_after_clean_full`
- `RHEA:20689` `chlorogenate -> trans-caffeate`: `do_not_import_existing_positive_fix_generator_path`
- `RHEA:69683` `daidzein 7-O-beta-D-glucoside -> daidzein(1-)`: `do_not_import_existing_positive_already_in_clean_candidates`

## Next Actions

- P1 `trace_clean_full_to_clean_candidates_loss` (2 pair): BHTRKEVKTKCXOH__RUDATBOHQWOJDD;RFDAIACWWDREDC__BHQCQFFYRZLCQQ
- P2 `promote_or_wire_source_traceable_rule_then_rerun_clean_candidates` (1 pair): CWVRJTMFETXNAD__QAIPRVGONGVQAS
- P4 `provenance_enrichment_or_source_manifest_only` (1 pair): KYQZWONCHDNPDP__ZQSIJRDFPHDXIC

## Files Created

- `data/curation/round15_compound_structure_mapping.csv`
- `data/curation/round15_candidate_pair_gate.csv`
- `data/curation/round15_generator_recall_check.csv`
- `data/curation/round15_training_decisions.csv`
- `data/curation/round15_next_actions.csv`
- `data/curation/chebi_round15_raw/*.json`
- `data/curation/pubchem_round15_raw/*.json`

## Claim Boundary

Allowed after Round15:

> Round15 confirmed that the strongest Rhea candidates are already represented in the normalized reaction table, and that production progress depends on generator recall/path retention rather than duplicate positive import.

Forbidden:

> Round15 expanded the training set or proved production readiness.
