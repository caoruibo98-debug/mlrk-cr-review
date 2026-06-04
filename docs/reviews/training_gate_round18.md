# Round18 Training Gate Audit

## Working conclusion

No new sample is allowed into training yet. Round18 turns Round17 leads into gates and keeps `training_allowed_round18=False` across the board.

The main reason is not lack of raw rows. It is missing production-grade evidence routing:

1. P0 families still have route blockers.
2. P1 families need family-balanced gold evaluation before retraining.
3. External Rhea/PubMed leads need exact-pair extraction, ChEBI/PubChem mapping, local de-duplication, generator-route checks, and leakage-safe split assignment.
4. Negative labels remain the biggest risk area: hard decoys and assay negatives must not be mixed.

## P0/P1 gate snapshot

```text
                          priority_family priority_level                  route_gate              source_exactness_gate                                       negative_gate                                       decision
bile_acid_deconjugation_and_lipid_context             P0                     blocked      no_round18_external_candidate ranking_decoys_conditional; assay_negatives_missing                 do_not_retrain_fix_route_first
      bile_acid_secondary_transformations             P1 pass_for_route_but_eval_gap needs_manual_exact_pair_extraction ranking_decoys_conditional; assay_negatives_missing do_not_retrain_build_family_balanced_gold_eval
        glucuronide_sulfate_deconjugation             P1 pass_for_route_but_eval_gap needs_manual_exact_pair_extraction ranking_decoys_conditional; assay_negatives_missing do_not_retrain_build_family_balanced_gold_eval
hydroxycinnamate_reduction_and_hydrolysis             P0                     blocked      no_round18_external_candidate ranking_decoys_conditional; assay_negatives_missing                 do_not_retrain_fix_route_first
                  polyphenol_ring_fission             P1 pass_for_route_but_eval_gap      no_round18_external_candidate ranking_decoys_conditional; assay_negatives_missing do_not_retrain_build_family_balanced_gold_eval
             tryptophan_indole_metabolism             P1 pass_for_route_but_eval_gap needs_manual_exact_pair_extraction ranking_decoys_conditional; assay_negatives_missing do_not_retrain_build_family_balanced_gold_eval
         tma_choline_carnitine_metabolism             P1 pass_for_route_but_eval_gap needs_manual_exact_pair_extraction ranking_decoys_conditional; assay_negatives_missing do_not_retrain_build_family_balanced_gold_eval
```

## Gate counts

Training decisions:

```text
{
  "('do_not_retrain_build_family_balanced_gold_eval', False)": 5,
  "('do_not_retrain_fix_route_first', False)": 2,
  "('do_not_retrain_monitor_or_split_subfamily', False)": 11
}
```

Positive candidate exactness:

```text
{
  "('fail_for_positive_label', 'Rhea')": 3,
  "('needs_manual_extraction', 'PubMed')": 43,
  "('partial_pass', 'Rhea')": 17
}
```

Negative evidence:

```text
{
  "('blocked', 'explicit_assay_no_conversion')": 13,
  "('conditional', 'hard_decoy_rule_generated_nontruth')": 13,
  "('forbidden', 'database_absence')": 1,
  "('forbidden', 'generator_miss')": 1
}
```

ChEBI Round18 structure cache:

```text
{
  "mapped": 39,
  "missing_structure": 6,
  "lookup_failed_or_missing_raw": 1
}
```

## Next production step

Round19 should not retrain. It should pick one P0 route family and make a dry-run fix:

- bile-acid deconjugation: explain or repair final clean filtering so validated multiple true products are not silently lost.
- hydroxycinnamate reduction/hydrolysis: wire a source-traceable rule path for chlorogenate -> trans-caffeate and rerun a target-only clean generation dry run.

Only after those route gates pass should we promote P1 external leads into exact-pair curation.

## Written artifacts

- `data/curation/round18_training_readiness_gate.csv`
- `data/curation/round18_p0_route_repair_worklist.csv`
- `data/curation/round18_positive_candidate_gate.csv`
- `data/curation/round18_negative_evidence_gate.csv`
- `data/curation/round18_training_sample_schema_contract.csv`
- `data/curation/round18_next_iteration_workplan.csv`
- `data/curation/round18_chebi_structure_status.csv`

## Round18 external evidence update

External-source contract rows:

```text
                     source_name                                                                                                                                             foodgut_training_contract
ECREACT / rxn-biocatalysis-tools                          Every imported positive needs reaction SMILES or substrate/product structures, EC/source provenance, split trace, and source database label.
                       EnzymeMap                                                    Exact names must be resolved to structures before a label is admitted; generic reaction text is rule context only.
                          gapseq                                   For gut-microbe deployment claims, reaction evidence should carry organism/sequence/phenotype or literature support when available.
                            Rhea Rhea exact reactions can seed candidate positives only after ChEBI structure mapping, local de-duplication, generator-route check, and leakage-safe split assignment.
```

Live Rhea/PubMed rerun gates:

```text
source_type      source_id                        candidate_gate  training_allowed_round18
       Rhea     RHEA:28130           exact_reaction_partial_pass                     False
       Rhea     RHEA:28314           exact_reaction_partial_pass                     False
       Rhea NO_RHEA_RESULT                          not_a_sample                     False
       Rhea     RHEA:35095           exact_reaction_partial_pass                     False
       Rhea     RHEA:17369           exact_reaction_partial_pass                     False
     PubMed  PMID:40863168 manual_exact_pair_extraction_required                     False
     PubMed  PMID:35273169 manual_exact_pair_extraction_required                     False
```

Additional written artifacts:

- `data/curation/round18_external_source_contract_update.csv`
- `data/curation/round18_live_rhea_pubmed_rerun_triage.csv`

