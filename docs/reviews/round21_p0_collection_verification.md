# Round21 P0 Collection And Verification

## Working conclusion

Round21 confirms that the next audit should not simply add more rows to training. The P0 batch contains source-backed positives that are already in `reactions.parquet` but absent from final clean candidates. The right next move is exact evidence extraction plus route dry-run, then negative screening.

## Candidate worklist summary

```text
candidate_id      substrate_name                      product_name tier                       verification_status_round21                               route_gap_round21  training_allowed_round21
  R21-P0-001         urolithin C                       urolithin A gold  literature_verified_rhea_absent_not_training_yet reactions_to_clean_full_loss;final_clean_absent                     False
  R21-P0-002            dopamine                        m-tyramine gold      verified_exact_reaction_rhea_plus_literature reactions_to_clean_full_loss;final_clean_absent                     False
  R21-P0-003         urolithin C                    isourolithin A gold  literature_verified_rhea_absent_not_training_yet reactions_to_clean_full_loss;final_clean_absent                     False
  R21-P0-004         urolithin A                       urolithin B gold  literature_verified_rhea_absent_not_training_yet reactions_to_clean_full_loss;final_clean_absent                     False
  R21-P0-005      isourolithin A                       urolithin B gold  literature_verified_rhea_absent_not_training_yet reactions_to_clean_full_loss;final_clean_absent                     False
  R21-P0-006        urolithin M5                      urolithin M6 gold  literature_verified_rhea_absent_not_training_yet reactions_to_clean_full_loss;final_clean_absent                     False
  R21-P0-007        ellagic acid                      urolithin M5 gold  literature_verified_rhea_absent_not_training_yet reactions_to_clean_full_loss;final_clean_absent                     False
  R21-P0-008        urolithin M6                       urolithin C gold  literature_verified_rhea_absent_not_training_yet reactions_to_clean_full_loss;final_clean_absent                     False
  R21-P0-009 protocatechuic acid                          catechol gold blocked_pubmed_title_mismatch_requires_recuration reactions_to_clean_full_loss;final_clean_absent                     False
  R21-P0-010        urolithin M6                      urolithin M7 gold  literature_verified_rhea_absent_not_training_yet reactions_to_clean_full_loss;final_clean_absent                     False
  R21-P0-011        urolithin M7                       urolithin A gold  literature_verified_rhea_absent_not_training_yet reactions_to_clean_full_loss;final_clean_absent                     False
  R21-P0-012      corticosterone        11beta-hydroxyprogesterone gold             blocked_pending_external_verification reactions_to_clean_full_loss;final_clean_absent                     False
  R21-P0-013   hydrocaffeic acid 3-(3-hydroxyphenyl)propionic acid gold        literature_verified_rhea_partial_or_absent reactions_to_clean_full_loss;final_clean_absent                     False
  R21-P0-014         urolithin D                       urolithin C gold  literature_verified_rhea_absent_not_training_yet reactions_to_clean_full_loss;final_clean_absent                     False
```

## Source queries

```text
              query_id source                                                                                   result_summary                                                                    next_action
  R21-Q-RHEA-UROLITHIN   Rhea 0 returned records; use PubMed/literature as primary evidence for gut urolithin dehydroxylation. manual exact participant mapping against PubMed-curated urolithin paper tables
   R21-Q-RHEA-DOPAMINE   Rhea                                   Approved Rhea reaction: dopamine + AH2 = 3-tyramine + A + H2O.              use as exact positive verification and generator target test seed
   R21-Q-RHEA-CAFFEATE   Rhea                 0 returned records; hydrocaffeic/caffeic dehydroxylation remains literature-led.         use PubMed PMID 32067637 and exact structure check before route repair
R21-Q-PUBMED-UROLITHIN PubMed                                                     Returned PMIDs 41298472, 39856097, 37494568.                     extract exact substrate/product rows and sentence evidence
 R21-Q-PUBMED-CATECHOL PubMed                                                                          Returned PMID 32067637. extract substrate panel and exact products for catechol dehydroxylation family
```

## Family route actions

```text
 round          reaction_family  candidate_count  verified_or_literature_supported_count                                 route_gap_types                                   next_generator_gate                                                                     next_negative_gate  training_allowed_round21
    21 catechol_dehydroxylation                1                                       1 final_clean_absent;reactions_to_clean_full_loss exact target route dry-run before any training import generated decoys must be screened against Rhea/ECReact/EnzymeMap/current positive pool                     False
    21          decarboxylation                1                                       0 final_clean_absent;reactions_to_clean_full_loss exact target route dry-run before any training import generated decoys must be screened against Rhea/ECReact/EnzymeMap/current positive pool                     False
    21 functional_group_removal                5                                       4 final_clean_absent;reactions_to_clean_full_loss exact target route dry-run before any training import generated decoys must be screened against Rhea/ECReact/EnzymeMap/current positive pool                     False
    21     polyphenol_urolithin                7                                       7 final_clean_absent;reactions_to_clean_full_loss exact target route dry-run before any training import generated decoys must be screened against Rhea/ECReact/EnzymeMap/current positive pool                     False
```

## Written artifacts

- `data/curation/round21_p0_collection_worklist.csv`
- `data/curation/round21_source_query_manifest.csv`
- `data/curation/round21_route_repair_actions.csv`
