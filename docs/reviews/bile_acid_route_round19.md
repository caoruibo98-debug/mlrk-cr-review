# Round19 Bile Acid Route Audit

## Working conclusion

The P0 bile-acid issue is not a request for more positive rows. Both target reactions are already source-backed and present before the final clean path. The final route drops them for two related reasons: one target has no EC-bearing generating rule, and the other has EC-bearing rules in the full pool but they are not deterministically included in the sampled require-EC generator.

## Route loss

```text
 source_id                       pair_key clean_full_target_ec_class  clean_pair_rows clean_ec_classes_for_substrate                                                target_loss_cause_round19
RHEA:16309 BHTRKEVKTKCXOH__RUDATBOHQWOJDD                          0                0                              1 target_reaches_clean_full_as_ec0_then_absent_from_ec_bearing_final_clean
RHEA:19353 RFDAIACWWDREDC__BHQCQFFYRZLCQQ                          0                0                              2 target_reaches_clean_full_as_ec0_then_absent_from_ec_bearing_final_clean
```

## Rule variant dry run

```text
 source_id             rule_variant  rule_count  target_generated target_generated_ec_classes
RHEA:16309 sampled_800_no_ec_filter         800              True                           0
RHEA:16309   sampled_800_require_ec         563             False                            
RHEA:16309    all_pool_no_ec_filter        2359              True                           0
RHEA:16309      all_pool_require_ec        1702             False                            
RHEA:19353 sampled_800_no_ec_filter         800              True                           0
RHEA:19353   sampled_800_require_ec         563             False                            
RHEA:19353    all_pool_no_ec_filter        2359              True                       0;2;3
RHEA:19353      all_pool_require_ec        1702              True                         2;3
```

## Target-generating rule summary

```text
 source_id    rule_hash  rule_ecc                     source_datasets  pool_rows_for_rule  ec_repair_needed
RHEA:16309 8f2848281b71         0 microberx_reaction_rules_integrated                   3              True
RHEA:16309 6d2db3321cf5         0 microberx_reaction_rules_integrated                   2              True
RHEA:19353 646b603ae957         0 microberx_reaction_rules_integrated                   2              True
RHEA:19353 9d65a928f561         0 microberx_reaction_rules_integrated                   2              True
RHEA:19353 cfbf068c044b         0 microberx_reaction_rules_integrated                   1              True
RHEA:19353 6ccb1d8afabc         3 microberx_reaction_rules_integrated                   1             False
RHEA:19353 ca068e9ce07d         2 microberx_reaction_rules_integrated                   1             False
```

## External evidence

- RHEA:16309: taurochenodeoxycholate + H2O = chenodeoxycholate + taurine; Rhea page reports EC 3.5.1.74.
- RHEA:19353: glycocholate + H2O = cholate + glycine; Rhea/ENZYME reports EC 3.5.1.24.
- PubMed BSH leads are context/mechanism evidence until exact assay tables are manually extracted.

Rhea exact rows:

```text
 source_id rhea_status                                              rhea_equation  balanced expected_ec
RHEA:16309    approved taurochenodeoxycholate + H2O = chenodeoxycholate + taurine      True    3.5.1.74
RHEA:19353    approved                     glycocholate + H2O = cholate + glycine      True    3.5.1.24
```

PubMed leads:

```text
    pmid                                                                                                        title                                           training_use
38617281         Chemoproteomic profiling of substrate specificity in gut microbiota-associated bile salt hydrolases. mechanism_or_context_only_until_exact_assay_extraction
18757757 Functional and comparative metagenomic analysis of bile salt hydrolase activity in the human gut microbiome. mechanism_or_context_only_until_exact_assay_extraction
33526676     Lactobacillus bile salt hydrolase substrate specificity governs bacterial fitness and host colonization. mechanism_or_context_only_until_exact_assay_extraction
```

## Round20 gate

```text
 source_id                            production_blocker_round19                                                                        fix_action                                                                                                          success_criterion  training_allowed_round19
RHEA:16309 ec_provenance_gap_no_ec_bearing_rule_generates_target add_source_backed_ec_to_existing_bsh_rule_or_overlay_then_rerun_require_ec_dryrun target_generated=True under sampled require-EC or targeted source-backed overlay, then pair appears in final clean dry run                     False
RHEA:19353                 sampled_require_ec_rule_selection_gap         deterministically_include_source_backed_bsh_rules_in_require_ec_generator target_generated=True under sampled require-EC or targeted source-backed overlay, then pair appears in final clean dry run                     False
```

## Written artifacts

- `data/curation/round19_bile_acid_route_loss.csv`
- `data/curation/round19_bile_acid_rule_variant_dryrun.csv`
- `data/curation/round19_bile_acid_target_generating_rules.csv`
- `data/curation/round19_bile_acid_external_evidence.csv`
- `data/curation/round19_bile_acid_pubmed_leads.csv`
- `data/curation/round19_bile_acid_fix_recommendations.csv`
