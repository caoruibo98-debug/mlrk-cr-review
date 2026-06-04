# Round22 P0 Route Dry Run

## Working conclusion

Round22 tests whether Round21 source-backed positives are generator-route failures. Training remains blocked for every candidate. A candidate can only advance after exact evidence, target generation, and external decoy screening all pass.

## Decision counts

```text
                            decision_round22  count
   rule_family_missing_for_current_generator     12
do_not_route_repair_until_evidence_recurated      2
```

## Candidate decisions

```text
candidate_id      substrate_name                      product_name  any_target_generated  require_ec_target_generated  max_decoy_count                             decision_round22  training_allowed_round22
  R21-P0-001         urolithin C                       urolithin A                 False                        False                0    rule_family_missing_for_current_generator                     False
  R21-P0-002            dopamine                        m-tyramine                 False                        False               21    rule_family_missing_for_current_generator                     False
  R21-P0-003         urolithin C                    isourolithin A                 False                        False                0    rule_family_missing_for_current_generator                     False
  R21-P0-004         urolithin A                       urolithin B                 False                        False                0    rule_family_missing_for_current_generator                     False
  R21-P0-005      isourolithin A                       urolithin B                 False                        False                0    rule_family_missing_for_current_generator                     False
  R21-P0-006        urolithin M5                      urolithin M6                 False                        False                0    rule_family_missing_for_current_generator                     False
  R21-P0-007        ellagic acid                      urolithin M5                 False                        False                0    rule_family_missing_for_current_generator                     False
  R21-P0-008        urolithin M6                       urolithin C                 False                        False                0    rule_family_missing_for_current_generator                     False
  R21-P0-009 protocatechuic acid                          catechol                 False                        False               26 do_not_route_repair_until_evidence_recurated                     False
  R21-P0-010        urolithin M6                      urolithin M7                 False                        False                0    rule_family_missing_for_current_generator                     False
  R21-P0-011        urolithin M7                       urolithin A                 False                        False                0    rule_family_missing_for_current_generator                     False
  R21-P0-012      corticosterone        11beta-hydroxyprogesterone                 False                        False               25 do_not_route_repair_until_evidence_recurated                     False
  R21-P0-013   hydrocaffeic acid 3-(3-hydroxyphenyl)propionic acid                 False                        False               40    rule_family_missing_for_current_generator                     False
  R21-P0-014         urolithin D                       urolithin C                 False                        False                0    rule_family_missing_for_current_generator                     False
```

## Variant target generation

```text
          rule_variant  targets_generated  final_clean_eligible  candidates
    all_pool_all_rules                  0                     0          14
   all_pool_require_ec                  0                     0          14
       sampled_all_800                  0                     0          14
sampled_require_ec_800                  0                     0          14
```

## Written artifacts

- `data/curation/round22_p0_route_dryrun.csv`
- `data/curation/round22_target_generating_rules.csv`
- `data/curation/round22_decoy_screening_manifest.csv`
- `data/curation/round22_p0_route_decisions.csv`
