# Round20 BSH Rule Overlay Dry Run

## Working conclusion

A source-backed BSH EC repair overlay recovers both P0 bile-acid targets under the sampled require-EC generation path. This supports a targeted generator patch review, not retraining yet.

## Overlay manifest

```text
 source_id    rule_hash old_rule_ec_classes proposed_overlay_ec_class                        overlay_gate_status  overlay_allowed_round20
RHEA:16309 8f2848281b71                   0                         3         exact_bsh_rule_ec_repair_candidate                     True
RHEA:16309 6d2db3321cf5                   0                                        related_bsh_rule_context_only                    False
RHEA:19353 646b603ae957                   0                         3         exact_bsh_rule_ec_repair_candidate                     True
RHEA:19353 9d65a928f561                   0                                        related_bsh_rule_context_only                    False
RHEA:19353 cfbf068c044b                   0                                        related_bsh_rule_context_only                    False
RHEA:19353 6ccb1d8afabc                   3                                        related_bsh_rule_context_only                    False
RHEA:19353 ca068e9ce07d                   2                           blocked_non_bsh_or_wrong_direction_context                    False
```

## Dry-run result

```text
 source_id                        rule_variant  target_generated target_generated_ec_classes  decoy_count  final_clean_target_eligible_in_dryrun
RHEA:16309         baseline_sampled_require_ec             False                                        2                                  False
RHEA:16309 round20_exact_bsh_ec_repair_overlay              True                           3            3                                   True
RHEA:19353         baseline_sampled_require_ec             False                                        0                                  False
RHEA:19353 round20_exact_bsh_ec_repair_overlay              True                           3            1                                   True
```

## Decisions

```text
 source_id                                  decision                                                                                                         next_gate  training_allowed_round20
RHEA:16309 ready_for_targeted_generator_patch_review code review source-backed overlay integration; rerun clean candidate target-family test; then leakage split audit                     False
RHEA:19353 ready_for_targeted_generator_patch_review code review source-backed overlay integration; rerun clean candidate target-family test; then leakage split audit                     False
```

## Written artifacts

- `data/curation/round20_bsh_overlay_manifest.csv`
- `data/curation/round20_bsh_overlay_dryrun.csv`
- `data/curation/round20_bsh_clean_candidate_preview.csv`
- `data/curation/round20_bsh_overlay_decisions.csv`
