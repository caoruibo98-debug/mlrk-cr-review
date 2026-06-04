# Round20 Reaction-Type Gap And Sample Strategy

## Working conclusion

Ray's diagnosis is mostly right, but the production blocker is more precise than 'sample type is not enough'. The current model fails because reaction-family coverage, source-backed exactness, generator route retention, negative-label semantics, and evaluation size are all weak. Reaction-type gaps are one major cause, but route loss can erase positives that already exist in `reactions.parquet`.

## Highest-risk local families

```text
module                                  reaction_category  reactions_pairs  clean_full_pos  clean_final_pos  gold_pairs  silver_pairs clean_final_retention_vs_reactions                                                                                                                     gap_labels_round20
     A Flavonoid, isoflavonoid, and polyphenol metabolism               82              13                0           0             0                                0.0 no_gold_or_silver_eval_anchor;complete_final_clean_loss;clean_full_to_final_route_loss;missing_ec_mapping;rule_only_not_training_ready
     A           isomerization_epimerization_racemization               78               0                0           0            88                                0.0                                                                                                              complete_final_clean_loss
     A                          missing_reaction_category               49               0                0           0            44                                0.0                                                              complete_final_clean_loss;missing_ec_mapping;rule_only_not_training_ready
     A                                     Heme synthesis               21               0                0           0             7                                0.0                                                                                                              complete_final_clean_loss
     A                                Caffeine metabolism               20               0                0           0             0                                0.0                                no_gold_or_silver_eval_anchor;complete_final_clean_loss;missing_ec_mapping;rule_only_not_training_ready
     B Flavonoid, isoflavonoid, and polyphenol metabolism               24              19                1           0             0                             0.0417   no_gold_or_silver_eval_anchor;severe_final_clean_loss;clean_full_to_final_route_loss;missing_ec_mapping;rule_only_not_training_ready
     B          group_transfer_reconjugation_modification              119               9                5           0           124                              0.042                                                                                 severe_final_clean_loss;clean_full_to_final_route_loss
     C Flavonoid, isoflavonoid, and polyphenol metabolism               47              13                0           0             0                                0.0 no_gold_or_silver_eval_anchor;complete_final_clean_loss;clean_full_to_final_route_loss;missing_ec_mapping;rule_only_not_training_ready
     C                           functional_group_removal               23               1                0           0            25                                0.0                                                                               complete_final_clean_loss;clean_full_to_final_route_loss
     C                 Isoquinoline alkaloid biosynthesis               22              10                0           0             0                                0.0 no_gold_or_silver_eval_anchor;complete_final_clean_loss;clean_full_to_final_route_loss;missing_ec_mapping;rule_only_not_training_ready
     C                          missing_reaction_category               22               0                0           0            22                                0.0                                                                                           complete_final_clean_loss;missing_ec_mapping
     C          group_transfer_reconjugation_modification              161              13                2           0           168                             0.0124                                                                                 severe_final_clean_loss;clean_full_to_final_route_loss
```

## External source roles

```text
                               source_name                                       source_role priority                                                                                                                         guardrail
                                      Rhea                     primary positive verification       P0 Rhea hit is a reaction truth source, not proof that a specific gut species performs it unless enzyme/microbe evidence is attached
ECReact / rxn4chemistry biocatalysis-model       broad enzymatic reaction coverage benchmark       P1                                                   not microbiome-specific; source database and structure mapping must be retained
                                 EnzymeMap mapped and corrected enzymatic reaction expansion       P1                                         BRENDA-derived entries still require exact participant and direction checks before import
                                RetroRules                  candidate generation rule source       P1                            template existence is not biological occurrence and must not be imported as a positive label by itself
       gapseq / ModelSEED / KEGG / MetaCyc       microbial pathway and genome-route evidence       P1                                               genome-model absence is condition-specific and not a universal biochemical negative
                                       VMH  gut microbiome reaction and metabolite relevance       P1                     VMH context supports relevance, but exact substrate/product identity still needs reaction-source verification
                              CAZy / CAZac               carbohydrate-active enzyme coverage       P1                                        enzyme family annotation needs substrate specificity before creating exact reaction labels
                                  gutMGene       microbe-metabolite-gene association context       P2                                                                 association evidence is not substrate-product conversion evidence
```

## Negative sample logic

```text
                        negative_type                                                        acceptable_use                         not_allowed_claim                                                             condition_scope
ranker_decoy_generated_same_substrate                        within-substrate ranking negative / hard decoy         biochemically impossible reaction unlabeled unless screened; model-training weight should reflect uncertainty
      organism_context_absent_pathway                             condition-specific microbe-route negative universal negative for all gut microbiota                                             organism/strain/genome-specific
           assay_tested_no_conversion                      strongest negative, but assay-condition-specific                  reaction can never occur                                                enzyme/strain/assay-specific
  direction_or_cofactor_invalid_decoy structural/process decoy after Rhea/ChEBI direction and balance check     false across all biochemical contexts                                                    route-condition-specific
        invalid_structure_or_unmapped                                      exclude from training/evaluation                           negative sample                                                      data-quality exclusion
```

## Pipeline

```text
 step_order                    stage                                                                                 pass_gate                          output_artifact
          1                  collect                       gap family selected with explicit route-loss/sample-evidence reason          roundXX_collection_worklist.csv
          2             verify_truth exact participant identity plus curated reaction/source reference; no analogy-only import roundXX_verified_positive_candidates.csv
          3     modify_manifest_only                          training_allowed is false until route/generator/split gates pass      roundXX_positive_candidate_gate.csv
          4           negative_logic             unknown is not treated as false; negative class is named and condition-scoped      roundXX_negative_candidate_gate.csv
          5             route_dryrun     target_generated=True and decoy_count>0 without disabling anti-cheat filters globally                 roundXX_route_dryrun.csv
          6 summarize_and_next_round                                    CSV plus review doc committed before next round starts          docs/reviews/roundXX_summary.md
```

## Written artifacts

- `data/curation/round20_reaction_type_gap_matrix.csv`
- `data/curation/round20_external_training_source_candidates.csv`
- `data/curation/round20_negative_sample_logic.csv`
- `data/curation/round20_sample_acquisition_pipeline.csv`
