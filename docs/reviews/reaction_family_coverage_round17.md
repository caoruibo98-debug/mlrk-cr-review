# Round17 Reaction Family Coverage Audit

## Working conclusion

The production blocker is not just total sample count. It is a combined family-coverage and route-recall problem:

1. Some source-backed positives already exist locally but are not visible to final clean candidates.
2. Several gut-microbiome reaction families have external evidence but weak final clean/gold representation.
3. Current metrics are still too small for production interpretation (`clean2_metrics.csv` n_pts=[10]).

## Highest-priority families

```text
                                   family priority_level  priority_score                  dominant_gap_type  clean_pairs  clean_gold_pairs                                                                                        next_action
bile_acid_deconjugation_and_lipid_context             P0               8  route_repair_before_sample_import           96                 4                          repair clean route for existing positives; do not import duplicate labels
hydroxycinnamate_reduction_and_hydrolysis             P0               8  route_repair_before_sample_import            2                 0                          repair clean route for existing positives; do not import duplicate labels
      bile_acid_secondary_transformations             P1               4           gold_test_visibility_gap           14                 0   promote already curated gold/silver positives into family-balanced evaluation where leakage-safe
        glucuronide_sulfate_deconjugation             P1               4           gold_test_visibility_gap           62                 0   promote already curated gold/silver positives into family-balanced evaluation where leakage-safe
                  polyphenol_ring_fission             P1               4           gold_test_visibility_gap            0                 0   promote already curated gold/silver positives into family-balanced evaluation where leakage-safe
             tryptophan_indole_metabolism             P1               4           gold_test_visibility_gap           53                 3   promote already curated gold/silver positives into family-balanced evaluation where leakage-safe
         tma_choline_carnitine_metabolism             P1               4           gold_test_visibility_gap          122                 4   promote already curated gold/silver positives into family-balanced evaluation where leakage-safe
   azoreductase_nitroreductase_xenobiotic             P2               3 source_backed_family_undercoverage            0                 0 extract exact Rhea/PubMed substrate-product pairs, map structures, de-dup, then run generator gate
```

## External retrieval summary

Rhea Round17:

```text
{
  "approved_exact_or_specific_reaction_lead": 18,
  "generic_rule_template_context": 3,
  "no_direct_rhea_hit": 5
}
```

PubMed Round17:

```text
{
  "context_only": 21,
  "enzyme_or_mechanism_lead": 5,
  "literature_lead_for_manual_exact_pair_review": 41,
  "search_miss": 2
}
```

Rhea gives exact or generic reaction support for glucuronide/sulfate deconjugation, tryptophan/indole transformations, and choline/carnitine to TMA. PubMed gives many leads for tryptophan/indole and TMA biology, but most are association/pathway contexts and must not become labels without exact reaction extraction.

## Positive-label gate

Round17 keeps all new external leads as `training_allowed_round17=False`. A lead can move forward only after:

1. exact substrate and product are extracted,
2. ChEBI/PubChem/InChIKey mapping passes,
3. the pair is de-duplicated against the raw pool and `reactions.parquet`,
4. the clean generator can route the candidate or a source-traceable rule overlay is justified,
5. the label is assigned to a family-balanced split without product leakage.

## Negative-label gate

No-hit is not a biological negative. Use rule-generated non-truth products as ranking decoys only. Use assay negatives only when a paper or database explicitly reports no conversion under defined conditions.

## Written artifacts

- `data/curation/round17_reaction_family_coverage_matrix.csv`
- `data/curation/round17_rhea_query_triage.csv`
- `data/curation/round17_pubmed_query_triage.csv`
- `data/curation/round17_positive_expansion_queue.csv`
- `data/curation/round17_negative_sampling_logic.csv`
- `data/curation/round17_github_external_coverage_contract.csv`
- `data/curation/round17_training_decisions.csv`
