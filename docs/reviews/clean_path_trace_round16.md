# Round16 Clean Candidate Path Trace

## Working conclusion

Ray's suspicion is directionally right, but it needs a sharper split:

1. The model does have reaction-family coverage gaps.
2. The immediate production blocker is even more concrete: source-backed positives already present in `reactions.parquet` can be lost before final `clean_candidates.parquet`.
3. The current high metrics are not production evidence because `clean2_metrics.csv` evaluates tiny panels (`n_pts=10`).

## Local path result

Loss-stage counts:

```text
{
  "target_lost_between_clean_full_and_clean": 2,
  "absent_from_clean_full_but_fullrule_hit": 1,
  "present_in_final_clean": 1
}
```

The four Round15 Rhea candidates are not new labels. Three of them are repair targets:

- bile-acid deconjugation positives reach `clean_candidates_full` but are absent from final clean candidates.
- chlorogenate hydrolysis is in `reactions.parquet` and full-rule diagnosis can hit it, but deployment-like clean generation misses it.
- daidzein glycoside hydrolysis already reaches final clean candidates and should only receive provenance enrichment.

## Rule sensitivity result

The rule-sampling table shows whether a target appears only under broad/unannotated rules or under EC-bearing rules. A target that only appears under EC=0 rules is not safe for production unless rule provenance is recovered.

```text
                      pair_key  require_ec  target_generated
BHTRKEVKTKCXOH__RUDATBOHQWOJDD       False              True
BHTRKEVKTKCXOH__RUDATBOHQWOJDD        True             False
CWVRJTMFETXNAD__QAIPRVGONGVQAS       False             False
CWVRJTMFETXNAD__QAIPRVGONGVQAS        True             False
KYQZWONCHDNPDP__ZQSIJRDFPHDXIC       False              True
KYQZWONCHDNPDP__ZQSIJRDFPHDXIC        True              True
RFDAIACWWDREDC__BHQCQFFYRZLCQQ       False              True
RFDAIACWWDREDC__BHQCQFFYRZLCQQ        True              True
```

## External retrieval result

Rhea triage counts:

```text
{
  "approved_biochemical_context": 14,
  "needs_literature_or_other_database": 6,
  "source_backed_existing_positive": 4
}
```

PubMed triage counts:

```text
{
  "literature_lead_for_manual_exact_pair_review": 15,
  "context_only": 13,
  "search_miss": 2,
  "possible_assay_negative_or_cohort_context": 1
}
```

Rhea is strong for exact biochemical reactions such as bile-acid deconjugation and chlorogenate hydrolysis. It is weak for several gut-microbiome transformation families where the evidence lives in primary papers or curated microbiome databases rather than Rhea exact reactions.

## Production implication

Do not retrain yet. First fix route recall and evidence gates:

1. Preserve all true products for substrates that have multiple validated products, or document the final clean filtering rule.
2. Promote only source-traceable rules for known full-rule hits that are absent from deployment-like clean generation.
3. Expand gold/test sets by reaction family, not by random database duplicates.
4. Treat no-hit generator outputs as hard decoys, not biological negatives.
5. Add assay negatives only when a paper/database explicitly reports no conversion under defined microbe/condition/time.

## Written artifacts

- `data/curation/round16_clean_path_trace.csv`
- `data/curation/round16_rule_sampling_sensitivity.csv`
- `data/curation/round16_rhea_query_triage.csv`
- `data/curation/round16_pubmed_query_triage.csv`
- `data/curation/round16_external_model_position.csv`
- `data/curation/round16_training_decisions.csv`
- `data/curation/round16_collect_verify_modify_pipeline.csv`
