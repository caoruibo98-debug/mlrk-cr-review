# Sample Route Trace Round 3

## Purpose

Round 3 traces Round 2 positive candidates through the route:

`positive_pool -> reactions.parquet -> ltr_clean_candidates.py -> clean_candidates_full.parquet`

The goal is to separate true raw-sample gaps from candidate-generation and filtering gaps.

## Artifacts

- Trace script: `tools/trace_sample_route_round3.py`
- Trace CSV: `data/curation/sample_route_trace_round3.csv`
- Input manifest: `data/curation/reaction_gap_positive_manifest_round2.csv`

Run:

```powershell
python tools\trace_sample_route_round3.py
```

## Current Result

Round 3 traced 75 Round 2 rows across 12 reaction families.

| Route loss stage | Rows | Meaning |
| --- | ---: | --- |
| `holdout_leakage_guard` | 22 | Exact benchmark/holdout pairs; do not train. |
| `selected_substrate_missing_from_clean_output` | 13 | Substrate selected, source evidence exists, but no clean output row for that substrate. |
| `generation_miss_drops_substrate` | 9 | Pair is in fullrule audit but target product is not generated; substrate is dropped. |
| `no_true_product_generated_for_substrate` | 8 | No true product was generated for that substrate. |
| `already_reaches_clean_candidates` | 8 | Already present in clean candidates; use as control. |
| `positive_pool_to_reactions_filter` | 8 | Upstream evidence exists but `ltr_build.py` did not admit it into `reactions.parquet`. |
| `target_product_not_generated` | 2 | Same substrate has clean candidates, but this specific target product is missing. |
| `generated_in_audit_but_missing_from_clean` | 2 | Fullrule audit says generated, but clean candidates do not keep the target product. |
| `not_covered_by_fullrule_audit` | 2 | Needs targeted generation audit. |
| `external_curation_needed` | 1 | No local high-evidence positive-pool row after strict Round 2 filtering. |

## Code Gates

The main route gates are in `modular/ltr_clean_candidates.py`:

- line 59: `drop_duplicates("sb")` switches from reaction-level to substrate-level candidate generation.
- line 66: another substrate-level deduplication happens after priority/weak sampling.
- line 85: `if not gen_truth: continue` drops the substrate if no true product is generated.
- lines 110-111: the final table keeps only substrates with at least one positive and at least one decoy.

The upstream gates in `modular/ltr_build.py` include:

- line 40: `is_single_step` filter.
- line 42: `structure_feasibility_flag` must be parseable.
- line 90: reaction-level `(sb, pb)` deduplication and tier max.

## Interpretation

The user's core judgment is supported: important reaction types and real gut-microbial reaction families are not reaching production candidates. However, Round 3 shows that the next action is not simply "collect more samples." The first production blocker is route recovery:

1. Existing high-evidence reactions are lost before or during clean candidate generation.
2. The candidate generator drops entire substrates when no true product is generated.
3. Multi-step reactions such as urolithin, equol, enterolignan, and flavanol ring-fission need stage-aware candidate generation rather than endpoint shortcuts.
4. Only one Round 2 family, `isoflavone_c_glycoside_conversion`, currently needs external curation before local route recovery can begin.

## Next Priority

1. For `selected_substrate_missing_from_clean_output`, run targeted rule-generation traces and compare rule sets used by `fullrule_hit_miss.parquet` versus `clean_candidates_full.parquet`.
2. For `generation_miss_drops_substrate` and `target_product_not_generated`, add or import source-backed reaction rules before retraining.
3. For `positive_pool_to_reactions_filter`, inspect whether `is_single_step`, structure parsing, module mapping, or canonicalization removed valid reactions.
4. Keep `holdout_leakage_guard` rows out of training.
5. Curate independent C-glycoside examples from Dorea/puerarin literature before training that family.
