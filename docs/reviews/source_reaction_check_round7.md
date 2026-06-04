# Round7 Source-Reaction Check

Date: 2026-06-04

## Working Conclusion

Ray's judgment is mostly correct: the current blocker is not only "sample count is small"; it is **reaction-family and evidence coverage under deployment-faithful generation**.

The model starts with a large raw positive pool, but production depends on whether a real candidate generator can enumerate the true product and whether the positive/negative labels are biologically defensible. Current evidence says many high-value food-gut reactions are lost because the specific reaction family or rule provenance is weak, not because the LTR ranker itself is the only weak part.

The production problem has four linked causes:

1. Reaction-type gaps: important gut/food transformations are not consistently covered by deployment rules.
2. Rule provenance gaps: a fired SMARTS rule is often not traceable to a source reaction, EC, Rhea, KEGG, MetaCyc, or exact paper evidence.
3. Label gaps: rule-generated products are candidates, not positive labels.
4. Negative-label gaps: "not found in database" is unlabeled, not a true negative.

## Files Created

- `data/curation/source_reaction_check_round7.csv`
- `data/curation/exact_evidence_extraction_queue_round7.csv`
- `data/curation/literature_pubmed_summary_round7.csv`
- `data/curation/rhea_source_reaction_summary_round7.csv`
- `data/curation/negative_sample_acquisition_plan_round7.csv`
- `data/curation/external_model_reaction_coverage_round7.csv`
- `data/curation/ncbi_pubmed_esummary_round7_raw.json`
- `data/curation/rhea_round7_raw.json`

## Round7 Results

`source_reaction_check_round7.csv` audited 11 recovered high-priority rule candidates from Round6.

Status counts:

- `microberx_rule_needs_source_recovery`: 7
- `source_reaction_candidate_verified_not_imported`: 2
- `deprecated_mapping_unresolved_in_current_reac_prop`: 1
- `source_reaction_verified_but_pair_not_exact_positive`: 1

Interpretation:

- 2 RetroRules candidates have usable source-reaction traces through MetaNetX/Rhea/KEGG/MetaCyc, but they still need direction, stereochemistry, license, exact pair extraction, and holdout checks before training.
- 7 MicrobeRX candidates still need source recovery or replacement with a rule that has EC/database/literature provenance.
- Quercetin -> taxifolin should **not** be used as a training positive from the current evidence. Rhea confirmed a related oxidoreductase chemistry, but not this exact pair.
- Caffeic acid -> dihydrocaffeic acid has promising literature identity, but the selected RetroRules legacy reaction maps through deprecated MetaNetX IDs that did not resolve cleanly in `reac_prop`; it needs manual MetaNetX/source-reaction review.

## External Model Comparison

The comparison sources point to a consistent standard.

- ECREACT / RXN biocatalysis model uses Rhea, BRENDA, PathBank, and MetaNetX and spans all 7 first-level EC classes. It treats EC/source information as part of the reaction representation, not as an afterthought.
- RetroRules / RetroPathRL treats reaction SMARTS templates as candidate-generation rules with source reaction IDs, EC fields, rule direction, and rule usage. A rule firing alone is not a positive label.
- gapseq uses curated reaction databases and genome/pathway evidence to reconstruct bacterial metabolic networks. This shows why organism-context and conditional negatives matter.
- GutBugDB and gutSMASH are useful for gut-microbiome biological plausibility and organism/enzyme context, but their predictions should not be imported as gold substrate-product positives without exact evidence.

Practical consequence: this project should not try to "invent" new reaction rules for training. It should only import:

- exact positive substrate-product pairs from curated databases or papers;
- rules with traceable source reaction provenance;
- weak context evidence as weak/context features, not gold labels;
- hard decoys only as ranking negatives, not biological negatives.

## Positive Sample Rule

Allowed future positive labels must satisfy all of:

1. substrate and product mapped to stable chemical identifiers;
2. exact substrate-product relationship appears in a database or paper;
3. direction is explicit or defensible;
4. organism, enzyme, strain, community, or assay context is captured when available;
5. evidence tier is recorded;
6. source license and database version are recorded;
7. holdout and split leakage guard is passed.

Not allowed as positive labels:

- rule-fired products with no exact source pair;
- related family reactions;
- mechanism-support-only papers;
- database absence reversals;
- predicted EC/gene context alone.

## Negative Sample Rule

Good negative data is harder than positive data.

Allowed:

- `hard_decoy`: same-substrate generated candidate not in curated positives; only for LTR contrast, not a biological impossibility claim.
- `assay_negative`: paper explicitly reports no conversion under stated conditions. This is the strongest negative.
- `conditional_negative`: named organism/model lacks pathway/gene/reaction support under a defined database version and context.

Not allowed:

- "PubMed/Rhea/MetaNetX did not find it" as a global negative.
- "Rule did not generate it" as a biological negative.
- "No known gut microbe" without versioned organism/genome search.

## Next Expansion Pipeline

Round8 should continue as:

1. Select the 2 RetroRules candidates with source-reaction traces.
2. Manually extract exact pair text/table evidence for pinoresinol -> lariciresinol and lariciresinol -> secoisolariciresinol.
3. Resolve caffeic acid -> dihydrocaffeic acid MetaNetX deprecated source reaction mapping.
4. For the 7 MicrobeRX-only candidates, search replacement rules in Rhea/MetaNetX/ECREACT/RetroRules or exact papers before importing.
5. Build `positive_sample_expansion_round8.csv` only from verified exact pairs.
6. Build `negative_assay_candidates_round8.csv` only when papers report explicit no-conversion assays.
7. Retrain only after the new positives and negatives pass split/holdout checks.

## Production Readiness Verdict

Current status: **not production-ready**.

Reason: the ranking model can score candidates, but the candidate generator and label base are not yet strong enough to support reliable production claims across food-gut reaction families.

Best defensible current claim:

> This is a deployment-faithful rule-generated candidate ranking prototype with evidence-aware curation manifests. It can prioritize chemically plausible and source-supported candidate metabolites, but it is not yet an externally validated predictor of gut microbial metabolism.

Minimum next improvement:

> Convert the Round7 evidence queue into verified exact positive rows and context-aware negative rows, then rerun generation recall and LTR evaluation on leakage-controlled splits.

