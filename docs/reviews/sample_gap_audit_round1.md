# Sample Gap Audit Round 1

## Scope

This audit answers whether the current production gap is mainly caused by missing reaction types, weak sample coverage, or model ranking quality.

Current authoritative training/evaluation files:

- Raw positive pool: `D:/CRB/FoodGut/positive_sample_db/data/processed/foodgut_modular_positive_candidate_pool_v2.csv`
- Reaction rule master: `D:/CRB/Food models/Dataset/reaction_rule_master.csv`
- LTR positives: `D:/CRB/Food models/scripts/ssrf/ml_ranking_kernel/outputs/modular/ltr/reactions.parquet`
- LTR full candidates: `D:/CRB/Food models/scripts/ssrf/ml_ranking_kernel/outputs/modular/ltr/clean_candidates_full.parquet`
- Rule-generation audit: `D:/CRB/Food models/scripts/ssrf/ml_ranking_kernel/outputs/gen_gap/fullrule_hit_miss.parquet`
- LTR metrics: `D:/CRB/Food models/scripts/ssrf/ml_ranking_kernel/outputs/modular/ltr/clean2_metrics.csv`
- Challenge holdouts: `D:/CRB/Food models/scripts/ssrf/ml_ranking_kernel_prod_candidate/data/reaction_family_gap_seed_manifest.tsv`
- Round 1 expansion queue: `D:/CRB/Food models/scripts/ssrf/ml_ranking_kernel_prod_candidate/data/curation/reaction_gap_candidate_sources_round1.csv`

## Current Evidence

Raw positive pool contains 18,163 rows and 12,171 unique substrate/product InChIKey pairs. The final `reactions.parquet` contains 10,550 unique positive reactions: 80 gold, 3,934 silver, and 6,536 weak. The full LTR candidate table has 50,011 rows, including 2,647 positives and 47,364 negatives.

The high metrics are not production-grade evidence. `clean2_metrics.csv` evaluates only `n_pts=10` per module. The D module reaches `r@5=1.000`, but that is a tiny internal scaffold/product-disjoint evaluation, not external validation. Separately, `fullrule_hit_miss.parquet` reports rule-generation recall, not ML accuracy; B gold is 10/10 but is dominated by glycoside/oligosaccharide cleavage.

## Main Finding

The user's suspicion is correct, with one refinement: the blocker is not simply "we have no samples." The stronger diagnosis is:

1. Some production-critical reaction families are missing or underrepresented in the final candidate-generation path.
2. Some evidence already exists in the raw positive pool but disappears before `clean_candidates_full.parquet`.
3. Multi-step gut microbial pathways are being judged through single-step candidate-generation and small internal ranking tests.
4. Negative samples are mostly generated candidate non-matches; they should be treated as unlabeled or hard-negative only when a source supports that interpretation.

## Holdout Pair Trace

Important examples from the current manifest:

- `chlorogenic_acid_caffeic_acid`: present in positive pool and `reactions.parquet`, absent from `clean_candidates_full`; rule-generation hit exists.
- `ellagic_acid_urolithin_a`: present in positive pool, absent from `reactions.parquet` and candidates.
- `daidzein_equol`: present in positive pool, absent from `reactions.parquet` and candidates.
- `enterodiol_enterolactone`: present in positive pool and `reactions.parquet`, absent from candidates.
- `daidzein_dihydrodaidzein` and `genistein_dihydrogenistein`: local evidence exists, but clean candidates do not retain the positives.
- `puerarin_daidzein`, `punicalagin_ellagic_acid`, and direct catechin/epicatechin valerolactone holdouts are absent from the current positive pool.

## External Model Lessons

BioTransformer 3.0 combines rule-based and machine-learning metabolism prediction and includes human gut microbial transformations. MicrobeRX uses reaction rules plus evidence tables with reaction IDs, EC annotations, origin, PubMed, and microbial reaction matrices. GutBug/GutBugDB emphasize EC prediction, reaction centers, enzymes, strains, and population abundance context.

Therefore, the next production-oriented sample schema must include more than substrate/product:

- reaction family and transformation class
- single-step versus multistep path stage
- enzyme, EC, gene, or reaction database ID
- microbe/strain or community context
- source PMID/DOI/PMCID and evidence type
- representation status: SMILES, InChIKey, stereochemistry, fragment/polymer note
- split role and leakage policy

## Positive Sample Policy

Only promote a row to training-positive when it has a traceable substrate/product pair and at least one of:

- direct enzyme assay
- strain fermentation or human fecal fermentation with named product
- curated database reaction with reaction ID and EC/microbe evidence
- paper-supported pathway step with extractable product identity

Endpoint-only multi-step claims should stay as holdout or pathway evidence until intermediate steps are curated.

## Negative Sample Policy

Do not treat all unobserved reactions as true negatives. Prefer:

- class-contrast negatives: O-glycoside hydrolysis versus C-glycoside cleavage
- context negatives: non-producer strain/community evidence
- hard negatives: generated plausible products that are close in structure but not supported by a reviewed source
- unlabeled pool: all other generated non-matches

## Round 2 Work

1. For each row in `reaction_gap_candidate_sources_round1.csv`, extract exact independent substrate/product rows from the cited source.
2. Verify SMILES/InChIKey with PubChem/ChEBI/HMDB where possible.
3. Assign `training_allowed=true` only if the row is not a benchmark/challenge holdout and passes source verification.
4. Add a separate manifest for negatives with `negative_type` and `negative_evidence_basis`.
5. Re-run candidate-generation recall before retraining the ranker.
