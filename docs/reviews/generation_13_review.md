# Generation 13 Review

## Production question

Can substrate and reaction-family coverage be reported like a production KPI, instead of only reporting one blended score?

## Changes

1. Added `tools/reaction_family_kpis.py` and `scripts/reaction_family_kpis.py`.
2. Added reproducible outputs:
   - `outputs/appraisal/reaction_family_kpis.json`
   - `outputs/appraisal/reaction_family_kpis.csv`
3. Added reaction-family KPI contract tests.
4. Integrated the KPI report into `outputs/appraisal/production_scorecard.json`.
5. Updated documentation and the production iteration ledger.

## KPI definitions

- `top5_hit_rate`: strict expected-product top-5 recall. Connectivity-only matches with full-InChIKey mismatches are not counted.
- `candidate_pool_recall`: fraction of cases where the expected product appears anywhere in the generated candidate pool.
- `benchmark_traceability_rate`: fraction of cases with any traceable benchmark or model evidence source.
- `model_evidence_rate`: fraction of cases where matched output evidence is not limited to benchmark-panel labels.

The key distinction is:

> If `candidate_pool_recall` is zero, ranker tuning cannot fix the family. The candidate generator must be improved first.

## Results

- Contract tests: `PASSED 33`
- Core strict top-5 hit rate: `6 / 6`
- Challenge strict top-5 hit rate: `8 / 22`
- Reaction families reported: `19`
- Challenge families reported: `18`
- Highest-priority gap: `candidate_generation_blocked`
- Readiness: `internal_mvp_only`

## What the KPI report shows

Working internal benchmark niches:

- `glycoside_hydrolysis`: `6 / 6` strict top-5 in the core panel.
- `hydroxycinnamate_reduction`: `3 / 3` strict top-5 in the challenge panel.
- `gallate_decarboxylation`, `isoflavone_deglucosylation`, `stilbene_reduction`, `stilbene_dehydroxylation`, and `urolithin_dehydroxylation`: current challenge cases are strict top-5 hits.

Production blockers:

- `ellagitannin_hydrolysis`, `ellagitannin_urolithin_multistep`, `flavanol_ring_fission`, `flavanone_ring_fission`, `flavonol_ring_fission`, `isoflavone_reduction`, `enterolignan-style lignan metabolism`, and `phenolic_ester_hydrolysis` are still candidate-generation-blocked in the current panel.
- All working families are still evidence-integration-blocked because their matched products are traceable to benchmark labels, not independent model-output enzyme, microbe, strain, or literature evidence.

## Audit interpretation

The current system remains a hybrid of rule-based candidate generation, chemistry-only learned-to-rank scoring, evidence overlays, and strict appraisal logic. Generation 13 improves production observability, not biological validation.

The family KPI report makes the next production step sharper:

1. Do not tune the ranker for families where the expected product is absent from the candidate pool.
2. Add reaction-specific candidate generation for blocked families first.
3. Add independent enzyme, microbe, strain, gene, or paper evidence for families already hitting top-5.
4. Keep external benchmark and wet-lab claims out of the product boundary until real outputs or experiments are imported.

## Next round

Generation 14 should harden API behavior around safe error states, request validation, and timeouts so the prediction surface can support a small internal web app without overstating biological certainty.
