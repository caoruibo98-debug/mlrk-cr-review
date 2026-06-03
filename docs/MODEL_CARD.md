# Model Card

Generated from `outputs/appraisal/production_scorecard.json`.

## Status

- Readiness status: `internal_mvp_only`
- Generalization evidence level: `3` (internal benchmark evidence only; no external tool outputs or wet-lab validation imported)
- External result status: `awaiting_external_outputs`
- External AI/code review status: `timed_out` (claim allowed: `false`)
- Repository doctor status: `ready`
- Fresh-clone release-candidate status: `passed`

## System Type

This is a hybrid internal research system:

- rule-based candidate generation,
- chemistry-only learned-to-rank scoring,
- evidence overlay,
- contract-checked internal API.

It is not an externally validated biological occurrence predictor.

## Biochemical Quality Boundary

Candidate quality is a structural plausibility screen. High tier requires traceable enzyme, microbe, EC, PMID, or known-product evidence; structural-only candidates are capped at medium.

## Intended Use

Internal research MVP for food-polyphenol candidate generation/ranking, strongest on glycoside aglycone release.

## Not Intended Use

- Broad prediction of all food-derived gut microbial metabolites.
- Strain-aware or genome-aware gut microbiome metabolism prediction.
- Wet-lab occurrence probability.
- Consumer health, clinical, or treatment recommendation.

## Evaluation Snapshot

| Evaluation | Cases | Strict top-5 hits | Strict top-5 rate | Score |
| --- | ---: | ---: | ---: | ---: |
| Core food-glycoside panel | 6 | 6 | 1.000 | 3.88 |
| Production challenge panel | 22 | 8 | 0.364 | 2.81 |

Internal LTR_chem mean recall@5 is `0.94` versus random `0.532`, EC-only `0.525`, and Tanimoto `0.716`.

Reaction-family KPI count: `19`.

## Readiness Issues

- `external_benchmark_results_missing`: External benchmark inputs and import harness are ready, but real BioTransformer/MicrobeRX/GutBug-style outputs are not scored yet.
- `wet_lab_missing`: Wet-lab validation is not complete; public biological claims must stay limited.
- `public_web_not_allowed`: Only internal research MVP deployment is allowed by this manifest.

## Remaining Gaps

- Challenge-panel failures are dominated by candidate-generation recall, not ranking.
- Reaction-family coverage is still narrow for reductions, dehydroxylations, ring fission, decarboxylation, and multi-step gut microbial pathways.
- Reaction-family KPI reporting is now available, but KPI strength is internal benchmark evidence only.
- External tools have not yet been run on the exact same substrate panel.
- No gene/genome or strain-level abundance context is connected to predictions.
- Candidate generation still constrains the ceiling; ranking cannot recover products absent from the candidate set.
- Model-output evidence coverage remains incomplete even when benchmark traceability is present.

## Claim Boundary

Scores are ranking and prioritization signals over generated candidates. They are not wet-lab probabilities, occurrence probabilities, clinical signals, or consumer health recommendations.
