# Release Summary

## Final Position

- Production position: `internal_research_mvp`
- Readiness status: `internal_mvp_only`
- Public or clinical ready: `false`

The release candidate is suitable for internal research prioritization of rule-generated food-polyphenol metabolite candidates. It is not ready for public consumer, clinical, wet-lab probability, or broad gut microbiome metabolism claims.

## Internal And Ablation Metrics

- LTR_chem mean recall@5: `0.94`
- Random mean recall@5: `0.532`
- EC-only mean recall@5: `0.525`
- Tanimoto mean recall@5: `0.716`
- LTR_chem minus random recall@5: `0.408`
- LTR_chem minus EC-only recall@5: `0.415`
- LTR_chem minus Tanimoto recall@5: `0.224`
- Anti-cheat pass: `True`

## Panel Results

| Panel | Cases | Strict top-5 hits | Strict top-5 rate | Score |
| --- | ---: | ---: | ---: | ---: |
| Core | 6 | 6 | 1.0 | 3.88 |
| Challenge | 22 | 8 | 0.364 | 2.81 |

Challenge failure types: `{'expected_product_not_generated': 14, 'hit_top5_benchmark_only': 8}`.

## Cross-Family Coverage

- Reaction-family KPI status: `ready`
- Reaction families reported: `19`

## External Status

- External benchmark result status: `awaiting_external_outputs`
- External benchmark claim: external inputs and import harness are ready; real external outputs are not imported
- External AI/code review tool: `CodeRabbit`
- External AI/code review status: `timed_out`
- External AI/code review claim allowed: `false`

## Engineering Reproducibility

- Fresh-clone status: `passed`
- Fresh-clone commit: `fda7f6c8eff8ccd99eb3a6b7275037c2cfcacd2c`
- Repository doctor status: `ready`

## Remaining Blockers

- Challenge-panel failures are dominated by candidate-generation recall, not ranking.
- Reaction-family coverage is still narrow for reductions, dehydroxylations, ring fission, decarboxylation, and multi-step gut microbial pathways.
- Reaction-family KPI reporting is now available, but KPI strength is internal benchmark evidence only.
- External tools have not yet been run on the exact same substrate panel.
- No gene/genome or strain-level abundance context is connected to predictions.
- Candidate generation still constrains the ceiling; ranking cannot recover products absent from the candidate set.
- Model-output evidence coverage remains incomplete even when benchmark traceability is present.

## Readiness Issues

- `external_benchmark_results_missing`: External benchmark inputs and import harness are ready, but real BioTransformer/MicrobeRX/GutBug-style outputs are not scored yet.
- `wet_lab_missing`: Wet-lab validation is not complete; public biological claims must stay limited.
- `public_web_not_allowed`: Only internal research MVP deployment is allowed by this manifest.

## Final Claim Boundary

Use this release as an internal research MVP. Do not present it as an externally validated predictor, wet-lab occurrence model, strain-aware gut microbiome metabolism engine, clinical tool, or consumer health recommendation product.
