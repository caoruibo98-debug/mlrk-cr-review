# Generation 02 Review

## Frozen state

- Git tag: `generation_02`
- evaluated_commit: `72eab3d527e4dfba9dbf14b93a2bf65dda62dc0e`
- freeze_commit: `1ebf68cde316d9adf0e3642769b17cb72b375f41`
- Life-science application score: `4.32 / 5`
- Contract tests: `PASSED 9`

## CodeRabbit status

CodeRabbit CLI is installed in WSL (`0.5.3`). Authentication status returned exit code 0, but `coderabbit review --agent --base-commit b57c2ca` produced no NDJSON events, findings, or errors. This is not counted as a successful CodeRabbit review.

## Findings

1. The score now passes the requested `4/5` bar, but the evidence dimension mixes model-output evidence with benchmark-panel evidence. That is transparent in case rows, but the aggregate dimension should separate the two to avoid overstating model autonomy.
2. The evaluator now correctly detects structural matches for `(R)-naringenin` and CID-style hesperetin output. This is scientifically better than display-name matching and does not hide the remaining `quercitrin -> quercetin` generation miss.
3. Prediction payloads annotate each candidate with biochemical quality, but they do not yet provide a filtered "interpretation-ready" view. A web app could still display rejected lower-rank candidates too prominently.
4. The remaining true failure is candidate generation, not ranking: quercitrin does not generate quercetin in the top-10 output.
5. External head-to-head benchmark and wet-lab validation remain missing, so the status must stay `internal_mvp_only`.

## Required next iteration

1. Split model-output evidence from benchmark-panel traceability in the final score.
2. Add `interpretation_ready_top` and `quality_summary` to prediction payloads.
3. Add regression tests proving rejected candidates are excluded from the interpretation-ready view.
4. Freeze generation 03 and report whether the stricter score remains at or above `4/5`.
