# Generation 03 Review

## Frozen state

- Git tag: `generation_03`
- Commit: `7c939c7`
- Strict life-science application score: `4.15 / 5`
- Contract tests: `PASSED 10`
- Smoke test: `rutin` top-1 remains `quercetin`; `interpretation_ready_top` excludes rejected rank 4 when top-5 is requested.

## CodeRabbit status

CodeRabbit CLI is installed in WSL (`0.5.3`). A review was attempted with:

```text
coderabbit review --agent --base-commit 1ebf68c
```

The CLI emitted NDJSON auth events but did not review code:

- `awaiting_browser_auth`
- `automatic_login_failed`
- `authentication_failed`

Reason: browser authentication timed out. A user-side CodeRabbit login is required before CodeRabbit can provide issues.

## Findings

1. The stricter score remains above the requested bar: `4.15 / 5`.
2. The score is more honest than generation 02 because model-output evidence (`0.667`) is separated from benchmark traceability (`0.833`).
3. The model is now usable as an internal food-glycoside gut-metabolism prioritizer. It is not yet a general gut microbiome metabolism engine.
4. The remaining true panel miss is `quercitrin -> quercetin`; this is a candidate-generation/rule-coverage issue, not an evaluation naming issue.
5. Public production remains blocked by missing external BioTransformer/MicrobeRX/GutBug-style benchmark and missing wet-lab validation.

## Next work if continuing beyond 4/5

1. Fix the quercitrin rhamnoside hydrolysis rule path and prove it does not introduce new false positives.
2. Add a small external head-to-head runner for BioTransformer/MicrobeRX/GutBug-style baselines.
3. Add sample-specific input: taxa, EC, gene, or strain abundance should rerank or annotate candidates.
4. Complete CodeRabbit authentication and rerun review on `generation_01..generation_03`.

