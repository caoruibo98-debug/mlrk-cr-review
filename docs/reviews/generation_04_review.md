# Generation 04 Review

## Frozen state

- Target tag: `generation_04`
- Strict life-science application score: `4.43 / 5`
- Contract tests: `PASSED 12`
- Real food-glycoside panel: `6 / 6` expected products ranked in top 5.

## CodeRabbit status after user authentication

CodeRabbit authentication is now confirmed:

```text
authenticated=true
username=caoruibo98-debug
```

Two review attempts did not complete:

1. Full diff from `generation_01`:

```text
coderabbit review --agent --base cr_base_generation_01
```

Result: timed out after 10 minutes with no finding or error event.

2. Scoped review:

```text
coderabbit review --agent --base cr_base_generation_01 --dir mlrk_prod
```

Result: timed out after 10 minutes with no finding or error event.

Stored findings were checked and none were present. This generation therefore records CodeRabbit as authenticated but blocked by review timeout, not as a passed review.

## Model-quality changes

1. Added `mlrk_prod.glycoside_rescue.aromatic_o_glycoside_rescue_candidates`.
2. Added a conservative inference-time candidate rescue for aromatic O-glycosides.
3. The rescue is only applied to module `B` carbohydrate/glycoside inputs.
4. Rescue candidates are labeled with `candidate_source=curated_aromatic_o_glycoside_rescue`.
5. The rescue does not inject evidence and does not change model training.

## Evidence of improvement

- Before generation 04, `quercitrin -> quercetin` was a true generation miss.
- After generation 04, quercitrin ranks quercetin as top-1.
- The expected product is still counted as `benchmark_panel` evidence, not model-output evidence.
- Strict score improved from `4.15 / 5` to `4.43 / 5`.

## Remaining limitations

1. External BioTransformer/MicrobeRX/GutBug-style benchmark is still missing.
2. Wet-lab validation is still missing.
3. The model is strongest for food polyphenol glycoside aglycone prioritization, not for arbitrary gut microbiome metabolism.
4. CodeRabbit still needs a successful review run before it can be treated as an external code-quality gate.

