# Evaluation Protocol

## Current internal metrics

The current mainline metric file is:

```text
outputs/modular/ltr/clean2_metrics.csv
```

Primary internal metric:

```text
per-substrate recall@5 for LTR_chem against random, ec_only, and Tanimoto
```

## Required next benchmark

Run the same compound set through:

- BioTransformer gut microbial or SuperBio mode.
- MicrobeRX.
- A GutBug-style EC/enzyme baseline if executable access is available.
- Tanimoto and EC-only internal baselines.

## Report these separately

1. Generation recall: whether the true product is generated.
2. Ranking recall: where the true product ranks when generated.
3. Evidence quality: whether enzyme, microbe, gene, or literature evidence is attached.
4. Failure reason: no rule, invalid structure, out-of-domain, or ranking miss.

## Minimum publication-safe claim

The system improves ranking of rule-generated candidates in the tested internal candidate universe. It does not prove biological occurrence or complete metabolite generation.

