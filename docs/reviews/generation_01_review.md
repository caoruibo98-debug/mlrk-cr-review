# Generation 01 Review

## Frozen state

- Git tag: `generation_01`
- Commit: `b57c2ca`
- Life-science application score: `3.65 / 5`
- Contract tests: `PASSED 8`

## CodeRabbit status

CodeRabbit CLI review was attempted but did not run.

- `coderabbit` was not present on Windows.
- WSL is available.
- The official installer started but failed because WSL lacks `unzip`.
- Installing `unzip` through `sudo apt-get` is blocked because WSL requires a sudo password.

This generation therefore uses local AI review only. It must not be reported as a CodeRabbit result.

## Findings

1. The gold-panel evaluator matches expected products by display name only. This creates false misses for stereochemical names and PubChem/CID-style names, for example `(R)-naringenin` versus `naringenin`.
2. The evaluator does not resolve expected products to canonical structures or InChIKey blocks, so it cannot distinguish a true generation miss from a naming failure.
3. Evidence coverage is undercounted when the model produces the expected structure but the internal evidence lookup misses the substrate-product pair.
4. The new biochemical sanity gate correctly catches phosphorus introduction, tiny byproducts, large carbon gains, and missing evidence, but it is not yet used to produce a filtered "safe for interpretation" view.
5. Readiness remains `internal_mvp_only` because external head-to-head benchmark and wet-lab validation are still missing.

## Required next iteration

1. Resolve panel expected products to SMILES/InChIKey and match by structure before falling back to names.
2. Add a curated real-panel evidence overlay for known benchmark reactions, clearly labeled as benchmark evidence and not used for ranking.
3. Add a quality-filtered candidate view so UI/API users do not accidentally interpret rejected candidates as credible metabolites.
4. Re-run the real panel and freeze generation 02 only after contract tests and appraisal pass.

