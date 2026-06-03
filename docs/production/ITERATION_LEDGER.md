# Production Iteration Ledger

This ledger tracks the minimum 20-round path from research prototype to a GitHub-ready internal production candidate.

## Current Production Claim

The model is currently strongest as an internal research MVP for ranking rule-generated food polyphenol metabolite candidates, especially glycoside-to-aglycone transformations. It is not yet a broad predictor of all food-derived gut microbial metabolites.

## Round Status

| Round | Freeze tag | Status | Main question | Main result | Remaining gap |
| --- | --- | --- | --- | --- | --- |
| 01 | `generation_01` | complete | Can the original model be wrapped with a readiness gate? | Added initial appraisal and readiness framing. | Score and evidence logic were too coarse. |
| 02 | `generation_02` | complete | Can structure-aware matching avoid name-display errors? | Added InChIKey block-1 matching and real panel. | Evidence autonomy was overstated. |
| 03 | `generation_03` | complete | Can evidence sources and interpretation-ready outputs be separated? | Split model-output evidence from benchmark traceability. | Quercitrin still exposed a generation miss. |
| 04 | `generation_04` | complete | Can the hard quercitrin glycoside miss be rescued conservatively? | Added aromatic O-glycoside rescue and reached `4.43 / 5`. | CodeRabbit review still needed. |
| 05 | `generation_05` | complete | Can CodeRabbit findings be fixed without inflating the score? | Fixed appraisal false positives, invalid SMILES handling, rank parsing, and freeze traceability. | Fresh-clone reproducibility failed in mirror. |
| 06 | `generation_06` | complete | Can a clean GitHub mirror reproduce tests and model artifacts? | Versioned deployment artifacts and made prediction contract self-generating. | Evaluation scope still too narrow. |
| 07 | `generation_07` | complete | Can evaluation reveal the next production gap beyond glycoside hydrolysis? | Added challenge panel, production scorecard, no-candidate payload fix, reviewer README, dependency declarations, and `scripts/` entrypoints. | Challenge panel top-5 hit rate was `0 / 6`; next round needed sample expansion and failure taxonomy. |
| 08 | `generation_08` | complete | Can the challenge panel be expanded to at least 20 literature-backed cases and split generation failure from ranking failure? | Expanded challenge panel to `22` cases, added candidate-pool summaries, failure taxonomy, and scorecard failure counts. | Challenge top-5 is `1 / 22`; `21 / 22` failures are generation coverage failures. |
| 09 | `generation_09` | complete | Can CodeRabbit's generation 08 review findings be fixed for reproducibility and strict identity scoring? | Pinned dependencies, added full-InChIKey appraisal fields, split strict vs connectivity top-5, warned on missing freeze paths, removed hidden evidence-pool default, and scoped RDKit logging. | Core score is now `3.88 / 5` because unversioned local evidence is no longer used by default. |
| 10 | `generation_10` | complete | Can candidate generation recover first-wave products for challenge families without changing the scoring rubric? | Added generic microbial rescue rules for hydroxycinnamate reduction, stilbene reduction, aromatic decarboxylation, and phenol dehydroxylation. Challenge top-5 improved from `1 / 22` to `8 / 22`. | Remaining `14 / 22` failures need isoflavone reduction, enterolignan, ring-fission, and ellagitannin/urolithin logic. |
| 11 | `generation_11` | complete | Can external tool inputs be exported reproducibly? | Added BioTransformer TSV/SDF, MicrobeRX query table, GutBug PubChem-ID table, unified expected-label table, manifest, tests, and scorecard status. | External tools still need to be executed and imported for true external validation. |
| 12 | `generation_12` | complete | Can external tool outputs be imported and scored case-by-case? | Added external product/EC result import templates, scoring harness, tests, external-result scorecard, and readiness wording. | Real BioTransformer/MicrobeRX/GutBug outputs still need to be run and imported. |
| 13 | `generation_13` | complete | Can substrate/reaction-family coverage be reported like a product KPI? | Added reaction-family KPI JSON/CSV, tests, scorecard integration, and production-gap recommendations. | Real external results and independent model-output evidence are still missing. |
| 14 | `generation_14` | complete | Can API contracts expose safe error states and timeouts? | Added stable API error codes, timeout status, schema-validation handler, safe failed-job responses, API contract docs, and endpoint-level negative tests. | API is still process-local and not a public web service. |
| 15 | `generation_15` | complete | Can repository layout be standardized without breaking legacy scripts? | Added repository guide, repo doctor command/report, script/tool pairing checks, layout tests, and manifest repository contract. | Legacy code is still intentionally preserved; deeper module migration remains future work. |
| 16 | `generation_16` | complete | Can README become reviewer-first and GitHub-ready? | Reworked README around reviewer summary, fast review path, evidence, outputs, API, limitations, and reviewer docs; added README contract test. | Still needs final release packaging and external benchmark outputs. |
| 17 | pending | planned | Can model cards and claim boundaries be tied to scorecard outputs? | Not started. | Needs automated scorecard-to-doc update path. |
| 18 | pending | planned | Can CodeRabbit or equivalent review be run on the production branch after fixes? | Not started. | May be rate-limited; PR path preferred. |
| 19 | pending | planned | Can a release candidate pass fresh clone, tests, appraisal, and scorecard? | Not started. | Needs clean clone verification. |
| 20 | pending | planned | Can the final GitHub version justify its production boundary? | Not started. | Needs final internal/external/cross/ablation report. |

## Completion Rule For Each Future Round

Each round must include:

1. one specific production-readiness question,
2. one material code, data, evaluation, or documentation improvement,
3. contract tests,
4. relevant appraisal or scorecard run,
5. a review note under `docs/reviews/`,
6. a git commit,
7. a freeze manifest and tag.
