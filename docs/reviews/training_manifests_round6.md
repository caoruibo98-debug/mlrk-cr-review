# Training Manifests Round 6

## Purpose

Round 6 converts Round 5 rule recovery into production-gated training manifests.

The key rule remains:

```text
Reaction rule recovery is candidate generation.
It is not a positive biological label by itself.
```

No model training rows were promoted in this round. The manifests define what
must be verified before a row can enter training.

## Artifacts

- `tools/build_round6_training_manifests.py`
- `data/curation/rule_import_manifest_round6.csv`
- `data/curation/positive_sample_expansion_manifest_round6.csv`
- `data/curation/negative_label_policy_round6.csv`
- `data/curation/round6_pipeline_queue.csv`

## Rule Import Manifest

`rule_import_manifest_round6.csv` has 11 rows, one per recovered
substrate-product pair.

Status counts:

- `candidate_after_source_reaction_check`: 4
- `source_recovery_required_before_import`: 3
- `blocked_until_external_evidence_or_rule_provenance`: 4

All rows have:

```text
import_allowed_round6 = false
training_allowed_round6 = false
```

This is intentional. The next gate is source-reaction verification, not model
training.

Interpretation:

- `candidate_after_source_reaction_check`: RetroRules/EC provenance exists, but
  direction, license, source reaction, and exact biological evidence must still
  be checked.
- `source_recovery_required_before_import`: PubMed evidence may support the
  pair, but the recovered rule is MicrobeRX-only or lacks EC/source-reaction
  provenance.
- `blocked_until_external_evidence_or_rule_provenance`: either exact evidence or
  rule provenance is still too weak.

## Positive Sample Manifest

`positive_sample_expansion_manifest_round6.csv` has 13 evidence rows.

Status counts:

- `eligible_after_exact_extraction_and_rule_provenance_recovery`: 7
- `eligible_after_exact_extraction_and_rule_source_check`: 3
- `mechanism_support_only_not_training_positive`: 2
- `not_training_positive_in_round6`: 1

All rows have:

```text
training_allowed_round6 = false
```

This prevents a common mistake: a PubMed title or rule hit may support a family,
but training requires exact table/text extraction, compound mapping, rule
provenance, holdout status, and bilingual name fields.

Minimum missing fields before training:

- exact table or text extraction
- compound mapping confirmation
- source reaction or rule provenance
- holdout check
- Chinese name for later product/web-app use

## Negative Label Policy

`negative_label_policy_round6.csv` defines five label types.

Allowed:

- `hard_decoy`: allowed for learning-to-rank contrast, but it is not a true
  biological negative.
- `conditional_negative`: allowed only in the same organism, strain, model, or
  community context.
- `assay_negative`: strongest negative type, allowed when a paper reports no
  conversion under explicit conditions.

Not allowed:

- `unlabeled`: not trainable as a negative.
- `prohibited_not_found_negative`: forbidden. "Not found in PubMed/Rhea/VMH" is
  not a negative label.

This is the production-safe label boundary. The model can use hard decoys for
ranking, but claims must say "ranked against generated alternatives", not
"experimentally false alternatives".

## Pipeline Queue

`round6_pipeline_queue.csv` has 24 queued actions:

- 11 rule-import actions
- 13 positive-curation actions

High priority work:

1. Verify RetroRules source reactions for the 4 `candidate_after_source_reaction_check` rows.
2. Recover source reaction/EC provenance for MicrobeRX-only rows with PubMed support.
3. Extract exact substrate-product evidence from PubMed tables/text.
4. Keep holdout pairs blocked from training.
5. Only after those gates pass, rebuild `clean_candidates_full.parquet` with deterministic evidence-tiered rule inclusion.

## Current Answer To Ray's Question

Yes, sample type and reaction type coverage are real blockers, but the sharper
diagnosis is:

```text
The model has many rows, but too few production-eligible rows per real food-gut
reaction family, and the clean rule-generation route drops many supported
families before the ML ranker can learn from them.
```

The next implementation change should be a rule-import and candidate-generation
change, not blind sample expansion or retraining.

## Production Gate

A row becomes trainable only when all are true:

1. exact substrate-product pair is source-backed;
2. compound structures and InChIKeys are confirmed;
3. reaction family/type is assigned;
4. rule/provenance is retained if used for generation;
5. holdout leakage guard passes;
6. label type is explicit;
7. negative labels are not inferred from database absence.

Until then, the row stays in curation or candidate-generation, not supervised
training.
