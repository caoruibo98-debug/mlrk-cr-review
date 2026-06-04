# Round13 Overlay Dry Run

Date: 2026-06-04

## Purpose

Round13 executes the single dry-run candidate allowed by Round12:

- Pair: `pinoresinol -> lariciresinol`
- Rule: `retrorules_RR-02-0b00df2834267ac5-06-F`
- Source legacy ID: `MNXR113586_MNXM59759`

This round does not change training labels and does not promote the rule to production. It only tests whether the source-traceable overlay rule can generate the target product block.

## Files Created

- `data/curation/rule_overlay_dryrun_round13.csv`
- `data/curation/round13_generator_repair_decision.csv`
- `tools/run_round13_overlay_dryrun.py`

## Result

The dry run succeeded.

`rule_overlay_dryrun_round13.csv` reports:

- `target_generated=true`
- `dryrun_status=target_generated_by_overlay_rule`
- generated target InChIKey: `MHXCIKYXNYCMHY-AUSJPIAWSA-N`

This matches the target product block for lariciresinol.

## Interpretation

This proves a narrow but important point:

> At least one known-gold clean-generation miss can be repaired by a source-traceable recovered RetroRules rule without adding a duplicate positive label.

The current model's blocker is therefore actionable:

1. the label exists;
2. the source-supported rule exists;
3. the rule can generate the target;
4. the current clean generator does not include that path.

## What Is Still Not Proven

Round13 does not prove:

- production deployment readiness;
- global rule safety;
- correct direction across all substrates;
- leakage-safe evaluation improvement;
- wet-lab validity;
- new training labels.

## Round14 Direction

Round14 should wire this rule into a temporary clean-generator overlay for the pinoresinol substrate only and rerun the clean candidate output.

Success criterion:

> `clean_candidates_full` equivalent output contains `HGXBRUKMWQGOIE__MHXCIKYXNYCMHY` without duplicate-importing the positive label.

After that, the same overlay logic can be generalized only after source/direction/license checks and family-level negative controls.

## Claim Boundary

Allowed after Round13:

> A source-traceable RetroRules dry-run overlay can generate the known gold product lariciresinol from pinoresinol.

Forbidden:

> The production generator has been fixed, or the model is production-ready.

