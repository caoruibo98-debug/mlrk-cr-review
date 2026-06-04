# Round12 Rule Promotion Overlay

Date: 2026-06-04

## Purpose

Round12 converts the Round11 generator-miss diagnosis into a gated rule-promotion manifest.

This round still does **not** add training labels and does **not** modify the production generator. It asks:

> Which recovered rules are safe enough for a temporary generator dry run, and which must remain blocked because provenance is not strong enough?

## Files Created

- `data/curation/rule_promotion_overlay_round12.csv`
- `data/curation/rule_promotion_dryrun_plan_round12.csv`
- `data/curation/round12_rule_promotion_summary.csv`
- `tools/build_round12_rule_promotion_overlay.py`

## Main Result

Round12 reviewed 10 recovered rule candidates for the four known-gold generation misses.

Status counts:

- `blocked_microberx_source_recovery_required`: 7 rules / 3 pairs
- `blocked_manual_source_direction_license_review`: 2 rules / 1 pair
- `dryrun_candidate_source_traceable_rule`: 1 rule / 1 pair

Only one rule is allowed for dry-run overlay:

| pair | rule | source | status |
|---|---|---|---|
| pinoresinol -> lariciresinol | `retrorules_RR-02-0b00df2834267ac5-06-F` | RetroRules / source legacy `MNXR113586_MNXM59759` | dry-run candidate only |

Deployment promotion remains `false` for all rules.

## Why Only One Rule Passed Dry Run

The pinoresinol rule is tied to a RetroRules source reaction family that Round7 already marked as `source_reaction_candidate_verified_not_imported`, and Round11 found exact pair support from PMID 12736449.

It is still only a dry-run candidate because direction, stereochemistry, source license, and final generator behavior must be checked before production promotion.

## Blocked Rules

The following known-gold pairs still have generating rules, but those rules are not allowed into the generator overlay:

- `urolithin C -> urolithin A`: MicrobeRX-only source recovery required.
- `2'-fucosyllactose -> L-fucose`: MicrobeRX-only source recovery required and 2'-FL structure mapping remains under review.
- `isoxanthohumol -> 8-prenylnaringenin`: MicrobeRX-only source recovery required.

This is the correct conservative behavior. A rule that can generate the right product is not automatically production-safe.

## Production Meaning

Round12 separates three concepts that were previously too easy to blur:

1. Existing gold label: the pair is already known positive.
2. Rule can generate target: a SMARTS transformation can produce the product.
3. Rule is production-safe: source, direction, license, exactness, and split/holdout gates have passed.

Only the third one belongs in a production generator.

## Round13 Direction

Round13 should perform a targeted dry run for the single approved overlay candidate:

1. Add `retrorules_RR-02-0b00df2834267ac5-06-F` to a temporary generator overlay.
2. Rerun clean candidate generation for `pinoresinol` only.
3. Confirm that `lariciresinol` appears in the generated candidate set.
4. Keep the label unchanged as existing gold.
5. Record whether the generator recall repair works before considering broader overlay promotion.

For the other three pairs, Round13 should search replacement source-traceable rules in Rhea/MetaNetX/RetroRules/BRENDA or exact papers before any overlay.

## Claim Boundary

Allowed after Round12:

> One RetroRules source-traceable rule was identified as a dry-run generator overlay candidate for pinoresinol -> lariciresinol. Seven MicrobeRX-only rules remain blocked by source recovery requirements.

Forbidden:

> The production generator has been fixed, or the training set has new labels.

