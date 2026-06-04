# Master Rule Recovery Round 5

## Question

Round 4 found 15 route-loss rows where `fullrule_pair_hit=1`, but the clean
local positive-pool rule set did not generate the target product.

Round 5 asks:

1. Which exact `reaction_rule_master.csv` rules can regenerate those target pairs?
2. Which recovered rules are safe candidates for rule import?
3. Which pairs still require literature/database verification before becoming positives?

## Artifacts

- `tools/extract_master_rule_recovery_round5.py`
- `tools/build_rule_recovery_decisions_round5.py`
- `data/curation/master_rule_recovery_candidates_round5.csv`
- `data/curation/rule_recovery_decisions_round5.csv`
- `data/curation/literature_verification_round5.csv`

## Master Rule Recovery Result

The 15 Round 4 rows collapse to 11 unique `(module, substrate, product)` pairs.
All 11 pairs had at least one master rule that regenerated the target product.

Recovered rule rows:

- total recovered rule rows: 41
- unique pairs: 11
- MicrobeRX rules: 21
- RetroRules rules: 20

This confirms that a major blocker is rule routing/import, not only missing
positive labels.

## Import Decisions

`rule_recovery_decisions_round5.csv` separates candidates into two classes:

- `priority_rule_import_candidate`: 4 pairs
- `rule_candidate_needs_source_recovery`: 7 pairs

Priority rule-import candidates have RetroRules source plus EC/source-reaction
metadata:

- caffeic acid -> dihydrocaffeic acid
- pinoresinol -> lariciresinol
- lariciresinol -> secoisolariciresinol
- quercetin -> taxifolin

These still need source-reaction, direction, license, and exact biological
evidence checks before training. The rule can generate the product; it does not
by itself prove a gut microbial reaction.

Source-recovery candidates are MicrobeRX-only in the master rule table and lack
EC/source-reaction metadata:

- ethyl ferulate -> ferulic acid
- urolithin C -> urolithin A
- urolithin A -> urolithin B
- urolithin C -> isourolithin A
- 2'-fucosyllactose -> L-fucose
- isoxanthohumol -> 8-prenylnaringenin
- phloroglucinol carboxylic acid -> phloroglucinol

These must not be imported into the production rule pool until their source
reaction or a replacement database/literature rule is recovered.

## Literature Verification

`literature_verification_round5.csv` records pair-level evidence status.

Strong or promising exact-positive candidates:

- ethyl ferulate -> ferulic acid: PMID:19502437 supports cinnamoyl esterases
  from a Lactobacillus johnsonii stool isolate. Exact substrate table should be
  extracted before training.
- pinoresinol -> lariciresinol: PMID:12736449 directly supports human intestinal
  microflora and Enterococcus faecalis PDG-1 transformation.
- urolithin C -> urolithin A: local Round2 cites PMID:39856097 with direct
  enzyme evidence; re-fetch DOI/PubMed metadata before final import.
- 2'-fucosyllactose -> L-fucose: PMID:31138818, PMID:19520709, and related HMO
  papers support Bifidobacterium fucosylated oligosaccharide metabolism; this
  should be moved to the carbohydrate/glycan family instead of lignan context.

Mechanism/support-only candidates:

- quercetin -> taxifolin: PMID:36432010 supports taxifolin/Eubacterium ramulus
  enzyme context but not yet exact quercetin-to-taxifolin evidence.
- isoxanthohumol -> 8-prenylnaringenin: PMID:20397197 and PMID:41571201 support
  microbiota-derived metabolite context, but exact conversion evidence still
  needs extraction.
- phloroglucinol carboxylic acid -> phloroglucinol: the exact PubMed query
  returned zero records in Round 5; this is not a negative label.

## Negative Sample Boundary

Round 5 did not add biological negatives.

Recovered rules can create same-substrate hard decoys only when:

1. the generated product is not a curated positive for that substrate;
2. it is labeled as `hard_decoy`, not `true_negative`;
3. source, rule, and candidate-generation fields are retained.

Allowed negative-like labels remain:

- `hard_decoy`: generated non-positive candidate for ranking contrast.
- `conditional_negative`: defined organism/context lacks pathway/gene/model
  support.
- `assay_negative`: a paper reports no conversion under explicit conditions.
- `unlabeled`: default for unknowns.

## Production Implication

Round 5 narrows the production blocker:

```text
production blocker =
  missing evidence-backed reaction families
  + missing rule provenance/import for existing master rules
  + random/default clean rule sampling that drops recoverable reactions
```

The next production-aligned code change should be a rule-import manifest, not
model retraining:

1. Build a `rule_import_manifest_round6.csv` with rule source, source reaction,
   EC, direction, evidence tier, and holdout guard.
2. Import only `priority_rule_import_candidate` rows after source-reaction
   verification.
3. Keep MicrobeRX-only rows in `needs_source_recovery` until EC/source reaction
   provenance is recovered.
4. Replace random rule sampling with deterministic evidence-tiered rule
   inclusion.
