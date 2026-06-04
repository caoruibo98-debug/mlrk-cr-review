# Round8 Exact-Pair Evidence Extraction

Date: 2026-06-04

## Purpose

Round8 moves from "the paper exists" to "does the source support this exact substrate-product pair?"

This round still does **not** add rows to the training set. It creates a gated manifest of candidate positives and conditional negatives that can be imported only after chemical mapping, exact source extraction, license checks, and split/holdout guards.

## Files Created

- `data/curation/exact_pair_evidence_round8.csv`
- `data/curation/positive_sample_expansion_round8.csv`
- `data/curation/negative_evidence_candidates_round8.csv`
- `data/curation/round8_pipeline_queue.csv`
- `data/curation/ncbi_pubmed_efetch_round8_raw.xml`
- `data/curation/ncbi_pubmed_efetch_round8_extra_raw.xml`
- `data/curation/pubmed_phloroglucinol_round8_raw.json`
- `data/curation/rhea_phloroglucinol_round8_raw.json`

## Main Result

Round8 reviewed 14 source-pair rows.

Decision counts:

- `exact_positive_candidate_pending_mapping_holdout`: 4
- `exact_positive_candidate_better_source_available`: 1
- `fulltext_table_required_before_positive`: 1
- `fulltext_scheme_required_before_positive`: 1
- `context_evidence_not_global_positive`: 1
- `do_not_train_positive_current_evidence`: 4
- `do_not_use_as_primary_training_source`: 1
- `reject_as_positive_current_evidence`: 1

The final positive expansion table keeps only 4 unique best-source exact positive candidates. All remain `training_allowed_round8=false`.

## Exact Positive Candidates

These are the best Round8 positive candidates, pending chemical mapping and leakage controls:

| substrate | product | PMID | evidence class | status |
|---|---|---:|---|---|
| urolithin C | urolithin A | 39856097 | abstract exact pair strong | candidate, not training yet |
| 2'-fucosyllactose | L-fucose | 31138818 | abstract exact pair strong | candidate, not training yet |
| pinoresinol | lariciresinol | 12736449 | abstract exact pair strong | candidate, not training yet |
| isoxanthohumol | 8-prenylnaringenin | 16772450 | abstract exact pair strong | candidate, not training yet |

Important: these are not final training rows yet. They still need:

- exact source extraction from the paper text or table;
- substrate/product SMILES and InChIKey mapping;
- direction and stereochemistry;
- organism/enzyme/assay context;
- source license note;
- holdout/split guard.

## Rejected Or Blocked Positives

- `caffeic acid -> dihydrocaffeic acid`: current source supports a rosmarinic-acid pathway yielding DHCA, not the direct CA -> DHCA edge.
- `urolithin A -> urolithin B`: current PMID supports urolithin/isourolithin context, not this exact pair.
- `urolithin C -> isourolithin A`: current PMID supports isourolithin A from ellagic acid, not the exact intermediate pair.
- `quercetin -> taxifolin`: current source is taxifolin binding/conversion chemistry, not quercetin -> taxifolin.
- `phloroglucinol carboxylic acid -> phloroglucinol`: Round8 PubMed and Rhea narrow searches returned zero records. This is not negative proof, but it is not training evidence.

## Full-Text/Table Needed

- `ethyl ferulate -> ferulic acid`: PubMed abstract supports ethyl ferulate as substrate and ferulic-acid esterase activity, but the exact product mapping should be extracted from the enzyme substrate table/full text.
- `lariciresinol -> secoisolariciresinol`: PubMed abstract says a metabolic scheme exists, and Rhea has a related lignan redox reaction. Direction and stereochemistry must be reconciled before positive import.

## Negative Evidence

Round8 found 2 conditional negative candidates, but zero global negatives.

- PMID 19520709: enzyme-context non-activity for one fucosidase class on alpha1,2-fucosyl substrates. This is a conditional assay-negative candidate, not a global statement that 2'-FL cannot release fucose.
- PMID 40528807: rosmarinic acid remained stable in defined gut/antibiotic contexts. This is a context limit, not a caffeic-acid negative.

Rule:

> A true negative must be tied to a tested substrate, product or product class, organism/enzyme, assay condition, and detection limit. Database absence remains unlabeled.

## Why This Supports Ray's Diagnosis

The blocker is not merely total sample count. The project has many raw positives, but production depends on whether the exact reaction family can be generated and whether the label is defensible.

Round8 shows that after strict evidence filtering, only 4 source-backed exact-positive candidates are ready for the next curation step. That confirms the practical gap:

- many candidate pairs are family/pathway/context evidence only;
- MicrobeRX rules often need source recovery;
- exact positive labels require manual extraction;
- negatives are much harder than decoys.

## Round9 Pipeline

Round9 should do:

1. For the 4 exact candidates, map substrate/product structures to canonical SMILES and InChIKey.
2. Extract one source-backed row per pair with PMID, organism/enzyme, assay context, direction, and stereochemistry.
3. Check whether each pair is already in `reactions.parquet`, `clean_candidates_full.parquet`, or held out.
4. For ethyl ferulate and lariciresinol, extract full-text/table evidence before any label decision.
5. Build `positive_sample_import_round9.csv` with `training_allowed=false` until split guards pass.
6. Keep negative candidates separate as `conditional_negative`, never global negatives.

## Production Claim Boundary

Allowed claim after Round8:

> We identified four literature-supported exact substrate-product candidate positives for further curation and found two conditional negative evidence candidates. No new training labels were imported without source extraction and leakage controls.

Forbidden claim:

> The model now has validated production-level positive and negative training samples.

