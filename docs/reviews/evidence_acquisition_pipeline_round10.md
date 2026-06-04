# Round10 Evidence Acquisition Pipeline

Date: 2026-06-04

## Direct Answer

Ray's diagnosis is right, but it needs one refinement:

> Reaction-type coverage is a real blocker, but the sharper production blocker is deployment-faithful candidate generation plus source-traceable evidence. Some gold positives already exist, yet the clean candidate generator fails to enumerate them.

That means simply adding more rows is not enough. The next production step must improve:

1. reaction-family coverage;
2. source-reaction provenance;
3. exact positive labels;
4. conditional negative logic;
5. generator recall before LTR retraining.

## What Was Checked

Local evidence:

- `outputs/gen_gap/fullrule_hit_miss.parquet`
- `data/production_challenge_panel.csv`
- `data/curation/positive_sample_import_round9.csv`

External standards:

- ECREACT / RXN biocatalysis model
- RetroPathRL / RetroRules
- gapseq
- gutSMASH

Life-science database screens:

- PubMed ESearch/ESummary via NCBI Entrez
- Rhea reaction search

## Files Created

- `data/curation/reaction_family_gap_priorities_round10.csv`
- `data/curation/external_repository_reaction_standards_round10.csv`
- `data/curation/pubmed_candidate_screen_round10.csv`
- `data/curation/rhea_reaction_screen_round10.csv`
- `data/curation/sample_acquisition_policy_round10.csv`
- `data/curation/evidence_expansion_pipeline_round10.csv`
- `data/curation/round10_next_evidence_queue.csv`
- `data/curation/round10_iteration_log.csv`
- `data/curation/ncbi_round10_raw/*.json`
- `data/curation/rhea_round10_raw/*.json`

## Priority Reaction Families

`reaction_family_gap_priorities_round10.csv` ranks nine families:

1. `polyphenol_ring_fission`
2. `urolithin_dehydroxylation`
3. `hydroxycinnamate_reduction_and_hydrolysis`
4. `isoflavone_reductive_and_glycoside_metabolism`
5. `lignan_redox_and_deglycosylation`
6. `prenylflavonoid_o_demethylation`
7. `glycoside_and_hmo_hydrolysis`
8. `bile_acid_deconjugation_and_lipid_context`
9. `amino_acid_catabolism_and_redox`

Highest-risk local signals:

- `polyphenol_ring_fission`: full-rule recall 0.44, four production challenge cases.
- `urolithin_dehydroxylation`: full-rule recall 0.5455 and known gold pair `urolithin C -> urolithin A` absent from clean candidates.
- `lignan_redox_and_deglycosylation`: known gold pair `pinoresinol -> lariciresinol` absent from clean candidates.
- `prenylflavonoid_o_demethylation`: known gold pair `isoxanthohumol -> 8-prenylnaringenin` absent from clean candidates.
- `glycoside_and_hmo_hydrolysis`: known gold pair `2'-fucosyllactose -> L-fucose` absent from clean candidates.

## External Standards

The GitHub/model comparison supports a strict curation rule:

- ECREACT uses Rhea, BRENDA, PathBank, and MetaNetX and carries `rxn_smiles`, `ec`, and `source`.
- RetroPathRL uses RetroRules for candidate generation; the rule hit is not itself a biological label.
- gapseq uses curated reaction/pathway/genome resources and highlights versioned source and license constraints.
- gutSMASH supports gut anaerobe metabolic potential, but gene-cluster prediction is context evidence, not exact substrate-product proof.

Production implication:

> Do not invent reaction rules. Import only exact source-backed positives and source-traceable generator rules. Treat weak pathway/gene/database context as features or curation pointers, not gold labels.

## PubMed And Rhea Screen

Round10 did not import training labels. It created candidate screens:

- `pubmed_candidate_screen_round10.csv`: 35 PubMed candidate rows.
- `rhea_reaction_screen_round10.csv`: 46 Rhea screen rows.

Important findings:

- Rhea has strong source-rule candidates for chlorogenate/caffeate chemistry, daidzein glycoside chemistry, bile-acid deconjugation, and gallic-acid decarboxylation.
- Rhea is sparse or empty for key gut-specific gaps such as urolithin, equol, HMO hydrolysis, and flavanol/flavonol ring fission.
- PubMed returns many reviews, metabolomics studies, and context papers; these are useful for search direction, but not training labels until exact substrate-product rows are extracted.

## Positive Sample Rule

A future positive row can train only when it has:

- exact substrate and product;
- stable SMILES and InChIKey;
- source PMID/DOI/database/version;
- explicit or defensible direction;
- organism, enzyme, strain, community, or assay context when available;
- evidence tier;
- split/holdout guard;
- generator recall status.

Forbidden as positives:

- rule-fired products with no exact source pair;
- broad review statements;
- pathway-family evidence only;
- inferred "this should happen" chemistry;
- generated candidates with no database or paper support.

## Negative Sample Rule

Good negatives are not "not found in database".

Allowed:

- `assay_negative`: explicit no-conversion/no-product under defined conditions.
- `conditional_negative`: versioned organism/genome/pathway absence under a defined scope.
- `hard_decoy`: generated contrastive candidate for ranking only.

Forbidden:

- Rhea/PubMed no-hit as global negative.
- rule did not fire as biological negative.
- weak organism absence without database version and scope.

## Expansion Pipeline

`evidence_expansion_pipeline_round10.csv` defines the next cycle:

1. collect;
2. verify exactness;
3. map chemistry;
4. split and holdout guard;
5. generator recall check;
6. freeze and review;
7. retrain after gates.

The most important gate is stage 5:

> Do not retrain just because a pair exists. First confirm the deployment generator can enumerate the known product.

`round10_next_evidence_queue.csv` converts the priority families into Round11 work items with seed PMIDs/Rhea IDs/local pair keys, acceptance gates, and rejection rules.

## Next Production Move

Round11 should take one priority family at a time, starting with `polyphenol_ring_fission` and the four Round9 known-gold generation misses.

For each family:

1. extract exact source reactions from PubMed/Rhea/MetaNetX/BRENDA/RetroRules;
2. map structures;
3. recover or import source-traceable rules;
4. rerun clean candidate generation;
5. freeze a commit;
6. audit whether the true product is now generated;
7. only then retrain and evaluate.

## Claim Boundary

Allowed after Round10:

> We identified priority production-blocking reaction families, screened PubMed/Rhea/GitHub sources, and created a gated evidence-acquisition pipeline for exact positives, source-traceable rules, conditional negatives, hard decoys, and unlabeled candidates.

Forbidden:

> We have already expanded the training set with validated new positives or negatives.
