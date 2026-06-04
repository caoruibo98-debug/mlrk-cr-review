# Round25 rule92 overgeneration and urolithin extraction queue

## Working conclusion

Ray's current judgment is correct in the important sense: the production blocker is not the existence of only 22 samples. The 22-row file is a challenge panel. The blocking issue is earlier in the pipeline: source-backed reaction-family coverage and the label contract for generated negatives.

## Rule 92 audit

- Rule: `[#6:1]-[#6:2]-[#8:3]>>[#6:1].[#6:2]=[#8:3]` from EnzymeMap rule id 92 context (4.1.1.59/4.1.1.61/4.1.1.63).
- Module-substrate rows screened: 5332.
- Unique substrate blocks screened: 4698.
- Module-substrate rows with at least one generated product: 4408.
- Generated substrate-product rows: 39329.
- Generated rows already known as positive: 320 (0.0081 proxy rate; this is not precision because most rows are unlabeled).
- Generated rows overlapping current decoy space: 317.

Interpretation: rule 92 can generate the protocatechuate to catechol product, but as a broad template it creates many unlabeled products. These unlabeled products cannot be hardened as biological negatives.

## Urolithin evidence queue

- Urolithin extraction queue contains 5 high-value source rows; the targeted PubMed Enterocloster query directly returned four dehydroxylase-focused PMIDs.
- Rhea query for `urolithin dehydroxylase` returned zero hits in this round.
- PubChem mappings are present for urolithin A, C, and M6; urolithin G needs a synonym/structure-specific lookup.
- No urolithin row is training-allowed in round25; every row needs exact table extraction and atom mapping first.

## Negative sample policy

- Missing database evidence is not a negative label.
- Rule-generated unknowns can be ranking decoys only if the CSV says `label_scope=decoy_not_biological_negative`.
- Hard negatives require explicit no-conversion assay evidence under a named organism/enzyme/condition.

## Files

- `data/curation/round25_rule92_overgeneration_audit.csv`
- `data/curation/round25_rule92_candidate_decision.csv`
- `data/curation/round25_urolithin_exact_extraction_queue.csv`
- `data/curation/round25_negative_decoy_gate.csv`
- `data/curation/round25_external_query_manifest.csv`

## Sources recorded

- Rhea: https://www.rhea-db.org/rhea/22416 (approved exact protocatechuate decarboxylase reaction)
- Rhea: https://www.rhea-db.org/ (zero exact Rhea hits in round25 query)
- Rhea: https://www.rhea-db.org/ (approved exact dopamine + AH2 = 3-tyramine + A + H2O)
- PubMed: https://pubmed.ncbi.nlm.nih.gov/ (PMIDs 41797252, 41298472, 39856097, 37494568)
- PubChem: https://pubchem.ncbi.nlm.nih.gov/ (PubChem structures found for urolithin A, C, and M6; urolithin G name lookup failed)
- GitHub: https://github.com/rxn4chemistry/biocatalysis-model (ECReact is drawn from Rhea, BRENDA, PathBank, MetaNetX and covers all 7 EC classes.)
- GitHub: https://github.com/hesther/enzymemap (EnzymeMap processed_reactions.csv.gz is the current database file used for rule-level audit.)
- RetroRules: https://retrorules.org/docs (RetroRules uses AAM reaction-center templates and mono-substrate decomposition; not all plausible food-gut reactions are precomputed.)
