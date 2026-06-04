# Round24 Rule-level Coverage Review

## Working conclusion

Round24 moved from reaction-level evidence to rule-level evidence. The main production blocker still sits before ranking: most P0 families do not have a source-backed rule/template that current `run_reactants` can use.

## Rule-level findings

- EnzymeMap `processed_reactions.csv.gz` was downloaded into runtime and screened at exact substrate/product level.
- `protocatechuic acid -> catechol` has EnzymeMap exact pair coverage and rule 92 generates the target in current `run_reactants`.
- Rule 92 is broad aromatic acid decarboxylation and is not promoted because the gut-specific evidence and overgeneration risk are unresolved.
- Urolithin P0 reactions have no Rhea/ECReact/EnzymeMap exact rule-level coverage in this run despite strong recent PubMed support.
- `dopamine -> m-tyramine` has exact Rhea support but no EnzymeMap/ECReact exact pair or rule-level support in this run.
- `corticosterone -> 11beta-hydroxyprogesterone` has an EnzymeMap reverse/context hit from Bos taurus hydroxylase chemistry, but it did not generate the target and is not microbiome-appropriate.

## Training gate

No Round24 sample or rule is admitted to training. The one target-generating rule candidate must first pass biological re-curation, aromatic acid overgeneration review, and false-negative screening.

## Generated artifacts

* `data/curation/round24_external_rule_level_coverage.csv`
* `data/curation/round24_rule_promotion_candidates.csv`
* `data/curation/round24_enzymemap_p0_exact_pair_screen.csv`
* `data/curation/round24_enzymemap_rule_dryrun.csv`
* `data/curation/round24_enzymemap_keyword_screen.csv`
* `data/curation/round24_family_decision_matrix.csv`
* `data/curation/round24_negative_mining_update.csv`
* `data/curation/round24_external_query_manifest.csv`
