# Reaction Family Expansion Status

## Version Intent

- Branch intent: `reaction-family-expansion`
- Seed manifest: `data/reaction_family_gap_seed_manifest.tsv`
- KPI source: `outputs/appraisal/reaction_family_kpis.json`
- Challenge holdouts are training data: `false`

## Seed Summary

- Seed rows: `14`
- Seeded families: `12`
- Training-allowed rows: `0`
- Holdout rows: `14`
- Blocked families with no seed: `[]`
- Ready for retraining families: `[]`

## Family Status

| Reaction family | Holdout seeds | Training seeds | Target independent training pairs | Additional needed | Current status |
| --- | ---: | ---: | ---: | ---: | --- |
| ellagitannin_hydrolysis | 1 | 0 | 20 | 20 | seeded_holdout_only |
| ellagitannin_urolithin_multistep | 1 | 0 | 30 | 30 | seeded_holdout_only |
| flavanol_ring_fission | 2 | 0 | 30 | 30 | seeded_holdout_only |
| flavanone_ring_fission | 1 | 0 | 25 | 25 | seeded_holdout_only |
| flavonol_ring_fission | 1 | 0 | 25 | 25 | seeded_holdout_only |
| isoflavone_c_glycoside_conversion | 1 | 0 | 20 | 20 | seeded_holdout_only |
| isoflavone_reduction | 2 | 0 | 25 | 25 | seeded_holdout_only |
| isoflavone_reductive_metabolism | 1 | 0 | 30 | 30 | seeded_holdout_only |
| lignan_deglucosylation | 1 | 0 | 20 | 20 | seeded_holdout_only |
| lignan_multistep_demethylation_dehydroxylation | 1 | 0 | 30 | 30 | seeded_holdout_only |
| lignan_oxidation | 1 | 0 | 15 | 15 | seeded_holdout_only |
| phenolic_ester_hydrolysis | 1 | 0 | 20 | 20 | seeded_holdout_only |

## Policy

Seed rows resolve known blocked benchmark cases into machine-readable structures, but must be replaced by independent cases before training to avoid leakage.

## Next Steps

- Collect independent non-holdout positive pairs for each seeded family.
- Attach source IDs, English names, Chinese names, SMILES, InChIKey, PMIDs/PMCIDs/DOIs, enzymes, microbes, and evidence tiers.
- Add reaction-family candidate-generation rules only after at least one independent training seed and one holdout case exist for that family.
- Benchmark ChemBERTa-2 against MoLFormer and MolT5/ChemT5-style encoders after the new cases are versioned.
- Retrain only after data audit confirms no challenge-holdout leakage.
