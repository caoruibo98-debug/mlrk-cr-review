# Biochemical Logic Plugin Review

Branch: `biochem-logic-review`

## Scope

This review used the Life Science Research evidence lane for metabolomics and microbiome context. The question was whether the current model's biochemical logic is suitable for production movement, and whether any code change was needed before stronger claims.

## Evidence Check

The model's broad reaction families are biologically plausible for dietary polyphenol gut metabolism:

- Flavonoid conversion by intestinal bacteria includes O- and C-deglycosylation, demethylation, dehydroxylation, ester cleavage, C=C reduction, ring fission, chain modification, and decarboxylation.
- Hydroxycinnamate metabolism can include bacterial reduction of p-coumaric, caffeic, and ferulic acids to substituted phenylpropionic acids.
- Ellagic acid and ellagitannin metabolism to urolithins is a multi-step gut microbial process involving lactone-ring cleavage, decarboxylation, and successive dehydroxylations.
- Recent urolithin work supports that some dehydroxylation steps can be enzyme- and organism-specific, so generic SMARTS rules should be treated as hypotheses unless evidence is attached.

Representative sources checked:

- https://pmc.ncbi.nlm.nih.gov/articles/PMC4939924/
- https://pmc.ncbi.nlm.nih.gov/articles/PMC6052270/
- https://pmc.ncbi.nlm.nih.gov/articles/PMC3679724/
- https://pmc.ncbi.nlm.nih.gov/articles/PMC12771579/

## Finding

The reaction-space logic is directionally appropriate for an internal research MVP, but the quality-tier semantics were too optimistic:

- A structurally plausible candidate with no enzyme, microbe, EC, PMID, or known-product evidence could still receive tier `high`.
- Broad microbial rescue rules, especially dehydroxylation, were not explicitly marked as pathway hypotheses.
- This could make users confuse "chemically plausible generated candidate" with "biologically evidenced gut microbial metabolite."

## Change Made

- Structural-only candidates are now capped at tier `medium`.
- Evidence-backed candidates can still reach tier `high`.
- Prediction payloads now expose `evidence_level`:
  - `traceable_biochemical_evidence`
  - `structural_only`
- Broad microbial rescue candidates without evidence receive `heuristic_microbial_rescue`.
- Phenol dehydroxylation rescue candidates without evidence also receive `broad_dehydroxylation_hypothesis`.
- The model card now states this quality boundary explicitly.

## Production Meaning

This change does not improve rank metrics by hiding failures or changing the expected-product benchmark. It improves claim safety. The model remains suitable for internal prioritization and LC-MS target triage, but structural-only outputs must not be presented as verified biological occurrence.
