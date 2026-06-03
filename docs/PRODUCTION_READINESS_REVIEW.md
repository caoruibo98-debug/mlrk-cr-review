# Production Readiness Review

## Current production boundary

This clone is a production-candidate shell around the existing modular L-RCLSS ranker. It does not change the original model logic. The allowed use is internal research prioritization of rule-generated candidate metabolites.

## Current system type

Hybrid research system:

1. Rule-based candidate generation through reaction SMARTS and RDKit.
2. Machine-learning ranking through LTR_chem.
3. Evidence overlay from known enzyme, microbe, and literature records.

It is not a validated biological occurrence predictor.

## Strengths

- Candidate generation and ranking are separated.
- The final ranking path avoids EC and evidence features.
- Clean evaluation includes random, EC-only, and Tanimoto probes.
- The output already includes an honest note that scores are not wet-lab probabilities.

## Production blockers

1. External head-to-head benchmark is not complete.
2. Generation recall remains the hard ceiling.
3. Model artifacts and metrics were not previously bound by a single manifest.
4. Web/API error contracts are not yet formalized.
5. Security and timeout controls are not implemented.
6. No public claims should be made without wet-lab or independent validation.

## Readiness target

The next realistic target is an internal research MVP:

- single-compound query,
- job-based execution,
- ranked table output,
- JSON/CSV export,
- evidence overlay,
- explicit claim boundaries.

