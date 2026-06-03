# Web App MVP Specification

## Intended users

- Food metabolomics researchers.
- Gut microbiome researchers.
- Internal model developers reviewing candidate quality.

## Primary user flow

1. User submits compound name or SMILES.
2. App validates and canonicalizes the input.
3. App routes the compound to module A/B/C/D.
4. App generates rule-based candidates.
5. App ranks candidates using LTR_chem.
6. App displays candidates with score, structure, evidence, and claim boundary.
7. User exports JSON/CSV for wet-lab planning.

## Required views

- Query form.
- Job status panel.
- Ranked candidates table.
- Candidate detail drawer.
- Evidence ladder.
- Failure explanation for no-candidate cases.
- Download panel.

## Non-goals for MVP

- Consumer health recommendation.
- Clinical interpretation.
- Automatic causal claims.
- Public batch API.

## API shape

```text
POST /api/jobs
GET  /api/jobs/{job_id}
GET  /api/jobs/{job_id}/result
GET  /api/jobs/{job_id}/download.csv
```

Current implementation:

```text
GET  /health
GET  /readiness
POST /api/jobs
GET  /api/jobs/{job_id}
GET  /api/jobs/{job_id}/result
```

CSV download remains a next-step endpoint.

## Current local server

Start:

```bash
python -m mlrk_prod.cli serve --host 127.0.0.1 --port 8765
```

The current API is job-based. A typical rutin job completes in about 12-15 seconds on the tested local machine.

## Result fields

- query name and SMILES
- module and route reason
- number of rule candidates
- rank
- product name
- product SMILES
- normalized rank score
- evidence overlay
- honest note
