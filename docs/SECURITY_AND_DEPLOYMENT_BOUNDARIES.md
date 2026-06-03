# Security And Deployment Boundaries

## Allowed deployment

Internal research MVP only.

## Blocked deployment

Public web app is blocked until external validation, rate limiting, input hardening, and monitoring are complete.

## Required controls

- Maximum SMILES length.
- Maximum reaction enumeration time.
- Per-job timeout.
- Per-user rate limit.
- No arbitrary file paths from user input.
- No raw stack traces in API responses.
- Model and rule artifacts loaded from a pinned manifest.
- Structured logs without raw secrets or personal data.

## Scientific safety controls

- Every result must state that the score is not a wet-lab probability.
- Every no-candidate result must distinguish generation failure from ranking failure.
- Evidence overlay must be labeled as prior evidence, not proof of occurrence.
- Health-effect claims must be excluded from MVP output.

