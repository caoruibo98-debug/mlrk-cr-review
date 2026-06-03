# Generation 14 Review

## Production question

Can API contracts expose safe error states and timeouts so an internal web app can handle failures without parsing raw exception strings?

## Changes

1. Added `API_CONTRACT_VERSION` and returned it from `/health`, `/readiness`, and job submission.
2. Added a stable `ApiError` object with `code`, `message`, `retryable`, and optional `job_status`.
3. Added job id shape validation before lookup.
4. Added structured handling for:
   - invalid request semantics,
   - schema validation errors,
   - queued/running jobs,
   - failed jobs,
   - timed-out jobs.
5. Added API timeout status support with `prediction_timeout` and HTTP `504`.
6. Prevented raw internal exception strings from being exposed through failed job result responses.
7. Added `docs/API_CONTRACT.md`.
8. Expanded API contract tests.

## Results

- Contract tests: `PASSED 39`
- Production scorecard status: `internal_mvp_only`
- Core strict top-5 hit rate: `6 / 6`
- Challenge strict top-5 hit rate: `8 / 22`
- Readiness: `internal_mvp_only`

## Endpoint behavior now covered

| Scenario | HTTP | Code |
| --- | --- | --- |
| Missing name/smiles | `400` | `invalid_request` |
| Schema violation such as `topn=0` | `422` | `invalid_request_schema` |
| Malformed job id | `400` | `invalid_job_id` |
| Queued or running result | `202` | `job_not_ready` |
| Worker failure | `500` | `prediction_failed` |
| Worker timeout | `504` | `prediction_timeout` |

## Audit interpretation

This generation improves the deployment surface but does not improve biological model capability. The correct claim is:

> The API is safer for internal productization because clients can distinguish invalid input, pending jobs, failed jobs, and timeout states without parsing internal exceptions.

The current production boundary remains unchanged: internal research MVP only.

## Next round

Generation 15 should standardize repository layout and entrypoints without breaking legacy scripts, so the project looks and behaves more like a maintained human-authored repository.
