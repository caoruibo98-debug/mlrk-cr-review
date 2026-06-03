# Internal API Contract

The API is an internal research MVP surface. It is not approved for public consumer, clinical, or wet-lab probability claims.

## Contract version

Current contract version:

```text
2026-06-03.generation14
```

The version appears in `/health`, `/readiness`, and job submission responses.

## Endpoints

```text
GET  /health
GET  /readiness
POST /api/jobs
GET  /api/jobs/{job_id}
GET  /api/jobs/{job_id}/result
```

## Submit request

```json
{
  "name": "rutin",
  "topn": 5
}
```

or:

```json
{
  "smiles": "O=C1...",
  "topn": 5
}
```

Exactly one of `name` or `smiles` is required.

## Submit response

```json
{
  "job_id": "32_hex_characters",
  "status": "queued",
  "api_contract_version": "2026-06-03.generation14",
  "poll_url": "/api/jobs/{job_id}",
  "result_url": "/api/jobs/{job_id}/result",
  "timeout_seconds": 180
}
```

## Stable error object

API errors return a structured `detail` object:

```json
{
  "code": "job_not_ready",
  "message": "job result is not ready yet.",
  "retryable": true,
  "job_status": "queued"
}
```

Known error codes:

| HTTP | Code | Meaning |
| --- | --- | --- |
| 400 | `invalid_request` | Request is syntactically valid but violates model input rules. |
| 400 | `invalid_job_id` | Job id is not a 32-character hexadecimal id. |
| 404 | `job_not_found` | Job id is valid-shaped but not known to the process. |
| 202 | `job_not_ready` | Job exists but has not succeeded yet. Retry polling. |
| 422 | `invalid_request_schema` | Request body failed API schema validation. |
| 500 | `prediction_failed` | Internal worker failed. Do not expose raw exception text to users. |
| 504 | `prediction_timeout` | Prediction exceeded the API timeout. Retryable for internal users. |

## Product boundary

`/readiness` always returns:

```json
{
  "public_web_allowed": false,
  "claim_boundary": "Internal research MVP only; not a consumer, clinical, or wet-lab probability service."
}
```

The web app should display predictions as ranked research candidates, not as occurrence probabilities, health effects, or clinical advice.
