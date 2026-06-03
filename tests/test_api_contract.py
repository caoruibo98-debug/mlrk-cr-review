from __future__ import annotations

import time
import warnings

warnings.filterwarnings("ignore", message="Using `httpx` with `starlette.testclient` is deprecated.*")
from fastapi.testclient import TestClient

from mlrk_prod.api import ApiError, JobRecord, SubmitPrediction, _jobs, _lock, _set_job, app
from mlrk_prod.predict_runner import PredictionRequest, validate_request


client = TestClient(app)


def reset_jobs() -> None:
    with _lock:
        _jobs.clear()


def test_prediction_request_requires_exactly_one_identifier() -> None:
    try:
        validate_request(PredictionRequest())
    except ValueError:
        pass
    else:
        raise AssertionError("empty request should fail")


def test_prediction_request_accepts_name() -> None:
    validate_request(PredictionRequest(name="rutin", topn=8))


def test_prediction_request_rejects_large_topn() -> None:
    try:
        validate_request(PredictionRequest(name="rutin", topn=100))
    except ValueError:
        pass
    else:
        raise AssertionError("topn > 50 should fail")


def test_api_rejects_missing_identifier_with_stable_error_code() -> None:
    reset_jobs()
    response = client.post("/api/jobs", json={"topn": 5})
    assert response.status_code == 400
    assert response.json()["detail"]["code"] == "invalid_request"


def test_api_schema_validation_uses_stable_error_code_without_raw_input() -> None:
    reset_jobs()
    response = client.post("/api/jobs", json={"name": "rutin", "topn": 0})
    payload = response.json()
    assert response.status_code == 422
    assert payload["detail"]["code"] == "invalid_request_schema"
    assert payload["validation_errors"]
    assert "input" not in payload["validation_errors"][0]


def test_api_rejects_malformed_job_id_before_lookup() -> None:
    response = client.get("/api/jobs/not-a-real-job/result")
    assert response.status_code == 400
    assert response.json()["detail"]["code"] == "invalid_job_id"


def test_api_pending_result_uses_retryable_202_contract() -> None:
    reset_jobs()
    now = time.time()
    job_id = "a" * 32
    _set_job(
        JobRecord(
            job_id=job_id,
            status="queued",
            created_at=now,
            updated_at=now,
            request=SubmitPrediction(name="rutin", topn=5),
        )
    )
    response = client.get(f"/api/jobs/{job_id}/result")
    payload = response.json()["detail"]
    assert response.status_code == 202
    assert payload["code"] == "job_not_ready"
    assert payload["retryable"] is True
    assert payload["job_status"] == "queued"


def test_api_failed_result_does_not_expose_internal_exception_text() -> None:
    reset_jobs()
    now = time.time()
    job_id = "b" * 32
    _set_job(
        JobRecord(
            job_id=job_id,
            status="failed",
            created_at=now,
            updated_at=now,
            request=SubmitPrediction(name="rutin", topn=5),
            error=ApiError(
                code="prediction_failed",
                message="prediction failed inside the internal worker; inspect server logs before making biological claims.",
                retryable=False,
                job_status="failed",
            ),
        )
    )
    response = client.get(f"/api/jobs/{job_id}/result")
    payload = response.json()["detail"]
    assert response.status_code == 500
    assert payload["code"] == "prediction_failed"
    assert "Traceback" not in payload["message"]
    assert "D:\\" not in payload["message"]


def test_api_timed_out_result_uses_gateway_timeout_contract() -> None:
    reset_jobs()
    now = time.time()
    job_id = "c" * 32
    _set_job(
        JobRecord(
            job_id=job_id,
            status="timed_out",
            created_at=now,
            updated_at=now,
            request=SubmitPrediction(name="rutin", topn=5),
            error=ApiError(
                code="prediction_timeout",
                message="prediction exceeded the 180 second API timeout.",
                retryable=True,
                job_status="timed_out",
            ),
        )
    )
    response = client.get(f"/api/jobs/{job_id}/result")
    payload = response.json()["detail"]
    assert response.status_code == 504
    assert payload["code"] == "prediction_timeout"
    assert payload["retryable"] is True
