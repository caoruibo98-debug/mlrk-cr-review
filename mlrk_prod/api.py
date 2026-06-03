from __future__ import annotations

import threading
import time
from dataclasses import asdict
from subprocess import TimeoutExpired
from typing import Literal
from uuid import uuid4

from fastapi import FastAPI, HTTPException, Request
from fastapi.exceptions import RequestValidationError
from fastapi.responses import JSONResponse
from pydantic import BaseModel, Field

from .manifest import readiness_status, validate_readiness
from .predict_runner import DEFAULT_TIMEOUT_SECONDS, PredictionRequest, run_prediction, validate_request


API_CONTRACT_VERSION = "2026-06-03.generation14"
JOB_ID_LENGTH = 32
JOB_TIMEOUT_SECONDS = DEFAULT_TIMEOUT_SECONDS


class SubmitPrediction(BaseModel):
    name: str | None = Field(default=None, max_length=120)
    smiles: str | None = Field(default=None, max_length=1500)
    topn: int = Field(default=12, ge=1, le=50)


class ApiError(BaseModel):
    code: str
    message: str
    retryable: bool = False
    job_status: str | None = None


class JobRecord(BaseModel):
    job_id: str
    status: Literal["queued", "running", "succeeded", "failed", "timed_out"]
    created_at: float
    updated_at: float
    request: SubmitPrediction
    result: dict | None = None
    error: ApiError | None = None
    elapsed_seconds: float | None = None


app = FastAPI(
    title="ML Ranking Kernel Internal MVP",
    version="0.1.0",
    description="Internal research MVP for ranking rule-generated food microbiome metabolite candidates.",
)

_jobs: dict[str, JobRecord] = {}
_lock = threading.Lock()


def _api_error(
    status_code: int,
    code: str,
    message: str,
    *,
    retryable: bool = False,
    job_status: str | None = None,
) -> HTTPException:
    detail = ApiError(code=code, message=message, retryable=retryable, job_status=job_status)
    return HTTPException(status_code=status_code, detail=detail.model_dump(exclude_none=True))


def _validate_job_id(job_id: str) -> None:
    is_hex = all(c in "0123456789abcdef" for c in job_id.lower())
    if len(job_id) != JOB_ID_LENGTH or not is_hex:
        raise _api_error(400, "invalid_job_id", "job_id must be a 32-character hexadecimal id.")


def _schema_errors(exc: RequestValidationError) -> list[dict[str, str]]:
    errors = []
    for item in exc.errors():
        errors.append(
            {
                "loc": ".".join(str(part) for part in item.get("loc", [])),
                "message": str(item.get("msg", "invalid value")),
                "type": str(item.get("type", "validation_error")),
            }
        )
    return errors


@app.exception_handler(RequestValidationError)
async def request_validation_exception_handler(_request: Request, exc: RequestValidationError) -> JSONResponse:
    detail = ApiError(
        code="invalid_request_schema",
        message="request body failed API schema validation.",
        retryable=False,
    ).model_dump(exclude_none=True)
    return JSONResponse(status_code=422, content={"detail": detail, "validation_errors": _schema_errors(exc)})


def _set_job(job: JobRecord) -> None:
    with _lock:
        _jobs[job.job_id] = job


def _get_job(job_id: str) -> JobRecord:
    _validate_job_id(job_id)
    with _lock:
        job = _jobs.get(job_id)
    if job is None:
        raise _api_error(404, "job_not_found", "job not found.")
    return job


def _run_job(job_id: str) -> None:
    job = _get_job(job_id)
    job.status = "running"
    job.updated_at = time.time()
    _set_job(job)
    try:
        req = PredictionRequest(name=job.request.name, smiles=job.request.smiles, topn=job.request.topn)
        run = run_prediction(req, timeout_seconds=JOB_TIMEOUT_SECONDS)
        job.status = "succeeded"
        job.result = {
            "payload": run.payload,
            "output_path": str(run.output_path),
            "elapsed_seconds": run.elapsed_seconds,
        }
        job.elapsed_seconds = run.elapsed_seconds
        job.updated_at = time.time()
    except TimeoutExpired:
        job.status = "timed_out"
        job.error = ApiError(
            code="prediction_timeout",
            message=f"prediction exceeded the {JOB_TIMEOUT_SECONDS} second API timeout.",
            retryable=True,
            job_status="timed_out",
        )
        job.updated_at = time.time()
    except Exception:
        job.status = "failed"
        job.error = ApiError(
            code="prediction_failed",
            message="prediction failed inside the internal worker; inspect server logs before making biological claims.",
            retryable=False,
            job_status="failed",
        )
        job.updated_at = time.time()
    _set_job(job)


@app.get("/health")
def health() -> dict:
    return {"status": "ok", "api_contract_version": API_CONTRACT_VERSION}


@app.get("/readiness")
def readiness() -> dict:
    issues = validate_readiness()
    return {
        "status": readiness_status(issues),
        "issues": [asdict(i) for i in issues],
        "api_contract_version": API_CONTRACT_VERSION,
        "public_web_allowed": False,
        "claim_boundary": "Internal research MVP only; not a consumer, clinical, or wet-lab probability service.",
    }


@app.post("/api/jobs")
def submit_prediction(req: SubmitPrediction) -> dict:
    try:
        validate_request(PredictionRequest(name=req.name, smiles=req.smiles, topn=req.topn))
    except ValueError as exc:
        raise _api_error(400, "invalid_request", str(exc)) from exc
    job_id = uuid4().hex
    now = time.time()
    job = JobRecord(job_id=job_id, status="queued", created_at=now, updated_at=now, request=req)
    _set_job(job)
    thread = threading.Thread(target=_run_job, args=(job_id,), daemon=True)
    thread.start()
    return {
        "job_id": job_id,
        "status": job.status,
        "api_contract_version": API_CONTRACT_VERSION,
        "poll_url": f"/api/jobs/{job_id}",
        "result_url": f"/api/jobs/{job_id}/result",
        "timeout_seconds": JOB_TIMEOUT_SECONDS,
    }


@app.get("/api/jobs/{job_id}")
def get_job(job_id: str) -> dict:
    return _get_job(job_id).model_dump()


@app.get("/api/jobs/{job_id}/result")
def get_job_result(job_id: str) -> dict:
    job = _get_job(job_id)
    if job.status == "timed_out":
        detail = job.error or ApiError(
            code="prediction_timeout",
            message="prediction exceeded the API timeout.",
            retryable=True,
            job_status="timed_out",
        )
        raise HTTPException(status_code=504, detail=detail.model_dump(exclude_none=True))
    if job.status == "failed":
        detail = job.error or ApiError(
            code="prediction_failed",
            message="prediction failed inside the internal worker.",
            retryable=False,
            job_status="failed",
        )
        raise HTTPException(status_code=500, detail=detail.model_dump(exclude_none=True))
    if job.status != "succeeded":
        raise _api_error(
            202,
            "job_not_ready",
            "job result is not ready yet.",
            retryable=True,
            job_status=job.status,
        )
    return job.result or {}
