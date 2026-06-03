from __future__ import annotations

import threading
import time
from dataclasses import asdict
from typing import Literal
from uuid import uuid4

from fastapi import FastAPI, HTTPException
from pydantic import BaseModel, Field

from .manifest import readiness_status, validate_readiness
from .predict_runner import PredictionRequest, run_prediction, validate_request


class SubmitPrediction(BaseModel):
    name: str | None = Field(default=None, max_length=120)
    smiles: str | None = Field(default=None, max_length=1500)
    topn: int = Field(default=12, ge=1, le=50)


class JobRecord(BaseModel):
    job_id: str
    status: Literal["queued", "running", "succeeded", "failed"]
    created_at: float
    updated_at: float
    request: SubmitPrediction
    result: dict | None = None
    error: str | None = None
    elapsed_seconds: float | None = None


app = FastAPI(
    title="ML Ranking Kernel Internal MVP",
    version="0.1.0",
    description="Internal research MVP for ranking rule-generated food microbiome metabolite candidates.",
)

_jobs: dict[str, JobRecord] = {}
_lock = threading.Lock()


def _set_job(job: JobRecord) -> None:
    with _lock:
        _jobs[job.job_id] = job


def _get_job(job_id: str) -> JobRecord:
    with _lock:
        job = _jobs.get(job_id)
    if job is None:
        raise HTTPException(status_code=404, detail="job not found")
    return job


def _run_job(job_id: str) -> None:
    job = _get_job(job_id)
    job.status = "running"
    job.updated_at = time.time()
    _set_job(job)
    try:
        req = PredictionRequest(name=job.request.name, smiles=job.request.smiles, topn=job.request.topn)
        run = run_prediction(req)
        job.status = "succeeded"
        job.result = {
            "payload": run.payload,
            "output_path": str(run.output_path),
            "elapsed_seconds": run.elapsed_seconds,
        }
        job.elapsed_seconds = run.elapsed_seconds
        job.updated_at = time.time()
    except Exception as exc:
        job.status = "failed"
        job.error = str(exc)
        job.updated_at = time.time()
    _set_job(job)


@app.get("/health")
def health() -> dict:
    return {"status": "ok"}


@app.get("/readiness")
def readiness() -> dict:
    issues = validate_readiness()
    return {"status": readiness_status(issues), "issues": [asdict(i) for i in issues]}


@app.post("/api/jobs")
def submit_prediction(req: SubmitPrediction) -> dict:
    try:
        validate_request(PredictionRequest(name=req.name, smiles=req.smiles, topn=req.topn))
    except ValueError as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc
    job_id = uuid4().hex
    now = time.time()
    job = JobRecord(job_id=job_id, status="queued", created_at=now, updated_at=now, request=req)
    _set_job(job)
    thread = threading.Thread(target=_run_job, args=(job_id,), daemon=True)
    thread.start()
    return {"job_id": job_id, "status": job.status}


@app.get("/api/jobs/{job_id}")
def get_job(job_id: str) -> dict:
    return _get_job(job_id).model_dump()


@app.get("/api/jobs/{job_id}/result")
def get_job_result(job_id: str) -> dict:
    job = _get_job(job_id)
    if job.status == "failed":
        raise HTTPException(status_code=500, detail=job.error)
    if job.status != "succeeded":
        raise HTTPException(status_code=202, detail=job.status)
    return job.result or {}

