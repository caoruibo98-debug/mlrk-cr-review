from __future__ import annotations

from mlrk_prod.predict_runner import PredictionRequest, validate_request


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

