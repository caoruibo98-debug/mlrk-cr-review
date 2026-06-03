from __future__ import annotations

import argparse
import json
from dataclasses import asdict

from .manifest import readiness_status, validate_readiness
from .predict_runner import PredictionRequest, run_prediction


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="mlrk-prod")
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("validate-readiness", help="Validate production-candidate artifacts and gates.")
    pred = sub.add_parser("predict", help="Run one internal MVP prediction.")
    pred.add_argument("--name", default=None)
    pred.add_argument("--smiles", default=None)
    pred.add_argument("--topn", type=int, default=12)
    serve = sub.add_parser("serve", help="Start the internal FastAPI MVP.")
    serve.add_argument("--host", default="127.0.0.1")
    serve.add_argument("--port", type=int, default=8765)
    args = parser.parse_args(argv)

    if args.command == "validate-readiness":
        issues = validate_readiness()
        status = readiness_status(issues)
        print(json.dumps({"status": status, "issues": [asdict(i) for i in issues]}, indent=2))
        return 1 if status == "blocked" else 0
    if args.command == "predict":
        run = run_prediction(PredictionRequest(name=args.name, smiles=args.smiles, topn=args.topn))
        print(json.dumps(run.payload, ensure_ascii=False, indent=2))
        return 0
    if args.command == "serve":
        import uvicorn
        uvicorn.run("mlrk_prod.api:app", host=args.host, port=args.port)
        return 0
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
