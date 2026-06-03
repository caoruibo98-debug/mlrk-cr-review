from __future__ import annotations

import importlib
import inspect
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))


TEST_MODULES = [
    "tests.test_manifest_contract",
    "tests.test_prediction_contract",
    "tests.test_api_contract",
    "tests.test_biosanity",
    "tests.test_appraisal_matching",
]


def main() -> int:
    failures: list[str] = []
    total = 0
    for module_name in TEST_MODULES:
        module = importlib.import_module(module_name)
        for name, fn in inspect.getmembers(module, inspect.isfunction):
            if not name.startswith("test_"):
                continue
            total += 1
            try:
                fn()
            except Exception as exc:
                failures.append(f"{module_name}.{name}: {type(exc).__name__}: {exc}")
    if failures:
        print(f"FAILED {len(failures)}/{total} contract tests")
        for item in failures:
            print(item)
        return 1
    print(f"PASSED {total} contract tests")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
