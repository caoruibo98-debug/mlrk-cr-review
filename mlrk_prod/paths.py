from __future__ import annotations

import os
from pathlib import Path


def kernel_root() -> Path:
    """Return the production-candidate kernel root.

    MLRK_KERNEL_ROOT can override this for tests or service deployment.
    """
    override = os.environ.get("MLRK_KERNEL_ROOT")
    if override:
        return Path(override).resolve()
    return Path(__file__).resolve().parents[1]


def resolve_kernel_path(relative_path: str | Path) -> Path:
    p = Path(relative_path)
    return p if p.is_absolute() else kernel_root() / p

