"""
kio.py — L-RCLSS 共享 IO / 化学小工具（fresh，无规则特征耦合）。

设计原则：
  - RDKit 仅在**离线**数据准备里用（canonical SMILES / InChIKey / Murcko scaffold）。
    内核**打分路径**不调 RDKit（见 encode_molecules.py / train_kernel.py）。
  - 写出强制落在 ml_ranking_kernel/outputs/ 内（安全契约）。
"""
from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, Iterable

import pandas as pd
import yaml

# === 路径锚点 ===
# src/kio.py → parents[0]=src, [1]=ml_ranking_kernel, [2]=ssrf, [3]=scripts, [4]=Food models
KERNEL_ROOT = Path(__file__).resolve().parents[1]
PROJECT_ROOT = Path(__file__).resolve().parents[4]
OUTPUTS = KERNEL_ROOT / "outputs"
CONFIG_PATH = KERNEL_ROOT / "configs" / "kernel_config.yaml"

DATASETS_DIR = OUTPUTS / "datasets"
FEATURES_DIR = OUTPUTS / "features"
SPLITS_DIR = OUTPUTS / "splits"
MODELS_DIR = OUTPUTS / "models"
METRICS_DIR = OUTPUTS / "metrics"
REPORTS_DIR = OUTPUTS / "reports"
PREDICTIONS_DIR = OUTPUTS / "predictions"


def setup_logging(verbose: bool = False) -> logging.Logger:
    logging.basicConfig(
        level=logging.DEBUG if verbose else logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
        datefmt="%H:%M:%S",
    )
    return logging.getLogger("l_rclss")


log = logging.getLogger("l_rclss")


def load_config(path: str | Path | None = None) -> dict:
    p = Path(path) if path else CONFIG_PATH
    with open(p, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def resolve_input(rel_or_abs: str | Path) -> Path:
    """输入路径相对 PROJECT_ROOT 解析（绝对路径原样返回）。"""
    p = Path(rel_or_abs)
    return p if p.is_absolute() else (PROJECT_ROOT / p)


def _is_within(path: Path, parent: Path) -> bool:
    try:
        path.resolve().relative_to(parent.resolve())
        return True
    except ValueError:
        return False


def safe_output_path(rel_or_abs: str | Path) -> Path:
    """强制写出落在 OUTPUTS 内（安全契约）。相对路径相对 OUTPUTS 解析。"""
    p = Path(rel_or_abs)
    p = p if p.is_absolute() else (OUTPUTS / p)
    if not _is_within(p, OUTPUTS):
        raise ValueError(f"Refusing to write outside kernel outputs: {p} (root={OUTPUTS})")
    p.parent.mkdir(parents=True, exist_ok=True)
    return p


def read_table(path: str | Path) -> pd.DataFrame:
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(path)
    sfx = path.suffix.lower()
    if sfx == ".parquet":
        return pd.read_parquet(path)
    if sfx in {".csv", ".txt"}:
        return pd.read_csv(path, low_memory=False)
    if sfx == ".tsv":
        return pd.read_csv(path, sep="\t", low_memory=False)
    raise ValueError(f"Unsupported table format: {path}")


def write_table(df: pd.DataFrame, rel_or_abs: str | Path, force: bool = False) -> Path:
    path = safe_output_path(rel_or_abs)
    if path.exists() and not force:
        raise FileExistsError(f"{path} exists. Re-run with --force.")
    sfx = path.suffix.lower()
    if sfx == ".parquet":
        df.to_parquet(path, index=False)
    elif sfx == ".csv":
        df.to_csv(path, index=False)
    else:
        raise ValueError(f"Unsupported output format: {path}")
    log.info("Wrote %s rows -> %s", len(df), path)
    return path


# ============================================================
# 化学小工具（离线；RDKit 仅在这里）
# ============================================================

def normalize_inchikey(value: Any) -> str:
    return "" if pd.isna(value) else str(value).strip().upper()


def inchikey_block1(value: Any) -> str:
    t = normalize_inchikey(value)
    return t.split("-")[0][:14] if t else ""


_MOL_CACHE: dict[str, Any] = {}


def _mol(smiles: str):
    from rdkit import Chem  # 局部 import：保证打分路径文件不会顶层依赖 rdkit
    if smiles in _MOL_CACHE:
        return _MOL_CACHE[smiles]
    m = Chem.MolFromSmiles(smiles) if isinstance(smiles, str) and smiles else None
    _MOL_CACHE[smiles] = m
    return m


def canonical_smiles(smiles: Any) -> str | None:
    if not isinstance(smiles, str) or not smiles.strip():
        return None
    from rdkit import Chem
    m = _mol(smiles.strip())
    if m is None:
        return None
    return Chem.MolToSmiles(m)


def smiles_to_inchikey(smiles: Any) -> str | None:
    m = _mol(smiles.strip()) if isinstance(smiles, str) and smiles.strip() else None
    if m is None:
        return None
    from rdkit import Chem
    try:
        return Chem.MolToInchiKey(m)
    except Exception:
        return None


def murcko_scaffold(smiles: Any) -> str | None:
    """Bemis-Murcko scaffold SMILES（无环骨架退化为空时回退原分子 canonical）。"""
    m = _mol(smiles.strip()) if isinstance(smiles, str) and smiles.strip() else None
    if m is None:
        return None
    from rdkit import Chem
    from rdkit.Chem.Scaffolds import MurckoScaffold
    try:
        scaf = MurckoScaffold.GetScaffoldForMol(m)
        s = Chem.MolToSmiles(scaf)
        return s if s else Chem.MolToSmiles(m)
    except Exception:
        return None


def coalesce_columns(df: pd.DataFrame, candidates: Iterable[str]) -> pd.Series | None:
    for c in candidates:
        if c in df.columns:
            return df[c]
    return None
