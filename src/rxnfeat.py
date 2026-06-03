"""
rxnfeat.py — 反应特征装配（共享给 train/evaluate/predict）

反应表示 r = [z_s, z_p, z_s - z_p, z_s ⊙ z_p]，z 来自冻结化学 LM（mol_embeddings.parquet）。
**仅** 分子 embedding 衍生 → 规则无关（不含 rule_id/EC/RCLSS/category）。
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import kio


def load_embeddings(suffix: str = "") -> tuple[dict[str, np.ndarray], int]:
    p = kio.FEATURES_DIR / f"mol_embeddings{suffix}.parquet"
    df = pd.read_parquet(p)
    emb = {s: np.asarray(v, dtype=np.float32) for s, v in zip(df["smiles"], df["embedding"])}
    dim = len(next(iter(emb.values()))) if emb else 0
    return emb, dim


def _l2(v: np.ndarray) -> np.ndarray:
    n = np.linalg.norm(v, axis=-1, keepdims=True)
    return v / np.clip(n, 1e-9, None)


# 特征 mode → 用哪些 block（decoy-bias ablation 用）
#   full   = [zs, zp, diff, prod]   完整反应特征
#   zp     = [zp]                   纯产物（测"只学产物像不像真代谢物"=decoy bias）
#   zs     = [zs]                   纯底物
#   diff   = [diff]                 纯变换方向 z_s - z_p
#   zs_zp  = [zs, zp]               底物+产物，无显式交互
MODE_BLOCKS = {
    "full": ("zs", "zp", "diff", "prod"),
    "zp": ("zp",),
    "zs": ("zs",),
    "diff": ("diff",),
    "zs_zp": ("zs", "zp"),
}


def mode_ndim(mode: str, dim: int) -> int:
    return len(MODE_BLOCKS[mode]) * dim


def build_features(
    sub_smiles,
    prod_smiles,
    emb: dict[str, np.ndarray],
    dim: int,
    l2_normalize: bool = True,
    mode: str = "full",
) -> tuple[np.ndarray, np.ndarray]:
    """返回 (X[n, ndim], ok_mask[n])。两端 embedding 缺失的行 ok=False。"""
    blocks = MODE_BLOCKS[mode]
    sub_smiles = list(sub_smiles)
    prod_smiles = list(prod_smiles)
    n = len(sub_smiles)
    X = np.zeros((n, len(blocks) * dim), dtype=np.float32)
    ok = np.zeros(n, dtype=bool)
    for i, (s, p) in enumerate(zip(sub_smiles, prod_smiles)):
        zs = emb.get(s)
        zp = emb.get(p)
        if zs is None or zp is None:
            continue
        if l2_normalize:
            zs = _l2(zs)
            zp = _l2(zp)
        parts = {"zs": zs, "zp": zp, "diff": zs - zp, "prod": zs * zp}
        X[i] = np.concatenate([parts[b] for b in blocks])
        ok[i] = True
    return X, ok


def feature_block_names(dim: int, mode: str = "full") -> list[str]:
    """特征列名（用于规则无关断言 + provenance）。全部 emb 衍生，无规则字段。"""
    names = []
    for blk in MODE_BLOCKS[mode]:
        names += [f"{blk}_{i}" for i in range(dim)]
    return names
