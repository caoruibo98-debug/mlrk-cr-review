"""
fuse_hybrid.py — Phase 5

混合融合（用户"SMARTS 稳 + 分子ML 拓宽"）。三个排序变体：
  - ml          : 纯学习内核分（headline，规则无关）
  - hybrid_gate : SMARTS gate × ML（只在规则生成的候选里用 ML 排序；机制有效性闸）
  - hybrid_fuse : 每底物 min-max 归一后 α·RCLSS + (1-α)·ML（late fusion）
对照基线：
  - rclss       : 规则先验分（rclss_prior）
"""
from __future__ import annotations

import numpy as np
import pandas as pd


def _minmax_per_group(df: pd.DataFrame, col: str, group: str) -> pd.Series:
    def f(s):
        v = s.astype(float)
        lo, hi = v.min(), v.max()
        return (v - lo) / (hi - lo) if hi > lo else v * 0.0
    return df.groupby(group)[col].transform(f)


def add_fusion_scores(df: pd.DataFrame, alpha: float = 0.5, group: str = "substrate_block1") -> pd.DataFrame:
    """加 hybrid_gate / hybrid_fuse 分列。要求 df 含 ml_kernel_score, rclss_prior。"""
    out = df.copy()
    out["rclss_prior"] = pd.to_numeric(out["rclss_prior"], errors="coerce")
    has_rule = out["rclss_prior"].notna()

    # gate：规则没生成的候选（rclss 缺）打到 -inf（不参与 top 排序）
    out["score_gate"] = np.where(has_rule, out["ml_kernel_score"], -np.inf)

    # fuse：每底物归一后线性融合；规则缺的候选 rclss 归一记 0
    norm_ml = _minmax_per_group(out, "ml_kernel_score", group)
    out["_rclss_filled"] = out["rclss_prior"].fillna(out.groupby(group)["rclss_prior"].transform("min"))
    norm_rc = _minmax_per_group(out, "_rclss_filled", group).fillna(0.0)
    out["score_fuse"] = alpha * norm_rc + (1 - alpha) * norm_ml
    out = out.drop(columns=["_rclss_filled"])
    return out


SYSTEM_SCORE_COL = {
    "rclss": "rclss_prior",
    "ml": "ml_kernel_score",
    "hybrid_gate": "score_gate",
    "hybrid_fuse": "score_fuse",
}
