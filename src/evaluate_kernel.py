"""
evaluate_kernel.py — Phase 6

读 train_kernel 的 OOF 打分候选全集，算 {rclss, ml, hybrid_gate, hybrid_fuse} × {scaffold, random}
的 within-substrate 排序指标（TopK recall / MRR / top1-precision）+ AUROC/AUPRC，
并判定成功门（标准式）：
  ML(scaffold) > RCLSS(scaffold)（TopK & MRR） 且 |ML(scaffold)-ML(random)| ≤ tol。

用法：
  python scripts/ssrf/ml_ranking_kernel/src/evaluate_kernel.py [--alpha 0.5] [--force]
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import kio
import fuse_hybrid


def competition_ranks(scores: np.ndarray) -> np.ndarray:
    """**悲观** tie-breaking：并列项取该组最差名次 rank = #{score_j >= score_i}。
    这样"所有候选打同分"不会让真产物白白得 rank 1（避免 zs-only 假满分这类伪信号）。
    NaN/-inf 自然排末。"""
    s = np.where(np.isnan(scores), -np.inf, scores)
    order = np.argsort(-s, kind="mergesort")
    ranks = np.empty(len(s), dtype=float)
    sorted_s = s[order]
    ranks_sorted = np.zeros(len(s))
    i = 0
    while i < len(sorted_s):
        j = i
        while j < len(sorted_s) and sorted_s[j] == sorted_s[i]:
            j += 1
        for k in range(i, j):
            ranks_sorted[k] = j   # 悲观：该并列组的最末 1-indexed 名次
        i = j
    ranks[order] = ranks_sorted
    return ranks


def ranking_metrics(df: pd.DataFrame, score_col: str, topk: list[int]) -> dict:
    """per (substrate, true-product) case 的 rank → TopK recall / MRR / top1-precision."""
    case_ranks = []
    top1_true = 0
    n_sub = 0
    for sub, g in df.groupby("substrate_block1"):
        ranks = competition_ranks(g[score_col].to_numpy(dtype=float))
        y = g["y_eval"].to_numpy()
        n_sub += 1
        if (ranks == 1).any():
            top1_true += int(y[ranks == 1].max() == 1)
        for r, yy in zip(ranks, y):
            if yy == 1:
                case_ranks.append(r)
    case_ranks = np.array(case_ranks, dtype=float)
    out = {"n_cases": int(len(case_ranks)), "n_substrates": int(n_sub),
           "mrr": float(np.mean(1.0 / case_ranks)) if len(case_ranks) else 0.0,
           "top1_precision": float(top1_true / n_sub) if n_sub else 0.0}
    for k in topk:
        out[f"recall@{k}"] = float(np.mean(case_ranks <= k)) if len(case_ranks) else 0.0
    return out


def auc_metrics(df: pd.DataFrame, score_col: str) -> dict:
    from sklearn.metrics import roc_auc_score, average_precision_score
    s = pd.to_numeric(df[score_col], errors="coerce").replace([np.inf, -np.inf], np.nan)
    sub = df.assign(_s=s).dropna(subset=["_s"])
    y = sub["y_eval"].to_numpy()
    if len(np.unique(y)) < 2:
        return {"auroc": None, "auprc": None, "n": int(len(sub))}
    return {"auroc": float(roc_auc_score(y, sub["_s"])),
            "auprc": float(average_precision_score(y, sub["_s"])),
            "n": int(len(sub))}


def evaluate_file(path: Path, split_name: str, topk: list[int], alpha: float,
                  clean_only: bool = False) -> pd.DataFrame:
    df = pd.read_parquet(path)
    if clean_only:
        # 排除 High-leakage 正样本案例（仅评估这些底物的干净度）
        bad_subs = set(df.loc[(df.y_eval == 1) & (df["leakage_flag"].astype(str).str.lower() == "high"),
                              "substrate_block1"])
        df = df[~df["substrate_block1"].isin(bad_subs)]
    df = fuse_hybrid.add_fusion_scores(df, alpha=alpha)
    rows = []
    for system, col in fuse_hybrid.SYSTEM_SCORE_COL.items():
        rk = ranking_metrics(df, col, topk)
        au = auc_metrics(df, col)
        rows.append({"split": split_name, "system": system, "clean_only": clean_only,
                     **rk, **{f"auc_{k}": v for k, v in au.items()}})
    return pd.DataFrame(rows)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--alpha", type=float, default=0.5)
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()
    kio.setup_logging()
    cfg = kio.load_config()
    topk = cfg["evaluate"]["topk"]
    gap_tol = float(cfg["evaluate"].get("scaffold_random_gap_tol", 0.1))

    frames = []
    for split in ["scaffold", "random"]:
        p = kio.PREDICTIONS_DIR / f"oof_scored_edges.{split}.parquet"
        if not p.exists():
            kio.log.warning("missing %s", p)
            continue
        frames.append(evaluate_file(p, split, topk, args.alpha, clean_only=False))
        frames.append(evaluate_file(p, split, topk, args.alpha, clean_only=True))
    if not frames:
        raise SystemExit("No OOF predictions. Run train_kernel.py first.")
    metrics = pd.concat(frames, ignore_index=True)
    kio.write_table(metrics, kio.METRICS_DIR / "ranking_metrics.csv", force=args.force)

    # ===== 成功门（标准式）=====
    def get(split, system, metric, clean=False):
        m = metrics[(metrics.split == split) & (metrics.system == system) & (metrics.clean_only == clean)]
        return float(m[metric].iloc[0]) if len(m) and pd.notna(m[metric].iloc[0]) else None

    gate = {}
    for metric in ["recall@5", "recall@10", "mrr", "top1_precision"]:
        ml_s = get("scaffold", "ml", metric)
        rc_s = get("scaffold", "rclss", metric)
        ml_r = get("random", "ml", metric)
        gate[metric] = {
            "ml_scaffold": ml_s, "rclss_scaffold": rc_s, "ml_random": ml_r,
            "ml_beats_rclss_scaffold": (ml_s is not None and rc_s is not None and ml_s > rc_s),
            "scaffold_random_gap": (abs(ml_s - ml_r) if ml_s is not None and ml_r is not None else None),
            "gap_small": (ml_s is not None and ml_r is not None and abs(ml_s - ml_r) <= gap_tol),
        }
    beats = all(gate[m]["ml_beats_rclss_scaffold"] for m in ["recall@5", "recall@10", "mrr"])
    small_gap = all(gate[m]["gap_small"] for m in ["recall@5", "recall@10", "mrr"]
                    if gate[m]["scaffold_random_gap"] is not None)
    verdict = {
        "criterion": "standard: ML(scaffold) > RCLSS(scaffold) AND |ML(scaffold)-ML(random)| small",
        "gap_tol": gap_tol,
        "ml_beats_rclss_on_scaffold": beats,
        "small_scaffold_random_gap": small_gap,
        "PASS": bool(beats and small_gap),
        "per_metric": gate,
    }
    kio.safe_output_path("reports/success_gate.json").write_text(
        json.dumps(verdict, indent=2, ensure_ascii=False), encoding="utf-8")

    # 控制台摘要
    show = metrics[~metrics.clean_only][["split", "system", "n_cases", "mrr",
                                         "recall@5", "recall@10", "top1_precision", "auc_auprc"]]
    kio.log.info("== ranking metrics (all cases) ==\n%s", show.round(3).to_string(index=False))
    kio.log.info("== SUCCESS GATE: PASS=%s (beats_rclss=%s, small_gap=%s) ==",
                 verdict["PASS"], beats, small_gap)
    for m in ["recall@5", "recall@10", "mrr"]:
        g = gate[m]
        kio.log.info("  %-12s ml_scaf=%.3f rclss_scaf=%.3f ml_rand=%.3f gap=%.3f",
                     m, g["ml_scaffold"] or 0, g["rclss_scaffold"] or 0, g["ml_random"] or 0,
                     g["scaffold_random_gap"] or 0)


if __name__ == "__main__":
    main()
