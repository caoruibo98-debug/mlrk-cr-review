"""
analyze_report.py — 面向论文的诚实分析汇总

读 OOF 打分 + 两份 metrics（baseline / +netneg），输出：
  1) per-fold CV 方差（mean±std）——诚实标注 n 小、CI 宽
  2) reaction-category 分层（重点 deglycosylation / glycoside hydrolysis）
  3) clean subset（剔除 High-leakage 正样本）
  4) network-negatives 增益 delta
  5) decoy-bias 消融摘要（zp-only vs full）
写到 outputs/reports/analysis_report.md。
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import kio
from evaluate_kernel import ranking_metrics

TOPK = [1, 5, 10, 20]
GLYCO_KW = ("glycos", "deglyc", "galactosid", "glucosid", "glucuron")


def per_fold_ci(oof: pd.DataFrame, score_col: str) -> dict:
    vals = {"recall@10": [], "recall@5": [], "mrr": []}
    for f, g in oof.groupby("fold"):
        m = ranking_metrics(g, score_col, TOPK)
        for k in vals:
            vals[k].append(m[k])
    return {k: (float(np.mean(v)), float(np.std(v))) for k, v in vals.items()}


def subset_metrics(oof: pd.DataFrame, keep_subs: set, score_col: str) -> dict:
    sub = oof[oof["substrate_block1"].isin(keep_subs)]
    return ranking_metrics(sub, score_col, TOPK)


def main() -> None:
    kio.setup_logging()
    oof = pd.read_parquet(kio.PREDICTIONS_DIR / "oof_scored_edges.scaffold.parquet")
    pos = oof[oof.y_eval == 1].copy()
    lines = ["# L-RCLSS 诚实分析汇总 (scaffold split, +network negatives)\n"]

    # 1) per-fold CV 方差
    lines.append("## 1) Per-fold CV 方差 (n=92, ~18 cases/fold → CI 宽)\n")
    for name, col in [("ml", "ml_kernel_score"), ("rclss", "rclss_prior")]:
        ci = per_fold_ci(oof, col)
        lines.append(f"- **{name}**: recall@10 {ci['recall@10'][0]:.3f}±{ci['recall@10'][1]:.3f} | "
                     f"recall@5 {ci['recall@5'][0]:.3f}±{ci['recall@5'][1]:.3f} | "
                     f"mrr {ci['mrr'][0]:.3f}±{ci['mrr'][1]:.3f}")
    lines.append("")

    # 2) reaction-category 分层（glycoside/deglyc vs rest）
    cat = pos["reaction_category"].fillna("").str.lower()
    glyco_subs = set(pos.loc[cat.apply(lambda c: any(k in c for k in GLYCO_KW)), "substrate_block1"])
    rest_subs = set(pos["substrate_block1"]) - glyco_subs
    lines.append(f"## 2) Reaction-category 分层\n")
    lines.append(f"glycoside/deglyc 相关底物: {len(glyco_subs)}；其余: {len(rest_subs)}\n")
    for label, subs in [("glycoside/deglyc", glyco_subs), ("rest", rest_subs)]:
        for name, col in [("ml", "ml_kernel_score"), ("rclss", "rclss_prior")]:
            m = subset_metrics(oof, subs, col)
            lines.append(f"- {label:16s} {name:5s}: recall@10={m['recall@10']:.3f} recall@5={m['recall@5']:.3f} "
                         f"mrr={m['mrr']:.3f} (n={m['n_cases']})")
    lines.append("")

    # 3) clean subset（剔除 High-leakage 正样本底物）
    high_subs = set(pos.loc[pos["leakage_flag"].astype(str).str.lower() == "high", "substrate_block1"])
    clean_subs = set(pos["substrate_block1"]) - high_subs
    lines.append(f"## 3) Clean subset（剔除 {len(high_subs)} 个 High-leakage 底物 → {len(clean_subs)} 底物）\n")
    for name, col in [("ml", "ml_kernel_score"), ("rclss", "rclss_prior")]:
        m = subset_metrics(oof, clean_subs, col)
        lines.append(f"- {name:5s} clean: recall@10={m['recall@10']:.3f} recall@5={m['recall@5']:.3f} "
                     f"mrr={m['mrr']:.3f} (n={m['n_cases']})")
    lines.append("")

    # 4) network-negatives delta
    lines.append("## 4) Network-negatives 增益 (baseline → +netneg, scaffold ml)\n")
    base_p = kio.METRICS_DIR / "ranking_metrics.baseline.csv"
    cur_p = kio.METRICS_DIR / "ranking_metrics.csv"
    if base_p.exists():
        b = pd.read_parquet(base_p) if base_p.suffix == ".parquet" else pd.read_csv(base_p)
        c = pd.read_csv(cur_p)
        def get(df, metric):
            r = df[(df.split == "scaffold") & (df.system == "ml") & (~df.clean_only)]
            return float(r[metric].iloc[0]) if len(r) else float("nan")
        for metric in ["recall@5", "recall@10", "mrr", "auc_auprc"]:
            lines.append(f"- {metric:10s}: {get(b,metric):.3f} → {get(c,metric):.3f}  (Δ {get(c,metric)-get(b,metric):+.3f})")
    else:
        lines.append("- (no baseline metrics saved)")
    lines.append("")

    # 5) decoy-bias 消融摘要
    abl_p = kio.METRICS_DIR / "decoy_bias_ablation.csv"
    if abl_p.exists():
        a = pd.read_csv(abl_p)
        sc = a[a.split == "scaffold"].set_index("mode")
        lines.append("## 5) Decoy-bias 消融 (scaffold)\n")
        for mode in ["full", "zp", "diff", "zs_zp", "zs"]:
            if mode in sc.index:
                lines.append(f"- {mode:6s}: recall@10={sc.loc[mode,'recall@10']:.3f} recall@5={sc.loc[mode,'recall@5']:.3f} "
                             f"mrr={sc.loc[mode,'mrr']:.3f} top1={sc.loc[mode,'top1']:.3f}")
        if "full" in sc.index and "zp" in sc.index:
            lines.append(f"\n  → zp-only/full: recall@10={sc.loc['zp','recall@10']/sc.loc['full','recall@10']:.2f}, "
                         f"recall@5={sc.loc['zp','recall@5']/sc.loc['full','recall@5']:.2f}, "
                         f"top1={sc.loc['zp','top1']/max(sc.loc['full','top1'],1e-9):.2f} "
                         f"(产物主导→decoy bias；底物在 top 名次有小幅增益)")

    report = "\n".join(lines)
    kio.safe_output_path("reports/analysis_report.md").write_text(report, encoding="utf-8")
    print(report)


if __name__ == "__main__":
    main()
