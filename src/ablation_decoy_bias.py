"""
ablation_decoy_bias.py — decoy-bias 消融

问题：ML 是真学到"S→P 这条变换合理"，还是只学到"P 看起来像不像真代谢物"（decoy bias）？
做法：同一 scaffold-CV 流水线，只换反应特征 block：
  full   = [zs,zp,diff,prod]   完整
  zp     = [zp]                只看产物 → 若 ≈ full，说明主要靠产物像不像真代谢物（bias）
  zs     = [zs]                只看底物 → 应接近随机（无产物信息）
  diff   = [zs-zp]             只看变换方向
  zs_zp  = [zs,zp]             底物+产物，无显式交互
判读：full 显著 > zp ⇒ 模型用到了底物-产物关系（变换），bias 有限；
      zp ≈ full        ⇒ 模型主要靠产物可信度，decoy bias 强。

用法：python scripts/ssrf/ml_ranking_kernel/src/ablation_decoy_bias.py [--force]
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import kio
import rxnfeat
from train_kernel import run_split
from evaluate_kernel import ranking_metrics


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--modes", nargs="+", default=["full", "zp", "zs", "diff", "zs_zp"])
    ap.add_argument("--splits", nargs="+", default=["scaffold", "random"])
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()
    kio.setup_logging()
    cfg = kio.load_config()
    tr_cfg = dict(cfg["train"]); tr_cfg["l2_normalize"] = cfg["reaction_features"].get("l2_normalize", True)
    seed = int(cfg["splits"]["random_seed"]); topk = cfg["evaluate"]["topk"]

    edges = pd.read_parquet(kio.SPLITS_DIR / "edges_with_splits.parquet")
    cand = pd.read_parquet(kio.DATASETS_DIR / "substrate_candidates.parquet").rename(
        columns={"rclss_production_score": "rclss_prior"})
    cand["substrate_block1"] = cand["substrate_block1"].astype(str).str.upper()
    cand["product_block1"] = cand["product_block1"].astype(str).str.upper()
    emb, dim = rxnfeat.load_embeddings()

    split_cols = {"scaffold": "scaffold_fold", "random": "random_fold"}
    rows = []
    for split in args.splits:
        for mode in args.modes:
            oof = run_split(split_cols[split], split, edges, cand, emb, dim, tr_cfg, seed, mode=mode)
            m = ranking_metrics(oof, "ml_kernel_score", topk)
            rows.append({"split": split, "mode": mode, "n_feat_dims": rxnfeat.mode_ndim(mode, dim),
                         "mrr": m["mrr"], "recall@5": m["recall@5"], "recall@10": m["recall@10"],
                         "top1": m["top1_precision"]})
    res = pd.DataFrame(rows)
    kio.write_table(res, kio.METRICS_DIR / "decoy_bias_ablation.csv", force=args.force)
    kio.log.info("== decoy-bias ablation ==\n%s", res.round(3).to_string(index=False))
    # 自动判读（scaffold）
    sc = res[res.split == "scaffold"].set_index("mode")
    if "full" in sc.index and "zp" in sc.index:
        full_r10, zp_r10 = sc.loc["full", "recall@10"], sc.loc["zp", "recall@10"]
        ratio = zp_r10 / full_r10 if full_r10 else float("nan")
        kio.log.info("[verdict] zp-only / full recall@10 = %.2f  → %s",
                     ratio, "decoy-bias 强(产物主导)" if ratio > 0.85 else
                            "模型用到底物-产物变换关系(bias 有限)" if ratio < 0.7 else "中等")


if __name__ == "__main__":
    main()
