"""
make_splits.py — Phase 3

按 Bemis-Murcko scaffold（底物）做泄漏控制划分：
  - scaffold split：GroupKFold by substrate_scaffold（难，跨 scaffold 泛化）
  - random split  ：KFold by substrate_block1（易，同 scaffold 可跨 train/test 泄漏 = 乐观上界）
两套都在**底物粒度**划分，保证某底物的正样本+in-substrate decoy 同进同出。
PU 未标注（y=-1）恒为 train-only（fold=-1）。

泄漏隔离：leakage_flag==High 的正样本标 quarantine=True（评估时可排除出干净测试）。

断言：scaffold split 下任一 fold 的 scaffold 集合与其它 fold 不相交。

用法：
  python scripts/ssrf/ml_ranking_kernel/src/make_splits.py [--force]
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.model_selection import GroupKFold, KFold

sys.path.insert(0, str(Path(__file__).resolve().parent))
import kio


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()
    kio.setup_logging()
    cfg = kio.load_config()
    n_folds = int(cfg["splits"]["n_folds"])
    seed = int(cfg["splits"]["random_seed"])
    quarantine_level = str(cfg["splits"].get("quarantine_leakage", "High")).lower()

    edges = pd.read_parquet(kio.DATASETS_DIR / "kernel_edges.parquet")
    labeled = edges[edges["y"].isin([0, 1])].copy()

    # 每底物一行：scaffold（取正样本的 scaffold；decoy 同底物共享）
    sub = (labeled.sort_values("y", ascending=False)
           .drop_duplicates("substrate_block1")[["substrate_block1", "substrate_scaffold"]]
           .reset_index(drop=True))
    sub["substrate_scaffold"] = sub["substrate_scaffold"].fillna("").replace("", np.nan)
    # scaffold 缺失的，用自身 block1 当作独立 scaffold（不与他人合并）
    sub["scaffold_key"] = sub["substrate_scaffold"].fillna("NOSCAF_" + sub["substrate_block1"])

    n_sub = len(sub)
    k = min(n_folds, n_sub, sub["scaffold_key"].nunique())
    kio.log.info("substrates=%d unique_scaffolds=%d folds=%d", n_sub, sub["scaffold_key"].nunique(), k)

    # scaffold split：GroupKFold by scaffold
    gkf = GroupKFold(n_splits=k)
    sub["scaffold_fold"] = -1
    for f, (_, te) in enumerate(gkf.split(sub, groups=sub["scaffold_key"])):
        sub.loc[sub.index[te], "scaffold_fold"] = f

    # random split：KFold by substrate（打散 scaffold）
    kf = KFold(n_splits=k, shuffle=True, random_state=seed)
    sub["random_fold"] = -1
    order = sub.sample(frac=1.0, random_state=seed).index.to_numpy()
    for f, (_, te) in enumerate(kf.split(order)):
        sub.loc[order[te], "random_fold"] = f

    # 断言：scaffold split 各 fold scaffold 不相交
    scaf_by_fold = sub.groupby("scaffold_fold")["scaffold_key"].agg(set)
    for a in scaf_by_fold.index:
        for b in scaf_by_fold.index:
            if a < b:
                inter = scaf_by_fold[a] & scaf_by_fold[b]
                assert not inter, f"scaffold leak between fold {a},{b}: {list(inter)[:3]}"
    kio.log.info("[assert] scaffold folds disjoint: OK")

    # quarantine：High-leakage 正样本的底物
    pos = labeled[labeled["y"] == 1]
    quar_subs = set(pos.loc[pos["leakage_flag"].str.lower().eq(quarantine_level), "substrate_block1"])
    sub["quarantine"] = sub["substrate_block1"].isin(quar_subs)
    kio.log.info("quarantined (leakage=%s) substrates: %d", quarantine_level, int(sub["quarantine"].sum()))

    kio.write_table(sub, kio.SPLITS_DIR / "substrate_splits.parquet", force=args.force)

    # 同时落一个 edge-level 合并版（方便 train/eval 直接用）
    merged = edges.merge(sub[["substrate_block1", "scaffold_key", "scaffold_fold",
                              "random_fold", "quarantine"]],
                         on="substrate_block1", how="left")
    merged["scaffold_fold"] = merged["scaffold_fold"].fillna(-1).astype(int)   # PU/未覆盖 → train-only
    merged["random_fold"] = merged["random_fold"].fillna(-1).astype(int)
    merged["quarantine"] = merged["quarantine"].fillna(False)
    kio.write_table(merged, kio.SPLITS_DIR / "edges_with_splits.parquet", force=args.force)
    kio.log.info("fold sizes (scaffold): %s",
                 sub["scaffold_fold"].value_counts().sort_index().to_dict())


if __name__ == "__main__":
    main()
