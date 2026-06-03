"""
build_kernel_dataset.py — Phase 1

把多源独立正样本 + 三路负样本 + PU 未标注池，组装成一张统一边表
outputs/datasets/kernel_edges.parquet。

边 = (substrate → product)。列：
  edge_uid, substrate_smiles, substrate_inchikey, substrate_block1, substrate_scaffold,
  product_smiles, product_inchikey, product_block1,
  y (1 正 / 0 负 / -1 未标注), label_source, leakage_flag, confidence, weight,
  reaction_category, multistep, in_pool, rclss_prior, substrate_fanout, primary_pmid

不消费 rule_id/EC/RCLSS 作 ML 特征（rclss_prior 仅留作评估对照基线）。

用法：
  python scripts/ssrf/ml_ranking_kernel/src/build_kernel_dataset.py [--smoke] [--force]
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import kio


# ============================================================
# 归一化小工具
# ============================================================

def map_confidence(value) -> float:
    s = ("" if pd.isna(value) else str(value)).strip().lower()
    if not s:
        return 0.6
    if "medium-high" in s or "med-high" in s:
        return 0.85
    if "high" in s:
        return 1.0
    if "medium" in s or "med" in s:
        return 0.7
    if "low" in s:
        return 0.4
    # curated evidence-type 映射
    if "direct_enzyme_assay" in s:
        return 1.0
    if "strain_fermentation" in s:
        return 0.85
    if "comparative_genomics" in s:
        return 0.6
    return 0.6


def map_leakage(value) -> str:
    s = ("" if pd.isna(value) else str(value)).strip().lower()
    if s.startswith("high"):
        return "high"
    if "medium-high" in s or "med-high" in s:
        return "medium_high"
    if s.startswith("medium") or s.startswith("med"):
        return "medium"
    if s.startswith("low"):
        return "low"
    return s or "unknown"


def is_smiles_usable(value) -> bool:
    s = ("" if pd.isna(value) else str(value)).strip()
    if not s:
        return False
    if "tbd" in s.lower() or "?" == s:
        return False
    return True


# ============================================================
# 加载正样本
# ============================================================

def load_positive_source(src: dict) -> pd.DataFrame:
    path = kio.resolve_input(src["path"])
    df = kio.read_table(path)
    out = pd.DataFrame()
    out["substrate_smiles_raw"] = df[src["substrate_smiles_col"]]
    out["product_smiles_raw"] = df[src["product_smiles_col"]]
    cat_col = src.get("category_col")
    out["reaction_category"] = df[cat_col].astype(str) if cat_col and cat_col in df.columns else ""
    conf_col = src.get("confidence_col")
    out["confidence_raw"] = df[conf_col] if conf_col and conf_col in df.columns else ""
    pmid_col = src.get("pmid_col")
    out["primary_pmid"] = df[pmid_col].astype(str) if pmid_col and pmid_col in df.columns else ""
    out["label_source"] = src["name"]

    # leakage
    leak_col = src.get("leakage_col")
    if src.get("leakage_tier") == "from_column" and leak_col and leak_col in df.columns:
        out["leakage_raw"] = df[leak_col]
    else:
        out["leakage_raw"] = src.get("leakage_tier", "low")

    # review flag → needs_review
    review_col = src.get("review_flag_col")
    if review_col and review_col in df.columns:
        rf = df[review_col].astype(str).str.upper()
        out["needs_review"] = rf.str.contains("REVIEW") | rf.str.contains("NEEDS")
    else:
        out["needs_review"] = False

    # benchmark：单步适用性 → multistep / 排除 No*
    if "is_suitable_for_ssrf_rclss" in df.columns:
        suit = df["is_suitable_for_ssrf_rclss"].astype(str).str.lower()
        out["multistep"] = suit.str.startswith("partial")
        out["exclude_nonstep"] = suit.str.startswith("no")
    else:
        out["multistep"] = False
        out["exclude_nonstep"] = False
    return out


def load_all_positives(cfg: dict) -> tuple[pd.DataFrame, pd.DataFrame]:
    frames = [load_positive_source(s) for s in cfg["positive_sources"]]
    raw = pd.concat(frames, ignore_index=True)

    # canonical 化（RDKit 离线）
    raw["substrate_smiles"] = raw["substrate_smiles_raw"].map(
        lambda s: kio.canonical_smiles(s) if is_smiles_usable(s) else None)
    raw["product_smiles"] = raw["product_smiles_raw"].map(
        lambda s: kio.canonical_smiles(s) if is_smiles_usable(s) else None)

    bad = raw["substrate_smiles"].isna() | raw["product_smiles"].isna()
    review = raw["needs_review"] | raw["exclude_nonstep"]
    drop = bad | review
    needs_review = raw[drop].copy()
    needs_review["drop_reason"] = np.select(
        [bad, raw["needs_review"], raw["exclude_nonstep"]],
        ["smiles_unparseable_or_TBD", "curator_needs_review", "not_single_step_endpoint"],
        default="other",
    )[drop.values] if drop.any() else []

    pos = raw[~drop].copy()
    pos["substrate_inchikey"] = pos["substrate_smiles"].map(kio.smiles_to_inchikey)
    pos["product_inchikey"] = pos["product_smiles"].map(kio.smiles_to_inchikey)
    pos["substrate_block1"] = pos["substrate_inchikey"].map(kio.inchikey_block1)
    pos["product_block1"] = pos["product_inchikey"].map(kio.inchikey_block1)
    pos["substrate_scaffold"] = pos["substrate_smiles"].map(kio.murcko_scaffold)
    pos = pos[pos["substrate_block1"].ne("") & pos["product_block1"].ne("")].copy()

    # 排除退化边（substrate==product）
    pos = pos[pos["substrate_block1"].ne(pos["product_block1"])].copy()

    pos["confidence"] = pos["confidence_raw"].map(map_confidence)
    pos["leakage_flag"] = pos["leakage_raw"].map(map_leakage)
    pos["y"] = 1

    # 去重 by (substrate_block1, product_block1)：保留 confidence 最高、leakage 最低
    leak_rank = {"low": 0, "medium": 1, "medium_high": 2, "high": 3, "unknown": 1}
    pos["_leak_rank"] = pos["leakage_flag"].map(leak_rank).fillna(1)
    pos = pos.sort_values(["confidence", "_leak_rank"], ascending=[False, True])
    pos = pos.drop_duplicates(["substrate_block1", "product_block1"], keep="first").reset_index(drop=True)
    return pos, needs_review


# ============================================================
# 负样本 + PU
# ============================================================

def load_candidate_frame(cfg: dict) -> tuple[pd.DataFrame, str]:
    """优先用 generate_candidates 的 substrate_candidates.parquet（in-substrate，含 rclss 基线 +
    覆盖正样本底物）；否则回退到静态 V4 pool（bootstrap/smoke）。
    返回统一列：substrate_smiles/inchikey/block1, product_smiles/inchikey/block1, rclss_prior, same_product。
    """
    cand_p = kio.DATASETS_DIR / "substrate_candidates.parquet"
    if cand_p.exists():
        df = pd.read_parquet(cand_p)
        df = df.rename(columns={"rclss_production_score": "rclss_prior"})
        df["substrate_block1"] = df["substrate_block1"].astype(str).str.upper()
        df["product_block1"] = df["product_block1"].astype(str).str.upper()
        df["same_product"] = df["substrate_block1"].eq(df["product_block1"])
        keep = ["substrate_smiles", "substrate_inchikey", "substrate_block1",
                "product_smiles", "product_inchikey", "product_block1", "rclss_prior", "same_product"]
        return df[keep], "substrate_candidates"
    return load_pool(cfg), "v4_pool_fallback"


def load_pool(cfg: dict) -> pd.DataFrame:
    cp = cfg["candidate_pool"]
    df = kio.read_table(kio.resolve_input(cp["path"]))
    out = pd.DataFrame()
    out["substrate_smiles"] = df[cp["substrate_smiles_col"]]
    out["product_smiles"] = df[cp["product_smiles_col"]]
    out["substrate_inchikey"] = df[cp["substrate_inchikey_col"]].map(kio.normalize_inchikey)
    out["product_inchikey"] = df[cp["product_inchikey_col"]].map(kio.normalize_inchikey)
    out["substrate_block1"] = (df["substrate_block1"] if "substrate_block1" in df.columns
                               else out["substrate_inchikey"].map(kio.inchikey_block1)).astype(str).str.upper()
    out["product_block1"] = (df["predicted_product_block1"] if "predicted_product_block1" in df.columns
                             else out["product_inchikey"].map(kio.inchikey_block1)).astype(str).str.upper()
    out["rclss_prior"] = pd.to_numeric(df.get(cp["rclss_prior_col"]), errors="coerce")
    spf = cp.get("same_product_flag_col")
    if spf and spf in df.columns:
        out["same_product"] = df[spf].fillna(False).astype(str).str.lower().isin({"true", "1", "yes"})
    else:
        out["same_product"] = out["substrate_block1"].eq(out["product_block1"])
    return out


def build_decoys_and_pu(pos: pd.DataFrame, pool: pd.DataFrame, cfg: dict, rng: np.random.Generator):
    neg_cfg = cfg["negatives"]["in_substrate_decoy"]
    max_per = int(neg_cfg["max_per_substrate"])

    # substrate fanout（每底物候选数）
    fanout = pool.groupby("substrate_block1").size().rename("substrate_fanout")

    pos_subs = set(pos["substrate_block1"])
    # 每底物正产物集合（排除这些 product_block1 作为 decoy，避免假负）
    pos_prod_by_sub = pos.groupby("substrate_block1")["product_block1"].agg(set).to_dict()

    cand = pool[pool["substrate_block1"].isin(pos_subs)].copy()
    if neg_cfg.get("exclude_same_product", True):
        cand = cand[~cand["same_product"]]
    cand = cand[cand["substrate_block1"].ne(cand["product_block1"])]

    decoy_rows = []
    for sub, grp in cand.groupby("substrate_block1"):
        forbidden = pos_prod_by_sub.get(sub, set())
        g = grp[~grp["product_block1"].isin(forbidden)]
        g = g.drop_duplicates("product_block1")
        if len(g) > max_per:
            g = g.sample(n=max_per, random_state=int(rng.integers(0, 2**31 - 1)))
        decoy_rows.append(g)
    decoys = pd.concat(decoy_rows, ignore_index=True) if decoy_rows else pool.iloc[0:0].copy()
    decoys["y"] = 0
    decoys["label_source"] = "in_substrate_decoy"
    decoys["leakage_flag"] = "decoy"
    decoys["confidence"] = 1.0
    decoys["multistep"] = False
    decoys["primary_pmid"] = ""
    decoys["reaction_category"] = ""

    # PU 未标注池：pool 里既非正边也非 decoy 的边，随机采样
    labeled_keys = set(zip(pos["substrate_block1"], pos["product_block1"])) | \
                   set(zip(decoys["substrate_block1"], decoys["product_block1"]))
    pool_keys = list(zip(pool["substrate_block1"], pool["product_block1"]))
    mask = np.array([k not in labeled_keys for k in pool_keys])
    unl = pool[mask & ~pool["same_product"]].copy()
    n_pu = int(cfg["negatives"]["pu_unlabeled"]["sample_size"])
    if len(unl) > n_pu:
        unl = unl.sample(n=n_pu, random_state=int(rng.integers(0, 2**31 - 1)))
    unl["y"] = -1
    unl["label_source"] = "pu_unlabeled"
    unl["leakage_flag"] = "unlabeled"
    unl["confidence"] = 0.0
    unl["multistep"] = False
    unl["primary_pmid"] = ""
    unl["reaction_category"] = ""

    return decoys, unl, fanout


# ============================================================
# 主流程
# ============================================================

def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default=None)
    ap.add_argument("--smoke", action="store_true", help="小样本快跑（PU 池缩到 2000）")
    ap.add_argument("--with-network-negatives", action="store_true",
                    help="并入 build_network_negatives 的 AGREDA 隐式负（Phase 1b）")
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--verbose", action="store_true")
    args = ap.parse_args()

    kio.setup_logging(args.verbose)
    cfg = kio.load_config(args.config)
    rng = np.random.default_rng(int(cfg["splits"]["random_seed"]))
    if args.smoke:
        cfg["negatives"]["pu_unlabeled"]["sample_size"] = 2000

    kio.log.info("== Phase 1: build kernel dataset ==")
    pos, needs_review = load_all_positives(cfg)
    kio.log.info("positives (deduped): %d   needs_review/excluded: %d", len(pos), len(needs_review))
    kio.log.info("  by source: %s", pos["label_source"].value_counts().to_dict())
    kio.log.info("  leakage:   %s", pos["leakage_flag"].value_counts().to_dict())

    pool, pool_src = load_candidate_frame(cfg)
    kio.log.info("candidate source = %s : %d edges, %d substrates (cover pos=%d/%d)",
                 pool_src, len(pool), pool["substrate_block1"].nunique(),
                 len(set(pool["substrate_block1"]) & set(pos["substrate_block1"])),
                 pos["substrate_block1"].nunique())

    decoys, unl, fanout = build_decoys_and_pu(pos, pool, cfg, rng)
    kio.log.info("in-substrate decoys: %d (over %d substrates)",
                 len(decoys), decoys["substrate_block1"].nunique())
    kio.log.info("PU unlabeled: %d", len(unl))

    # 正样本 in_pool 标记 + 取 rclss_prior（评估基线）
    pool_key = pool.drop_duplicates(["substrate_block1", "product_block1"])[
        ["substrate_block1", "product_block1", "rclss_prior"]]
    pos = pos.merge(pool_key, on=["substrate_block1", "product_block1"], how="left")
    pos["in_pool"] = pos["rclss_prior"].notna()
    kio.log.info("positives present in V4 pool (RCLSS-scoreable): %d / %d",
                 int(pos["in_pool"].sum()), len(pos))

    # decoy/unlabeled 也需要 scaffold（仅 decoy；PU 标 train-only 无需 scaffold）
    decoys["substrate_scaffold"] = decoys["substrate_smiles"].map(kio.murcko_scaffold)
    unl["substrate_scaffold"] = ""

    # 统一列
    keep = ["substrate_smiles", "substrate_inchikey", "substrate_block1", "substrate_scaffold",
            "product_smiles", "product_inchikey", "product_block1",
            "y", "label_source", "leakage_flag", "confidence", "multistep",
            "reaction_category", "primary_pmid", "rclss_prior"]
    for df in (pos, decoys, unl):
        for c in keep:
            if c not in df.columns:
                df[c] = "" if df is not pos and c in ("substrate_scaffold",) else np.nan

    pos["in_pool"] = pos["in_pool"]
    decoys["in_pool"] = True
    unl["in_pool"] = True

    frames = [pos[keep + ["in_pool"]], decoys[keep + ["in_pool"]], unl[keep + ["in_pool"]]]

    # Phase 1b: AGREDA 网络隐式负
    netneg_p = kio.DATASETS_DIR / "network_negatives.parquet"
    if args.with_network_negatives and netneg_p.exists():
        nn = pd.read_parquet(netneg_p)
        nn["substrate_scaffold"] = nn["substrate_smiles"].map(kio.murcko_scaffold)
        nn["multistep"] = False
        nn["reaction_category"] = nn.get("reaction_category", "")
        nn["primary_pmid"] = ""
        nn["in_pool"] = True
        for c in keep:
            if c not in nn.columns:
                nn[c] = np.nan
        # 去掉已是正样本/decoy 的边
        existing = set(zip(pd.concat(frames)["substrate_block1"], pd.concat(frames)["product_block1"]))
        nn = nn[~nn.apply(lambda r: (r["substrate_block1"], r["product_block1"]) in existing, axis=1)]
        frames.append(nn[keep + ["in_pool"]])
        kio.log.info("merged %d AGREDA network-implicit negatives", len(nn))

    edges = pd.concat(frames, ignore_index=True)

    # weight：正=confidence*(0.6 if multistep)；decoy=1；unlabeled=0（训练里 PU 处理）
    edges["weight"] = np.where(
        edges["y"] == 1, edges["confidence"] * np.where(edges["multistep"], 0.6, 1.0),
        np.where(edges["y"] == 0, 1.0, 0.0))

    # substrate_fanout
    edges = edges.merge(fanout, on="substrate_block1", how="left")
    edges["substrate_fanout"] = edges["substrate_fanout"].fillna(0).astype(int)

    edges = edges.reset_index(drop=True)
    edges["edge_uid"] = (edges["substrate_block1"] + "__" + edges["product_block1"]
                         + "__" + edges["label_source"] + "__" + edges.index.astype(str))

    # 输出
    kio.write_table(edges, kio.DATASETS_DIR / "kernel_edges.parquet", force=args.force)
    if len(needs_review):
        kio.write_table(needs_review.drop(columns=[c for c in ["multistep", "exclude_nonstep"]
                                                   if c in needs_review.columns], errors="ignore"),
                        kio.DATASETS_DIR / "needs_review_queue.csv", force=args.force)

    # 摘要
    summary = (edges.assign(y_name=edges["y"].map({1: "positive", 0: "negative", -1: "unlabeled"}))
               .groupby(["y_name", "label_source"]).size().rename("n").reset_index())
    kio.write_table(summary, kio.DATASETS_DIR / "kernel_edges_summary.csv", force=args.force)
    kio.log.info("== summary ==\n%s", summary.to_string(index=False))
    kio.log.info("total edges: %d  (pos=%d neg=%d unl=%d)",
                 len(edges), int((edges.y == 1).sum()), int((edges.y == 0).sum()), int((edges.y == -1).sum()))


if __name__ == "__main__":
    main()
