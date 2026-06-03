"""
build_food_test.py — Phase 1（每模块的"唯一诚实测试集"）

聚合 benchmark+gold+curated 的**严格单步**食品正样本，按 module_router 路由到 A/B/C/D，
canonical/去重/标 leakage。这是每模块 recall 报告的唯一来源。

单步判定：
  - benchmark: is_suitable_for_ssrf_rclss == "Yes"（严格单步；Partial/No 排除）
  - gold_edges: 默认单步（手写金标准单反应），reaction_type 含 pathway/net/+ 的排除
  - curated_gap7: curator_flag==OK 且 reaction 不含多步标志

输出 outputs/modular/food_clean.parquet + 每模块计数。
用法：python scripts/ssrf/ml_ranking_kernel/modular/build_food_test.py [--force]
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import kio  # noqa: E402
import module_router as mr  # noqa: E402

MODULAR_OUT = kio.OUTPUTS / "modular"

_MULTISTEP_KW = ("pathway", "net", " + ", "multi", "sequential", "->", "→")


def _usable(s) -> bool:
    s = ("" if pd.isna(s) else str(s)).strip()
    return bool(s) and "tbd" not in s.lower()


def _single_step_text(cat) -> bool:
    c = ("" if pd.isna(cat) else str(cat)).lower()
    return not any(k in c for k in _MULTISTEP_KW)


def load_benchmark() -> pd.DataFrame:
    b = pd.read_csv(kio.resolve_input("analysis/benchmark_ssrf_rclss_2026-05-22.csv"))
    suit = b["is_suitable_for_ssrf_rclss"].astype(str).str.lower()
    b = b[suit.str.startswith("yes")].copy()          # 严格单步
    out = pd.DataFrame({
        "substrate_smiles_raw": b["substrate_smiles"], "product_smiles_raw": b["product_smiles"],
        "reaction_category": b["reaction_type"].astype(str),
        "leakage_raw": b["possible_data_leakage_risk"], "label_source": "benchmark"})
    return out


def load_gold() -> pd.DataFrame:
    g = pd.read_csv(kio.resolve_input("scripts/ssrf/known_reaction_gold_edges.csv"))
    g = g[g["reaction_type"].map(_single_step_text)].copy()
    return pd.DataFrame({
        "substrate_smiles_raw": g["compound_smiles"], "product_smiles_raw": g["known_product_smiles"],
        "reaction_category": g["reaction_type"].astype(str),
        "leakage_raw": "low", "label_source": "gold"})


def load_curated() -> pd.DataFrame:
    c = pd.read_csv(kio.resolve_input(
        "Dataset/curated_literature/2026_04_19_microbiome_reaction_gap_fill/curated_literature_gap7_import_ready_v1.csv"))
    c = c[c["curator_flag"].astype(str).str.upper().eq("OK")].copy()
    c = c[c["reaction_category"].map(_single_step_text)].copy()
    return pd.DataFrame({
        "substrate_smiles_raw": c["substrate_smiles"], "product_smiles_raw": c["product_smiles"],
        "reaction_category": c["reaction_category"].astype(str),
        "leakage_raw": "low", "label_source": "curated"})


def map_leakage(v) -> str:
    s = ("" if pd.isna(v) else str(v)).strip().lower()
    if s.startswith("high"):
        return "high"
    if "medium-high" in s:
        return "medium_high"
    if s.startswith("medium"):
        return "medium"
    return "low"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()
    kio.setup_logging()

    raw = pd.concat([load_benchmark(), load_gold(), load_curated()], ignore_index=True)
    raw = raw[raw["substrate_smiles_raw"].map(_usable) & raw["product_smiles_raw"].map(_usable)].copy()

    raw["substrate_smiles"] = raw["substrate_smiles_raw"].map(kio.canonical_smiles)
    raw["product_smiles"] = raw["product_smiles_raw"].map(kio.canonical_smiles)
    raw = raw.dropna(subset=["substrate_smiles", "product_smiles"]).copy()
    raw["substrate_inchikey"] = raw["substrate_smiles"].map(kio.smiles_to_inchikey)
    raw["product_inchikey"] = raw["product_smiles"].map(kio.smiles_to_inchikey)
    raw["substrate_block1"] = raw["substrate_inchikey"].map(kio.inchikey_block1)
    raw["product_block1"] = raw["product_inchikey"].map(kio.inchikey_block1)
    raw = raw[raw["substrate_block1"].ne("") & raw["product_block1"].ne("")]
    raw = raw[raw["substrate_block1"].ne(raw["product_block1"])].copy()
    raw["substrate_scaffold"] = raw["substrate_smiles"].map(kio.murcko_scaffold)
    raw["leakage_flag"] = raw["leakage_raw"].map(map_leakage)
    raw["module"] = raw["substrate_smiles"].map(lambda s: mr.classify(s)[0])

    # 去重：同 (substrate_block1, product_block1) 取低泄漏
    leak_rank = {"low": 0, "medium": 1, "medium_high": 2, "high": 3}
    raw["_lr"] = raw["leakage_flag"].map(leak_rank).fillna(1)
    food = raw.sort_values("_lr").drop_duplicates(["substrate_block1", "product_block1"], keep="first")
    keep = ["substrate_smiles", "substrate_inchikey", "substrate_block1", "substrate_scaffold",
            "product_smiles", "product_inchikey", "product_block1",
            "module", "label_source", "leakage_flag", "reaction_category"]
    food = food[keep].reset_index(drop=True)

    kio.write_table(food, MODULAR_OUT / "food_clean.parquet", force=args.force)
    kio.log.info("clean single-step food positives: %d (substrates=%d)",
                 len(food), food["substrate_block1"].nunique())
    kio.log.info("== per module (edges / unique substrates) ==")
    for m in mr.MODULES:
        sub = food[food["module"] == m]
        kio.log.info("  %s %-18s : %3d edges  %3d substrates  | leakage %s",
                     m, mr.MODULE_NAMES[m], len(sub), sub["substrate_block1"].nunique(),
                     sub["leakage_flag"].value_counts().to_dict())


if __name__ == "__main__":
    main()
