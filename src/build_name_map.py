"""
build_name_map.py — 结构(InChIKey block1) → 名称 索引（一次性，缓存）

合并本地库（unified 优先 → ChEBI → FooDB），按 InChIKey block1（连接层）建名称表，
供 predict_substrate 给预测产物加名字。block1 级=忽略立体，名称为 best-effort。

输出 outputs/features/inchikey_name_map.parquet (block1, name, source)
用法：python scripts/ssrf/ml_ranking_kernel/src/build_name_map.py [--force]
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import kio

SOURCES = [  # (path, inchikey_col, name_col, source_tag) —— 优先级从高到低
    ("Dataset/unified/compound_master_index.csv", "standard_inchikey", "preferred_name", "unified"),
    ("Dataset/chebi/chebi_compounds.csv", "inchikey", "name", "chebi"),
    ("Dataset/foodb/foodb_compounds.csv", "inchikey", "name", "foodb"),
]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()
    kio.setup_logging()

    best: dict[str, tuple[str, str]] = {}   # block1 -> (name, source)
    for rel, ikc, nmc, tag in SOURCES:
        p = kio.resolve_input(rel)
        if not p.exists():
            kio.log.warning("skip missing %s", p)
            continue
        df = pd.read_csv(p, usecols=lambda c: c in (ikc, nmc), dtype=str, low_memory=False).dropna()
        added = 0
        for ik, nm in zip(df[ikc], df[nmc]):
            b1 = kio.inchikey_block1(ik)
            nm = str(nm).strip()
            if not b1 or not nm or nm.lower() in ("nan", "none"):
                continue
            if b1 not in best:          # first-wins（高优先级源先填）
                best[b1] = (nm, tag)
                added += 1
        kio.log.info("%-8s %s rows → +%d new block1 (total %d)", tag, len(df), added, len(best))

    out = pd.DataFrame([(b, n, s) for b, (n, s) in best.items()], columns=["block1", "name", "source"])
    kio.write_table(out, kio.FEATURES_DIR / "inchikey_name_map.parquet", force=args.force)
    kio.log.info("name map: %d block1 entries", len(out))


if __name__ == "__main__":
    main()
