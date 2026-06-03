"""
prepare_gutmgene.py — ② 准备：从 gutMGene 抽独立 substrate→product 反应（尤其补 C 蛋白）

gutMGene 是文献策展的菌→代谢物（独立于 MicrobeRX/RetroRules），human 表里有
substrate + metabolite + metsmiles(产物SMILES) + pmid + strain。产物结构现成，只需解析底物名。

输出 modular/outputs/ltr/gutmgene_reactions.parquet（可作独立测试层），并报每模块产出。
用法：python scripts/ssrf/ml_ranking_kernel/modular/prepare_gutmgene.py [--limit N] [--force]
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd
from rdkit import RDLogger

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import kio  # noqa: E402
import module_router as mr  # noqa: E402
from ltr_build import LTR_OUT  # noqa: E402

RDLogger.DisableLog("rdApp.*")
GM = kio.PROJECT_ROOT / "Dataset" / "gutmgene"
FILES = ["gutmgene_microbe_metabolite_human.csv", "gutmgene_metabolite_gene_human.csv",
         "gutmgene_microbe_gene_human.csv", "gutmgene_gene_lookup.csv"]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()
    kio.setup_logging()
    from resolve import name_to_smiles

    frames = []
    for f in FILES:
        p = GM / f
        if p.exists():
            frames.append(pd.read_csv(p, low_memory=False))
    d = pd.concat(frames, ignore_index=True)
    keep = ["substrate", "substrate_pubchem_cid", "metabolite", "metsmiles", "pmid", "strain", "gene"]
    d = d[[c for c in keep if c in d.columns]].copy()
    d = d.dropna(subset=["substrate", "metabolite"])
    d = d.drop_duplicates(["substrate", "metabolite"])
    if args.limit:
        d = d.head(args.limit)
    kio.log.info("gutMGene candidate rows (dedup substrate×metabolite): %d", len(d))

    # 产物结构现成；底物名→SMILES（缓存）
    cache = {}
    def resolve(name):
        if name in cache:
            return cache[name]
        smi, _ = name_to_smiles(str(name))
        cache[name] = smi
        return smi

    rows = []
    for i, r in enumerate(d.itertuples(index=False)):
        ms = getattr(r, "metsmiles", None)
        psmi = kio.canonical_smiles(ms) if isinstance(ms, str) and len(str(ms)) > 2 else None
        if not psmi:                       # metsmiles 多为空 → 解析 metabolite 名
            ps = resolve(r.metabolite)
            psmi = kio.canonical_smiles(ps) if ps else None
        ssmi = resolve(r.substrate)
        ssmi = kio.canonical_smiles(ssmi) if ssmi else None
        if not psmi or not ssmi:
            continue
        sb = kio.inchikey_block1(kio.smiles_to_inchikey(ssmi))
        pb = kio.inchikey_block1(kio.smiles_to_inchikey(psmi))
        if not sb or not pb or sb == pb:
            continue
        rows.append({"substrate_name": r.substrate, "substrate_smiles": ssmi, "substrate_block1": sb,
                     "product_name": getattr(r, "metabolite", ""), "product_smiles": psmi, "product_block1": pb,
                     "module": mr.classify(ssmi)[0], "pmid": getattr(r, "pmid", ""),
                     "strain": getattr(r, "strain", ""), "source": "gutmgene"})
        if (i + 1) % 100 == 0:
            kio.log.info("  resolved %d/%d", i + 1, len(d))
    out = pd.DataFrame(rows)
    if out.empty:
        kio.log.warning("gutMGene yielded 0 usable reactions (resolution failed)"); return
    out = out.drop_duplicates(["substrate_block1", "product_block1"])
    kio.write_table(out, LTR_OUT / "gutmgene_reactions.parquet", force=args.force)
    kio.log.info("== gutMGene independent reactions: %d (substrates=%d) ==", len(out), out["substrate_block1"].nunique())
    kio.log.info("per module: %s", out.groupby("module").size().to_dict())
    kio.log.info("with PMID: %d", int(out["pmid"].astype(str).str.len().gt(2).sum()))


if __name__ == "__main__":
    main()
