"""
build_network_negatives.py — Phase 1b（阶段化，不阻塞首版训练）

从 AGREDA per-species GEM 推导**代谢网络未覆盖隐式负样本**：
  底物 S 在网络中存在（KEGG→InChIKey 映射到某网络代谢物），候选产物 P 也在网络中，
  但网络里**没有任何反应**实现 S→P 这条有向转化 → 视为隐式负（独立于规则源，泄漏干净）。

桥接：AGREDA metKEGGID → (kegg_conv_chebi) → ChEBI → InChIKey block1。

产出 outputs/datasets/network_negatives.parquet，供 build_kernel_dataset 以
--with-network-negatives 并入重训，并量化对 AUPRC/泛化 gap 的提升。

用法：
  python scripts/ssrf/ml_ranking_kernel/src/build_network_negatives.py [--max-per-substrate 40] [--force]
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.io as sio

sys.path.insert(0, str(Path(__file__).resolve().parent))
import kio

AGREDA_SPECIES = kio.PROJECT_ROOT / "Dataset" / "agreda" / "AGREDA_v1.0.0" / "AGREDA_v1.0.0" / "mat" / "species"
KEGG_CHEBI = kio.PROJECT_ROOT / "Dataset" / "kegg" / "kegg_conv_chebi_compound.tsv"
CHEBI_COMPOUNDS = kio.PROJECT_ROOT / "Dataset" / "chebi" / "chebi_compounds.csv"


def _clean_kegg(v) -> str:
    if isinstance(v, np.ndarray):
        v = v.item() if v.size == 1 else ""
    s = str(v).strip()
    return s if re.match(r"^C\d{5}$", s) else ""


def build_kegg_to_block1() -> dict[str, str]:
    """KEGG Cxxxxx → InChIKey block1，经 ChEBI 桥。"""
    if not (KEGG_CHEBI.exists() and CHEBI_COMPOUNDS.exists()):
        kio.log.warning("KEGG/ChEBI bridge files missing; KEGG->block1 map will be empty")
        return {}
    kc = pd.read_csv(KEGG_CHEBI, sep="\t", header=None, names=["kegg", "chebi"], dtype=str)
    kc["kegg"] = kc["kegg"].str.replace("cpd:", "", regex=False).str.strip()
    kc["chebi"] = kc["chebi"].str.replace("chebi:", "", regex=False, case=False).str.upper().str.strip()
    chebi = pd.read_csv(CHEBI_COMPOUNDS, dtype=str, low_memory=False)
    # 找 chebi id 列 + inchikey 列
    id_col = next((c for c in chebi.columns if c.lower() in ("chebi_id", "chebiid", "id")), None)
    ik_col = next((c for c in chebi.columns if "inchikey" in c.lower()), None)
    if id_col is None or ik_col is None:
        kio.log.warning("chebi_compounds.csv missing id/inchikey columns: %s", list(chebi.columns))
        return {}
    chebi["_cid"] = chebi[id_col].astype(str).str.upper().str.replace("CHEBI:", "", regex=False).str.strip()
    chebi["_b1"] = chebi[ik_col].map(kio.inchikey_block1)
    chebi = chebi[chebi["_b1"].ne("")]
    kc["chebi"] = kc["chebi"].str.replace("CHEBI:", "", regex=False)
    merged = kc.merge(chebi[["_cid", "_b1"]], left_on="chebi", right_on="_cid", how="inner")
    m = dict(zip(merged["kegg"], merged["_b1"]))
    kio.log.info("KEGG->block1 bridge: %d entries", len(m))
    return m


def aggregate_agreda(limit: int | None = None):
    """返回 (network_block1 set, directed (b1_sub, b1_prod) edge set)，需 kegg->block1。"""
    kegg2b1 = build_kegg_to_block1()
    files = sorted(AGREDA_SPECIES.glob("*.mat"))
    if limit:
        files = files[:limit]
    kio.log.info("aggregating %d AGREDA species models ...", len(files))
    net_b1: set[str] = set()
    edges: set[tuple[str, str]] = set()
    for i, f in enumerate(files):
        try:
            m = sio.loadmat(f, struct_as_record=False, squeeze_me=True)
            model = m["model"]
            met_kegg = [ _clean_kegg(x) for x in np.atleast_1d(model.metKEGGID) ]
            met_b1 = [ kegg2b1.get(k, "") for k in met_kegg ]
            S = model.S.tocoo()
            lb = np.atleast_1d(getattr(model, "lb", np.zeros(S.shape[1])))
            ub = np.atleast_1d(getattr(model, "ub", np.ones(S.shape[1])))
            for b1 in met_b1:
                if b1:
                    net_b1.add(b1)
            # per reaction: reactants (S<0), products (S>0)
            from collections import defaultdict
            rxn_neg = defaultdict(list)
            rxn_pos = defaultdict(list)
            for r, c, v in zip(S.row, S.col, S.data):
                b1 = met_b1[r] if r < len(met_b1) else ""
                if not b1:
                    continue
                (rxn_neg if v < 0 else rxn_pos)[c].append(b1)
            for c in set(list(rxn_neg) + list(rxn_pos)):
                # 方向：可逆(lb<0 且 ub>0)两向都加；正向只 reactant→product；逆向只 product→reactant
                fwd = (c >= len(lb)) or (ub[c] > 0)
                rev = (c < len(lb)) and (lb[c] < 0)
                for a in rxn_neg.get(c, []):
                    for b in rxn_pos.get(c, []):
                        if a == b:
                            continue
                        if fwd:
                            edges.add((a, b))
                        if rev:
                            edges.add((b, a))
        except Exception as e:
            kio.log.debug("skip %s: %s", f.name, e)
        if (i + 1) % 200 == 0:
            kio.log.info("  %d/%d species  (net_mets=%d edges=%d)", i + 1, len(files), len(net_b1), len(edges))
    kio.log.info("AGREDA network: %d metabolites, %d directed edges", len(net_b1), len(edges))
    return net_b1, edges


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-per-substrate", type=int, default=40)
    ap.add_argument("--limit-species", type=int, default=None)
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()
    kio.setup_logging()

    cand_p = kio.DATASETS_DIR / "substrate_candidates.parquet"
    edges_p = kio.DATASETS_DIR / "kernel_edges.parquet"
    if not cand_p.exists():
        raise SystemExit("substrate_candidates.parquet missing; run generate_candidates.py first.")
    cand = pd.read_parquet(cand_p).rename(columns={"rclss_production_score": "rclss_prior"})
    edges = pd.read_parquet(edges_p)
    pos = edges[edges["y"] == 1]
    true_prod = pos.groupby("substrate_block1")["product_block1"].agg(set).to_dict()

    net_b1, net_edges = aggregate_agreda(args.limit_species)
    if not net_b1:
        raise SystemExit("Empty AGREDA network (bridge failed). Check KEGG/ChEBI files.")

    rng = np.random.default_rng(int(kio.load_config()["splits"]["random_seed"]))
    rows = []
    for sub, g in cand.groupby("substrate_block1"):
        if sub not in net_b1:
            continue  # 底物不在网络 → 无法判隐式负
        forbid = true_prod.get(sub, set())
        cn = g[g["product_block1"].isin(net_b1)]            # 产物也在网络
        cn = cn[~cn["product_block1"].isin(forbid)]
        cn = cn[~cn.apply(lambda r: (sub, r["product_block1"]) in net_edges, axis=1)]  # 网络无此转化
        cn = cn.drop_duplicates("product_block1")
        if len(cn) > args.max_per_substrate:
            cn = cn.sample(n=args.max_per_substrate, random_state=int(rng.integers(0, 2**31 - 1)))
        for _, r in cn.iterrows():
            rows.append({"substrate_smiles": r["substrate_smiles"], "substrate_inchikey": r.get("substrate_inchikey"),
                         "substrate_block1": sub, "product_smiles": r["product_smiles"],
                         "product_inchikey": r.get("product_inchikey"), "product_block1": r["product_block1"],
                         "y": 0, "label_source": "network_implicit_negative", "leakage_flag": "network",
                         "confidence": 1.0, "rclss_prior": r.get("rclss_prior")})
    out = pd.DataFrame(rows)
    kio.write_table(out, kio.DATASETS_DIR / "network_negatives.parquet", force=args.force)
    kio.log.info("network-implicit negatives: %d over %d substrates",
                 len(out), out["substrate_block1"].nunique() if len(out) else 0)


if __name__ == "__main__":
    main()
