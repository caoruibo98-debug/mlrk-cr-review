"""
ltr_clean_eval.py — 杜绝作弊的最终评测

候选来自 clean_candidates(纯规则生成、EC 一致、无注入、无 vocab-NN)。
切分：底物 scaffold + 产物不相交。多 seed。
方法：
  - tanimoto         : 相似度基线(不学习)
  - random           : 随机下限
  - ec_only          : 只用 'EC!=0' 排序 —— **作弊检测器**：干净时应≈random，泄漏时会虚高
  - LTR_chem         : 只结构(ChemBERTa)
  - LTR_full         : 结构+EC+Tanimoto特征
  - ensemble         : LTR_full 与 tanimoto 的 rank 平均(稳健"两者之长")
判据：LTR 要**同时** > tanimoto 且 > ec_only 才算真信号。

用法：python scripts/ssrf/ml_ranking_kernel/modular/ltr_clean_eval.py [--seeds 2] [--force]
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit.Chem import DataStructs
from rdkit import RDLogger
from sklearn.model_selection import GroupKFold

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import kio  # noqa: E402
from ltr_build import LTR_OUT  # noqa: E402
from ltr_train_eval import encode_all  # noqa: E402
from ltr_v3_train_eval import fp, build_feats, train, within_sub_eval  # noqa: E402

RDLogger.DisableLog("rdApp.*")
WEAK_TRAIN_CAP = 700
VARIANTS = [("LTR_full", True, True), ("LTR_chem", False, False)]


def ranks_of(scores):
    s = np.where(np.isnan(scores), -1e9, scores)
    order = np.argsort(-s, kind="mergesort"); r = np.empty(len(s)); r[order] = np.arange(1, len(s) + 1)
    return r


def stat(v):
    return f"{np.mean(v):.3f}±{np.std(v):.3f}" if v else "—"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", type=int, default=2)
    ap.add_argument("--candidates", default="clean_candidates.parquet")
    ap.add_argument("--tag", default="clean")
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()
    kio.setup_logging()
    C = pd.read_parquet(LTR_OUT / args.candidates)
    if "tier" not in C.columns:
        C["tier"] = np.where(C["is_clean"], "gold", "weak")
    C["scaf"] = C["substrate_scaffold"].fillna("NS_" + C["sb"])
    fpc = {s: fp(s) for s in set(C.substrate_smiles) | set(C.product_smiles)}
    emb, dim = encode_all(set(C.substrate_smiles) | set(C.product_smiles))
    kio.log.info("clean candidates=%d subs=%d median/sub=%.0f", len(C), C.sb.nunique(), C.groupby("sb").size().median())

    def tan_score(g):
        sf = fpc.get(g.substrate_smiles.iloc[0])
        return np.array([DataStructs.TanimotoSimilarity(sf, fpc.get(p)) if (sf and fpc.get(p)) else 0.0
                         for p in g["product_smiles"]])
    def ec_score(g):
        return (g["ec_class"].to_numpy() != 0).astype(float)

    rows = []; imp_rows = []
    METHODS = ["tanimoto", "random", "ec_only", "LTR_chem", "LTR_full", "ensemble"]
    for m in ["A", "B", "C", "D"]:
        Cm = C[C.module == m].copy()
        testpool = Cm[Cm.tier.isin(["gold", "silver"])].drop_duplicates("sb")
        if testpool.sb.nunique() < 8:
            kio.log.info("module %s: too few test subs (%d) skip", m, testpool.sb.nunique()); continue
        kio.log.info("== module %s: test_subs=%d (gold=%d) ==", m, len(testpool), (testpool.tier == "gold").sum())
        coll = {meth: {"mrr": [], "r@1": [], "r@5": []} for meth in METHODS}
        imp_acc = None
        for seed in range(args.seeds):
            ts = testpool.sample(frac=1.0, random_state=seed).reset_index(drop=True)
            k = min(5, ts["scaf"].nunique()); ts["fold"] = -1
            for f, (_, te) in enumerate(GroupKFold(k).split(ts, groups=ts["scaf"])):
                ts.iloc[te, ts.columns.get_loc("fold")] = f
            fold_of = dict(zip(ts.sb, ts.fold)); Cm["fold"] = Cm.sb.map(fold_of)
            weak_subs = Cm.loc[Cm.tier == "weak", "sb"].unique()
            rng = np.random.default_rng(seed)
            keep_weak = set(rng.choice(weak_subs, min(WEAK_TRAIN_CAP, len(weak_subs)), replace=False)) if len(weak_subs) else set()
            for f in range(k):
                test_sb = set(ts.loc[ts.fold == f, "sb"]); test = Cm[Cm.sb.isin(test_sb)]
                test_true = set(test.loc[test.y == 1, "pb"])
                tr = Cm[(Cm["fold"].isna()) | (Cm["fold"] != f)]
                tr = tr[(tr.tier != "weak") | (tr.sb.isin(keep_weak))]
                tr = tr[~((tr.y == 1) & (tr.pb.isin(test_true)))]        # 产物不相交
                if test.sb.nunique() == 0 or (tr.y == 1).sum() < 5:
                    continue
                models = {}
                for name, ue, ut in VARIANTS:
                    Xtr, ok = build_feats(tr, emb, dim, fpc, ue, ut); trk = tr[ok]; qid = pd.factorize(trk.sb)[0]
                    mdl = train(Xtr[ok], trk.y.to_numpy(np.float32), qid); models[name] = (mdl, ue, ut)
                    if name == "LTR_full" and imp_acc is None:
                        imp_acc = mdl.feature_importances_
                def ml_score(g, key):
                    mdl, ue, ut = models[key]; Xt, okt = build_feats(g, emb, dim, fpc, ue, ut)
                    s = np.full(len(g), -1e9)
                    if okt.any(): s[okt] = mdl.predict(Xt[okt])
                    return s
                scorers = {
                    "tanimoto": tan_score, "ec_only": ec_score,
                    "random": (lambda g, _r=np.random.default_rng(seed * 9 + f): _r.random(len(g))),
                    "LTR_chem": lambda g: ml_score(g, "LTR_chem"),
                    "LTR_full": lambda g: ml_score(g, "LTR_full"),
                    "ensemble": lambda g: -(ranks_of(ml_score(g, "LTR_full")) + ranks_of(tan_score(g))),
                }
                for meth, sfn in scorers.items():
                    r = within_sub_eval(sfn, test)
                    for kk in ["mrr", "r@1", "r@5"]:
                        coll[meth][kk].append(r[kk])
        for meth, d in coll.items():
            if d["mrr"]:
                rows.append({"module": m, "method": meth, "n_pts": len(d["mrr"]),
                             "mrr": stat(d["mrr"]), "r@1": stat(d["r@1"]), "r@5": stat(d["r@5"])})
        if imp_acc is not None:
            n = 4 * dim
            imp_rows.append({"module": m, "structure": round(float(imp_acc[:n].sum()), 3),
                             "ec": round(float(imp_acc[n:n+7].sum()), 3), "tanimoto": round(float(imp_acc[n+7:].sum()), 3)})

    res = pd.DataFrame(rows); kio.write_table(res, LTR_OUT / f"{args.tag}_metrics.csv", force=True)
    imp = pd.DataFrame(imp_rows)
    if len(imp): kio.write_table(imp, LTR_OUT / f"{args.tag}_feature_importance.csv", force=True)
    print(f"\n=== 杜绝作弊评测 [{args.tag}] (scaffold+product-disjoint, multi-seed) ===")
    print(res.to_string(index=False))
    if len(imp): print("\n特征组重要性:\n", imp.to_string(index=False))


if __name__ == "__main__":
    main()
