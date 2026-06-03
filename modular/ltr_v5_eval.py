"""
ltr_v5_eval.py — ① 银标提升 + 多seed + 时间留出(新文献召回)

复用 v4 去污染候选(真实decoy + 后面 product-disjoint 切分)，把测试集从 gold(80) 扩到
gold+silver(~4000)，多 seed 收窄 CI，并做时间留出回答"能不能召回新文献代谢物"。

报告两套 scope：gold(纯文献金标) 和 gold+silver(加DB银标，标泄漏)；× {tanimoto, random, LTR_chem, LTR_full}。
+ 特征组重要性 + 时间留出(train<2023 → test≥2023)。

用法：python scripts/ssrf/ml_ranking_kernel/modular/ltr_v5_eval.py [--seeds 3] [--force]
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
from ltr_train_eval import encode_all, TOPK  # noqa: E402
from ltr_v3_train_eval import fp, build_feats, train, within_sub_eval  # noqa: E402
from ltr_v4_train_eval import build_realistic_candidates  # noqa: E402

RDLogger.DisableLog("rdApp.*")
WEAK_TRAIN_CAP = 700   # 每模块训练用 weak 底物上限（控算力）
VARIANTS = [("LTR_full", True, True), ("LTR_chem", False, False)]


def stat(vals):
    return f"{np.mean(vals):.3f}±{np.std(vals):.3f}" if vals else "—"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", type=int, default=3)
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()
    kio.setup_logging()
    cand = pd.read_parquet(LTR_OUT / "candidates.parquet")
    vocab = pd.read_parquet(LTR_OUT / "vocab.parquet")
    if "tier" not in cand.columns:
        cand["tier"] = np.where(cand["is_clean"], "gold", "weak")
    allsm = set(cand.substrate_smiles) | set(cand.product_smiles) | set(vocab.product_smiles)
    fpc = {s: fp(s) for s in allsm}
    emb, dim = encode_all(allsm)
    C = build_realistic_candidates(cand, vocab, fpc, np.random.default_rng(0))
    # tier 没被 build 带出来时从 cand 回填
    if "tier" not in C.columns:
        tmap = cand.drop_duplicates("sb").set_index("sb")["tier"].to_dict()
        C["tier"] = C["sb"].map(tmap).fillna("weak")
    C["scaf"] = C["substrate_scaffold"].fillna("NS_" + C["sb"])
    kio.log.info("realistic candidates=%d (median/sub=%.0f)", len(C), C.groupby("sb").size().median())

    rows = []; imp_rows = []
    def tan_score(g):
        sf = fpc.get(g.substrate_smiles.iloc[0])
        return np.array([DataStructs.TanimotoSimilarity(sf, fpc.get(p)) if (sf and fpc.get(p)) else 0.0
                         for p in g["product_smiles"]])

    for m in ["A", "B", "C", "D"]:
        Cm = C[C.module == m].copy()
        testpool = Cm[Cm.tier.isin(["gold", "silver"])].drop_duplicates("sb")
        if testpool.sb.nunique() < 8:
            kio.log.info("module %s: too few gold+silver subs (%d) skip", m, testpool.sb.nunique()); continue
        kio.log.info("== module %s: gold=%d silver=%d test_subs=%d ==",
                     m, (testpool.tier == "gold").sum(), (testpool.tier == "silver").sum(), len(testpool))
        # 收集：per method × scope(gold / gold+silver) 的 per-(seed,fold) 指标
        coll = {(name, sc): {"mrr": [], "r@1": [], "r@5": []} for name in
                ["tanimoto", "random"] + [v[0] for v in VARIANTS] for sc in ["gold", "gold+silver"]}
        imp_acc = None
        for seed in range(args.seeds):
            ts = testpool.sample(frac=1.0, random_state=seed).reset_index(drop=True)
            k = min(5, ts["scaf"].nunique())
            ts["fold"] = -1
            for f, (_, te) in enumerate(GroupKFold(k).split(ts, groups=ts["scaf"])):
                ts.iloc[te, ts.columns.get_loc("fold")] = f
            fold_of = dict(zip(ts.sb, ts.fold)); Cm["fold"] = Cm.sb.map(fold_of)
            weak_subs = Cm.loc[Cm.tier == "weak", "sb"].unique()
            rng = np.random.default_rng(seed)
            keep_weak = set(rng.choice(weak_subs, min(WEAK_TRAIN_CAP, len(weak_subs)), replace=False)) if len(weak_subs) else set()
            for f in range(k):
                test_sb = set(ts.loc[ts.fold == f, "sb"])
                test = Cm[Cm.sb.isin(test_sb)]
                gold_sb = set(ts.loc[(ts.fold == f) & (ts.tier == "gold"), "sb"])
                test_true = set(test.loc[test.y == 1, "pb"])
                tr = Cm[(Cm["fold"].isna() | (Cm["fold"] != f))]
                tr = tr[(tr.tier != "weak") | (tr.sb.isin(keep_weak))]
                tr = tr[~((tr.y == 1) & (tr.pb.isin(test_true)))]    # product-disjoint
                if test.sb.nunique() == 0 or (tr.y == 1).sum() < 5:
                    continue
                scorers = {"tanimoto": tan_score,
                           "random": (lambda g, _r=np.random.default_rng(seed * 10 + f): _r.random(len(g)))}
                for name, ue, ut in VARIANTS:
                    Xtr, ok = build_feats(tr, emb, dim, fpc, ue, ut); trk = tr[ok]; qid = pd.factorize(trk.sb)[0]
                    model = train(Xtr[ok], trk.y.to_numpy(np.float32), qid)
                    def sf(g, _m=model, _ue=ue, _ut=ut):
                        Xt, okt = build_feats(g, emb, dim, fpc, _ue, _ut); sc = np.full(len(g), -1e9)
                        if okt.any(): sc[okt] = _m.predict(Xt[okt])
                        return sc
                    scorers[name] = sf
                    if name == "LTR_full" and imp_acc is None:
                        imp_acc = model.feature_importances_
                for name, sfn in scorers.items():
                    for sc, sub in [("gold+silver", test), ("gold", test[test.sb.isin(gold_sb)])]:
                        if sub.sb.nunique() == 0:
                            continue
                        r = within_sub_eval(sfn, sub)
                        for kk in ["mrr", "r@1", "r@5"]:
                            coll[(name, sc)][kk].append(r[kk])
        for (name, sc), d in coll.items():
            if d["mrr"]:
                rows.append({"module": m, "method": name, "scope": sc, "n_folds": len(d["mrr"]),
                             "mrr": stat(d["mrr"]), "r@1": stat(d["r@1"]), "r@5": stat(d["r@5"])})
        if imp_acc is not None:
            n = 4 * dim
            imp_rows.append({"module": m, "structure": round(float(imp_acc[:n].sum()), 3),
                             "ec": round(float(imp_acc[n:n+7].sum()), 3),
                             "tanimoto": round(float(imp_acc[n+7:].sum()), 3)})

        # 时间留出（新文献召回，gold+silver，LTR_full）
        fut = testpool[pd.to_numeric(testpool["year"], errors="coerce") >= 2023]
        if fut.sb.nunique() >= 3:
            fut_sb = set(fut.sb); ft = Cm[Cm.sb.isin(fut_sb)]; ft_true = set(ft.loc[ft.y == 1, "pb"])
            tr = Cm[~Cm.sb.isin(fut_sb)]; tr = tr[~((tr.y == 1) & (tr.pb.isin(ft_true)))]
            Xtr, ok = build_feats(tr, emb, dim, fpc, True, True); trk = tr[ok]; qid = pd.factorize(trk.sb)[0]
            model = train(Xtr[ok], trk.y.to_numpy(np.float32), qid)
            def sf2(g, _m=model):
                Xt, okt = build_feats(g, emb, dim, fpc, True, True); s = np.full(len(g), -1e9)
                if okt.any(): s[okt] = _m.predict(Xt[okt])
                return s
            r = within_sub_eval(sf2, ft)
            rt = within_sub_eval(tan_score, ft)
            rows.append({"module": m, "method": "LTR_full", "scope": "temporal>=2023",
                         "n_folds": r["n"], "mrr": round(r["mrr"], 3), "r@1": round(r["r@1"], 3), "r@5": round(r["r@5"], 3)})
            rows.append({"module": m, "method": "tanimoto", "scope": "temporal>=2023",
                         "n_folds": rt["n"], "mrr": round(rt["mrr"], 3), "r@1": round(rt["r@1"], 3), "r@5": round(rt["r@5"], 3)})

    res = pd.DataFrame(rows); kio.write_table(res, LTR_OUT / "v5_metrics.csv", force=True)
    imp = pd.DataFrame(imp_rows)
    if len(imp): kio.write_table(imp, LTR_OUT / "v5_feature_importance.csv", force=True)
    print("\n=== v5 (gold+silver, multi-seed, de-confounded) ===")
    print(res[res.scope.isin(["gold", "gold+silver"])].sort_values(["module", "scope", "method"]).to_string(index=False))
    tmp = res[res.scope.str.startswith("temporal")]
    if len(tmp):
        print("\n=== 新文献召回 (temporal≥2023) ===\n", tmp[["module", "method", "n_folds", "mrr", "r@1", "r@5"]].to_string(index=False))
    if len(imp):
        print("\n=== 特征组重要性 ===\n", imp.to_string(index=False))


if __name__ == "__main__":
    main()
