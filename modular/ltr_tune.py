"""
ltr_tune.py — 嵌套 CV 调参（诚实优化 LTR_chem 准确性）

每模块：外层 scaffold-CV(报 test)；外层训练集再内层切 train/val，grid 选参(按 val MRR)，
回训外层训练集，评外层 test。对比 {tuned LTR_chem, default LTR_chem, tanimoto}。
**调参只在 val 上做,从不看 test** → 不作弊。产物不相交切分照旧。

输出 tuned 参数(modular/ltr/tuned_params.json) + 对比表(tune_metrics.csv)。
用法：python scripts/ssrf/ml_ranking_kernel/modular/ltr_tune.py [--candidates clean_candidates_full.parquet] [--force]
"""
from __future__ import annotations

import argparse
import json
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
from ltr_v3_train_eval import fp, build_feats, within_sub_eval  # noqa: E402

RDLogger.DisableLog("rdApp.*")
DEFAULT = dict(n_estimators=300, max_depth=6, learning_rate=0.08, colsample_bytree=0.5, min_child_weight=2)
GRID = [
    dict(n_estimators=200, max_depth=4, learning_rate=0.10, colsample_bytree=0.6, min_child_weight=4),
    dict(n_estimators=300, max_depth=6, learning_rate=0.08, colsample_bytree=0.5, min_child_weight=2),
    dict(n_estimators=500, max_depth=6, learning_rate=0.05, colsample_bytree=0.4, min_child_weight=3),
    dict(n_estimators=400, max_depth=8, learning_rate=0.04, colsample_bytree=0.5, min_child_weight=5),
    dict(n_estimators=600, max_depth=5, learning_rate=0.05, colsample_bytree=0.6, min_child_weight=2),
]


def train_p(X, y, qid, params):
    import xgboost as xgb
    o = np.argsort(qid, kind="mergesort")
    m = xgb.XGBRanker(objective="rank:ndcg", subsample=0.8, n_jobs=12, tree_method="hist", **params)
    m.fit(X[o], y[o], qid=qid[o]); return m


def fit_eval(tr, te, emb, dim, fpc, params):
    Xtr, ok = build_feats(tr, emb, dim, fpc, False, False); trk = tr[ok]
    mdl = train_p(Xtr[ok], trk.y.to_numpy(np.float32), pd.factorize(trk.sb)[0], params)
    def sf(g):
        Xt, okt = build_feats(g, emb, dim, fpc, False, False); s = np.full(len(g), -1e9)
        if okt.any(): s[okt] = mdl.predict(Xt[okt])
        return s
    return within_sub_eval(sf, te), mdl


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--candidates", default="clean_candidates_full.parquet")
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()
    kio.setup_logging()
    C = pd.read_parquet(LTR_OUT / args.candidates)
    if "tier" not in C.columns:
        C["tier"] = np.where(C.is_clean, "gold", "weak")
    C["scaf"] = C.substrate_scaffold.fillna("NS_" + C.sb)
    fpc = {s: fp(s) for s in set(C.substrate_smiles) | set(C.product_smiles)}
    emb, dim = encode_all(set(C.substrate_smiles) | set(C.product_smiles))

    def tan(g):
        sf = fpc.get(g.substrate_smiles.iloc[0])
        return np.array([DataStructs.TanimotoSimilarity(sf, fpc.get(p)) if (sf and fpc.get(p)) else 0.0 for p in g.product_smiles])

    rows = []; tuned_params = {}
    for m in ["A", "B", "C", "D"]:
        Cm = C[C.module == m].copy()
        test_subs = Cm[Cm.tier.isin(["gold", "silver"])].drop_duplicates("sb").reset_index(drop=True)
        k = min(4, test_subs.scaf.nunique())
        if k < 2 or len(test_subs) < 8:
            continue
        test_subs["fold"] = -1
        for f, (_, te) in enumerate(GroupKFold(k).split(test_subs, groups=test_subs.scaf)):
            test_subs.iloc[te, test_subs.columns.get_loc("fold")] = f
        fold_of = dict(zip(test_subs.sb, test_subs.fold)); Cm["fold"] = Cm.sb.map(fold_of)
        wk = Cm.loc[Cm.tier == "weak", "sb"].unique()
        keep_wk = set(np.random.default_rng(0).choice(wk, min(700, len(wk)), replace=False)) if len(wk) else set()
        tun, dft, tn, chosen = [], [], [], {}
        for f in range(k):
            test = Cm[(Cm.fold == f) & Cm.tier.isin(["gold", "silver"])]
            ttrue = set(test.loc[test.y == 1, "pb"])
            otr = Cm[(Cm.fold.isna()) | (Cm.fold != f)]
            otr = otr[(otr.tier != "weak") | (otr.sb.isin(keep_wk))]
            otr = otr[~((otr.y == 1) & (otr.pb.isin(ttrue)))]   # 产物不相交(对 test)
            if test.sb.nunique() == 0 or (otr.y == 1).sum() < 8:
                continue
            # 内层：从 otr 的 gold+silver 底物切 val 选参
            inner = otr[otr.tier.isin(["gold", "silver"])].drop_duplicates("sb")
            ik = min(3, inner.scaf.nunique())
            best_p, best_v = DEFAULT, -1
            if ik >= 2:
                inner = inner.reset_index(drop=True); inner["ifold"] = -1
                for jf, (_, ite) in enumerate(GroupKFold(ik).split(inner, groups=inner.scaf)):
                    inner.iloc[ite, inner.columns.get_loc("ifold")] = jf
                vsb = set(inner.loc[inner.ifold == 0, "sb"])
                itr = otr[~otr.sb.isin(vsb)]; ival = otr[otr.sb.isin(vsb) & otr.tier.isin(["gold", "silver"])]
                ivtrue = set(ival.loc[ival.y == 1, "pb"]); itr = itr[~((itr.y == 1) & (itr.pb.isin(ivtrue)))]
                if ival.sb.nunique() and (itr.y == 1).sum() >= 8:
                    for p in GRID:
                        r, _ = fit_eval(itr, ival, emb, dim, fpc, p)
                        if r["mrr"] > best_v:
                            best_v, best_p = r["mrr"], p
            chosen[f] = best_p
            rt, _ = fit_eval(otr, test, emb, dim, fpc, best_p); tun.append(rt)
            rd, _ = fit_eval(otr, test, emb, dim, fpc, DEFAULT); dft.append(rd)
            tn.append(within_sub_eval(tan, test))
        def ag(L, key): return float(np.mean([x[key] for x in L])) if L else 0.0
        for name, L in [("LTR_tuned", tun), ("LTR_default", dft), ("tanimoto", tn)]:
            rows.append({"module": m, "method": name, "n_folds": len(L),
                         "mrr": round(ag(L, "mrr"), 3), "r@1": round(ag(L, "r@1"), 3), "r@5": round(ag(L, "r@5"), 3)})
        # 选最常被选中的参数作该模块 tuned 参数
        from collections import Counter
        cc = Counter(json.dumps(p, sort_keys=True) for p in chosen.values())
        tuned_params[m] = json.loads(cc.most_common(1)[0][0]) if cc else DEFAULT
        kio.log.info("module %s tuned params: %s", m, tuned_params[m])

    res = pd.DataFrame(rows); kio.write_table(res, LTR_OUT / "tune_metrics.csv", force=True)
    (LTR_OUT / "tuned_params.json").write_text(json.dumps(tuned_params, indent=2), encoding="utf-8")
    print("\n=== 调参对比 (nested-CV, LTR_chem) ===\n", res.to_string(index=False))
    print("\ntuned params:", json.dumps(tuned_params, indent=2))


if __name__ == "__main__":
    main()
