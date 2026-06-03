"""
ltr_v3_train_eval.py — 步骤1+2+3+4 合一：同底物候选判别 LTR（决定性实验）

任务变了：在**同底物候选集**（真产物 + 同底物 hard decoy）里把真产物排前面。
decoy 是 S 经规则 RunReactants 的错产物 → 结构也像 S → Tanimoto 难分 → 给 ML 留出空间。

特征：
  - 结构/反应：ChemBERTa [z_s,z_p,z_s−z_p,z_s⊙z_p]（1536）
  - 反应类型(步骤2)：生成规则的 EC class onehot（7）
  - Tanimoto(步骤3)：Morgan 相似度标量（1）——同时作显式特征 + 独立基线
对比：{tanimoto基线, random, LTR_chem, LTR_chem+ec, LTR_full(chem+ec+tan)}
评估：干净集 scaffold-CV，per-fold mean±std；+ 时间留出(步骤4)；+ 特征组重要性。
存终模型 → modular/outputs/ltr/models_v3/。

用法：python scripts/ssrf/ml_ranking_kernel/modular/ltr_v3_train_eval.py [--force]
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs
from rdkit import RDLogger

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import kio  # noqa: E402
from ltr_build import LTR_OUT  # noqa: E402
from ltr_train_eval import encode_all, rxn_feats, TOPK  # noqa: E402

RDLogger.DisableLog("rdApp.*")
BEST = dict(n_estimators=300, max_depth=6, learning_rate=0.08, colsample_bytree=0.5)


def fp(s):
    m = Chem.MolFromSmiles(s) if isinstance(s, str) else None
    return AllChem.GetMorganFingerprintAsBitVect(m, 2, 2048) if m else None


def build_feats(df, emb, dim, fpc, use_ec, use_tan):
    X, ok = rxn_feats(df["substrate_smiles"].tolist(), df["product_smiles"].tolist(), emb, dim)
    extra = []
    if use_ec:
        ec = df["ec_class"].to_numpy(int).clip(0, 6)
        oh = np.zeros((len(df), 7), np.float32); oh[np.arange(len(df)), ec] = 1.0
        extra.append(oh)
    if use_tan:
        tn = np.zeros((len(df), 1), np.float32)
        for i, (s, p) in enumerate(zip(df["substrate_smiles"], df["product_smiles"])):
            a, b = fpc.get(s), fpc.get(p)
            tn[i, 0] = DataStructs.TanimotoSimilarity(a, b) if (a and b) else 0.0
        extra.append(tn)
    if extra:
        X = np.hstack([X] + extra)
    return X, ok


def train(X, y, qid):
    import xgboost as xgb
    o = np.argsort(qid, kind="mergesort")
    m = xgb.XGBRanker(objective="rank:ndcg", subsample=0.8, min_child_weight=2,
                      n_jobs=12, tree_method="hist", **BEST)
    m.fit(X[o], y[o], qid=qid[o]); return m


def within_sub_eval(score_fn, test_df):
    """每底物在自己候选集里排，真产物 recall@k/MRR。score_fn(sub_df)->scores。
    并列分数用**随机** tie-break（加微小抖动），避免真产物因行序(注入在前)被高估。"""
    crs = []
    rng = np.random.default_rng(12345)
    for sb, g in test_df.groupby("sb"):
        g = g.reset_index(drop=True)
        sc = np.asarray(score_fn(g), dtype=float)
        sc = sc + rng.uniform(0, 1e-6, size=len(sc))   # 随机打破并列
        order = np.argsort(-sc, kind="mergesort")
        ranks = np.empty(len(sc)); ranks[order] = np.arange(1, len(sc) + 1)
        ytrue = g["y"].to_numpy()
        best = min((ranks[i] for i in range(len(g)) if ytrue[i] == 1), default=None)
        if best:
            crs.append(best)
    cr = np.array(crs, float)
    return {"n": len(cr), "mrr": float(np.mean(1/cr)) if len(cr) else 0,
            **{f"r@{k}": float(np.mean(cr <= k)) if len(cr) else 0 for k in TOPK}}


def agg(rows):
    out = {k: f"{np.mean([x[k] for x in rows]):.3f}±{np.std([x[k] for x in rows]):.3f}"
           for k in ["mrr", "r@1", "r@5", "r@10"]}
    out["n"] = int(sum(x["n"] for x in rows)); return out


def main() -> None:
    ap = argparse.ArgumentParser(); ap.add_argument("--force", action="store_true"); args = ap.parse_args()
    kio.setup_logging()
    cand = pd.read_parquet(LTR_OUT / "candidates.parquet")
    cand["scaf"] = cand["substrate_scaffold"].fillna("NS_" + cand["sb"])
    emb, dim = encode_all(set(cand["substrate_smiles"]) | set(cand["product_smiles"]))
    fpc = {s: fp(s) for s in set(cand["substrate_smiles"]) | set(cand["product_smiles"])}
    kio.log.info("candidates=%d emb_dim=%d", len(cand), dim)

    from sklearn.model_selection import GroupKFold
    rows = []; imp_rows = []
    models_dir = kio.safe_output_path("modular/ltr/models_v3/_m").parent
    VARIANTS = [("LTR_full", True, True), ("LTR_chem+ec", True, False), ("LTR_chem", False, False)]
    for m in ["A", "B", "C", "D"]:
        C = cand[cand["module"] == m].copy()
        clean = C[C["is_clean"]]
        csubs = clean.drop_duplicates("sb")
        k = min(5, csubs["scaf"].nunique())
        if k < 2 or len(csubs) < 5:
            kio.log.info("module %s: too few clean subs (%d), skip", m, len(csubs)); continue
        csubs = csubs.reset_index(drop=True); csubs["fold"] = -1
        for f, (_, te) in enumerate(GroupKFold(n_splits=k).split(csubs, groups=csubs["scaf"])):
            csubs.iloc[te, csubs.columns.get_loc("fold")] = f
        fold_of = dict(zip(csubs["sb"], csubs["fold"]))
        C["fold"] = C["sb"].map(fold_of).fillna(-1).astype(int)
        kio.log.info("== module %s: clean_subs=%d folds=%d median_cand=%.0f ==",
                     m, len(csubs), k, C[C.is_clean].groupby("sb").size().median())

        # 基线：Tanimoto / random（同候选集）
        tan_f, rnd_f = [], []
        per_var = {v[0]: [] for v in VARIANTS}
        imp_acc = None
        for f in range(k):
            test = C[(C["fold"] == f) & (C["is_clean"])]
            train_df = C[C["fold"] != f]  # 其它折干净 + 全 weak
            if test["sb"].nunique() == 0 or (train_df["y"] == 1).sum() < 5:
                continue
            # baselines
            tan_f.append(within_sub_eval(
                lambda g: np.array([DataStructs.TanimotoSimilarity(fpc.get(g.substrate_smiles.iloc[0]), fpc.get(p))
                                    if (fpc.get(g.substrate_smiles.iloc[0]) and fpc.get(p)) else 0.0
                                    for p in g["product_smiles"]]), test))
            rng = np.random.default_rng(f)
            rnd_f.append(within_sub_eval(lambda g: rng.random(len(g)), test))
            # LTR variants
            for name, use_ec, use_tan in VARIANTS:
                Xtr, ok = build_feats(train_df, emb, dim, fpc, use_ec, use_tan)
                trk = train_df[ok]; qid = pd.factorize(trk["sb"])[0]
                model = train(Xtr[ok], trk["y"].to_numpy(np.float32), qid)
                def sf(g, _m=model, _ue=use_ec, _ut=use_tan):
                    Xte, oke = build_feats(g, emb, dim, fpc, _ue, _ut)
                    sc = np.full(len(g), -1e9)
                    if oke.any():
                        sc[oke] = _m.predict(Xte[oke])
                    return sc
                per_var[name].append(within_sub_eval(sf, test))
                if name == "LTR_full" and imp_acc is None:
                    imp_acc = model.feature_importances_
        if tan_f:
            rows.append({"module": m, "method": "tanimoto_baseline", "eval": "scaffold_cv", **agg(tan_f)})
            rows.append({"module": m, "method": "random", "eval": "scaffold_cv", **agg(rnd_f)})
            for name in per_var:
                if per_var[name]:
                    rows.append({"module": m, "method": name, "eval": "scaffold_cv", **agg(per_var[name])})
        # 特征组重要性（LTR_full）
        if imp_acc is not None:
            n = 4 * dim
            imp_rows.append({"module": m, "imp_structure": float(imp_acc[:n].sum()),
                             "imp_ec": float(imp_acc[n:n+7].sum()), "imp_tanimoto": float(imp_acc[n+7:].sum())})
        # 时间留出（LTR_full）
        fut = clean[pd.to_numeric(clean["year"], errors="coerce") >= 2023]
        if fut["sb"].nunique() >= 2:
            tr = C[~C["sb"].isin(set(fut["sb"]))]
            Xtr, ok = build_feats(tr, emb, dim, fpc, True, True); trk = tr[ok]; qid = pd.factorize(trk["sb"])[0]
            model = train(Xtr[ok], trk["y"].to_numpy(np.float32), qid)
            def sf2(g, _m=model):
                Xte, oke = build_feats(g, emb, dim, fpc, True, True); sc = np.full(len(g), -1e9)
                if oke.any(): sc[oke] = _m.predict(Xte[oke])
                return sc
            r = within_sub_eval(sf2, C[C["sb"].isin(set(fut["sb"]))])
            rows.append({"module": m, "method": "LTR_full", "eval": "temporal>=2023", **r})
        # 终模型（全数据，LTR_full）
        Xall, ok = build_feats(C, emb, dim, fpc, True, True); pk = C[ok]; qid = pd.factorize(pk["sb"])[0]
        final = train(Xall[ok], pk["y"].to_numpy(np.float32), qid)
        final.save_model(str(models_dir / f"ltr_v3_{m}.json"))
        (models_dir / f"ltr_v3_{m}.meta.json").write_text(json.dumps(
            {"module": m, "params": BEST, "dim": dim, "feat": "chem+ec+tanimoto",
             "feat_dim": 4*dim+7+1}, indent=2), encoding="utf-8")

    res = pd.DataFrame(rows)
    kio.write_table(res, LTR_OUT / "v3_metrics.csv", force=True)
    imp = pd.DataFrame(imp_rows)
    if len(imp):
        kio.write_table(imp, LTR_OUT / "v3_feature_importance.csv", force=True)
    show = res[res["eval"] == "scaffold_cv"][["module", "method", "n", "mrr", "r@1", "r@5", "r@10"]]
    kio.log.info("== v3 same-substrate ranking (scaffold-CV mean±std) ==\n%s", show.to_string(index=False))
    print("\n", show.to_string(index=False))
    if len(imp):
        print("\n特征组重要性:\n", imp.round(3).to_string(index=False))
    tmp = res[res.eval.str.startswith("temporal")]
    if len(tmp):
        print("\n时间留出(新文献):\n", tmp[["module", "n", "mrr", "r@5", "r@10"]].round(3).to_string(index=False))


if __name__ == "__main__":
    main()
