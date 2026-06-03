"""
ltr_production.py — 把检索 LTR 推到"可报/可部署"：基线对照 + 调参 + per-fold CI + 存终模型

做四件事：
  1) 基线：random + Tanimoto(Morgan 相似度检索)  ——证明 LTR 比"挑最像底物的产物"强多少
  2) 调参：在 weak 数据上 grid（不碰干净测试），选 best XGBoost 参数
  3) 干净集 scaffold-CV：{random, tanimoto, LTR(tuned)} × recall@k/MRR，报 per-fold mean±std(CI)
  4) 训终模型：每模块用全部数据(weak+clean)+best 参数训练并保存 → 部署用

输出 modular/outputs/ltr/{production_metrics.csv, models/ltr_<module>.json + meta}
用法：python scripts/ssrf/ml_ranking_kernel/modular/ltr_production.py [--force]
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
from ltr_build import EV_COLS, LTR_OUT  # noqa: E402
from ltr_train_eval import encode_all, build_X, retrieval_eval, TOPK  # noqa: E402

RDLogger.DisableLog("rdApp.*")
PARAM_GRID = [
    dict(n_estimators=160, max_depth=5, learning_rate=0.12, colsample_bytree=0.5),
    dict(n_estimators=300, max_depth=6, learning_rate=0.08, colsample_bytree=0.5),
    dict(n_estimators=400, max_depth=4, learning_rate=0.10, colsample_bytree=0.6),
    dict(n_estimators=250, max_depth=8, learning_rate=0.06, colsample_bytree=0.4),
]


def fp(smiles):
    m = Chem.MolFromSmiles(smiles) if isinstance(smiles, str) else None
    return AllChem.GetMorganFingerprintAsBitVect(m, 2, nBits=2048) if m else None


def train_ranker_p(X, y, qid, params):
    import xgboost as xgb
    order = np.argsort(qid, kind="mergesort")
    m = xgb.XGBRanker(objective="rank:ndcg", subsample=0.8, min_child_weight=2,
                      n_jobs=12, tree_method="hist", **params)
    m.fit(X[order], y[order], qid=qid[order])
    return m


def tanimoto_eval(test_subs, vocab, sub_fp, prod_fp, true_prod):
    vB = vocab["product_block1"].tolist()
    crs = []
    for _, r in test_subs.iterrows():
        sf = sub_fp.get(r["substrate_smiles"]); truth = true_prod.get(r["sb"], set())
        if sf is None or not truth:
            continue
        sims = np.array([DataStructs.TanimotoSimilarity(sf, f) if f is not None else -1 for f in prod_fp])
        order = np.argsort(-sims, kind="mergesort")
        ranks = np.empty(len(sims)); ranks[order] = np.arange(1, len(sims) + 1)
        best = min((ranks[i] for i, b in enumerate(vB) if b in truth), default=None)
        if best:
            crs.append(best)
    cr = np.array(crs, float)
    return {"n": len(cr), "mrr": float(np.mean(1/cr)) if len(cr) else 0,
            **{f"r@{k}": float(np.mean(cr <= k)) if len(cr) else 0 for k in TOPK}}


def random_eval(test_subs, vocab, true_prod, rng):
    n = len(vocab); vB = vocab["product_block1"].tolist()
    crs = []
    for _, r in test_subs.iterrows():
        truth = true_prod.get(r["sb"], set())
        if not truth:
            continue
        order = rng.permutation(n)
        ranks = np.empty(n); ranks[order] = np.arange(1, n+1)
        best = min((ranks[i] for i, b in enumerate(vB) if b in truth), default=None)
        if best:
            crs.append(best)
    cr = np.array(crs, float)
    return {"n": len(cr), "mrr": float(np.mean(1/cr)) if len(cr) else 0,
            **{f"r@{k}": float(np.mean(cr <= k)) if len(cr) else 0 for k in TOPK}}


def tune(P_weak, vocab, true_prod, emb, dim, rng):
    """在 weak 上 train/val 选参（不碰干净测试）。"""
    subs = P_weak.loc[P_weak.y == 1, "sb"].unique()
    rng.shuffle(subs)
    val = set(subs[:min(150, len(subs)//5)])
    tr = P_weak[~P_weak.sb.isin(val)]
    Xtr, ok = build_X(tr, emb, dim, False); trk = tr[ok]
    qid = pd.factorize(trk["sb"])[0]
    val_df = P_weak[P_weak.sb.isin(val)].drop_duplicates("sb")[["sb", "substrate_smiles"]]
    best, best_r = PARAM_GRID[0], -1
    for params in PARAM_GRID:
        model = train_ranker_p(Xtr[ok], trk["y"].to_numpy(np.float32), qid, params)
        r = retrieval_eval(model, val_df, vocab, true_prod, emb, dim, False)
        kio.log.info("   tune %s -> weak-val r@10=%.3f", params, r["r@10"])
        if r["r@10"] > best_r:
            best_r, best = r["r@10"], params
    kio.log.info("  BEST params: %s (weak-val r@10=%.3f)", best, best_r)
    return best


def agg(rows):
    out = {}
    for k in ["mrr", "r@1", "r@5", "r@10", "r@20"]:
        v = [x[k] for x in rows]
        out[k] = f"{np.mean(v):.3f}±{np.std(v):.3f}"
    out["n"] = int(sum(x["n"] for x in rows))
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()
    kio.setup_logging()
    rng = np.random.default_rng(0)

    pairs = pd.read_parquet(LTR_OUT / "pairs.parquet")
    vocab_all = pd.read_parquet(LTR_OUT / "vocab.parquet")
    rxn = pd.read_parquet(LTR_OUT / "reactions.parquet")
    emb, dim = encode_all(set(pairs["substrate_smiles"]) | set(pairs["product_smiles"]) | set(vocab_all["product_smiles"]))

    rows = []
    models_dir = kio.safe_output_path("modular/ltr/models/_m").parent
    for m in ["A", "B", "C", "D"]:
        P = pairs[pairs["module"] == m].copy()
        V = vocab_all[vocab_all["module"] == m].copy()
        R = rxn[rxn["module"] == m]
        true_prod = R.groupby("sb")["pb"].agg(set).to_dict()
        # weak 子采样
        wsubs = P.loc[(P.scaffold_fold == -1) & (P.y == 1), "sb"].unique()
        if len(wsubs) > 2500:
            keep = set(rng.choice(wsubs, 2500, replace=False))
            P = P[(P.scaffold_fold >= 0) | (P.sb.isin(keep))].copy()
        clean = P[P.scaffold_fold >= 0].drop_duplicates("sb")
        folds = sorted(clean.scaffold_fold.unique())
        kio.log.info("== module %s: clean_subs=%d vocab=%d ==", m, clean.sb.nunique(), len(V))

        # FP 预算
        sub_fp = {s: fp(s) for s in clean["substrate_smiles"].unique()}
        prod_fp = [fp(s) for s in V["product_smiles"]]

        # 调参（weak）
        best = tune(P[P.scaffold_fold == -1], V, true_prod, emb, dim, np.random.default_rng(1))

        ltr_f, tan_f, rnd_f = [], [], []
        for f in folds:
            te = clean[clean.scaffold_fold == f][["sb", "substrate_smiles"]]
            tr = P[P.scaffold_fold != f]
            Xtr, ok = build_X(tr, emb, dim, False); trk = tr[ok]
            if (trk.y == 1).sum() < 3 or len(te) == 0:
                continue
            qid = pd.factorize(trk["sb"])[0]
            model = train_ranker_p(Xtr[ok], trk["y"].to_numpy(np.float32), qid, best)
            ltr_f.append(retrieval_eval(model, te, V, true_prod, emb, dim, False))
            tan_f.append(tanimoto_eval(te, V, sub_fp, prod_fp, true_prod))
            rnd_f.append(random_eval(te, V, true_prod, np.random.default_rng(7)))
        for name, fr in [("LTR_tuned", ltr_f), ("tanimoto", tan_f), ("random", rnd_f)]:
            if fr:
                rows.append({"module": m, "method": name, "eval": "scaffold_cv_perfold", **agg(fr)})

        # 时间留出
        future = clean[pd.to_numeric(clean.year if "year" in clean else np.nan, errors="coerce") >= 2023] \
            if "year" in clean.columns else clean.iloc[0:0]
        # 终模型：全部数据训练 + 保存
        Xall, ok = build_X(P, emb, dim, False); pk = P[ok]
        qid = pd.factorize(pk["sb"])[0]
        final = train_ranker_p(Xall[ok], pk["y"].to_numpy(np.float32), qid, best)
        final.save_model(str(models_dir / f"ltr_{m}.json"))
        (models_dir / f"ltr_{m}.meta.json").write_text(json.dumps(
            {"module": m, "params": best, "dim": dim, "feat": "chem_only_reaction_diff",
             "vocab_size": len(V), "n_train_pairs": int(len(pk))}, indent=2), encoding="utf-8")

    res = pd.DataFrame(rows)
    kio.write_table(res, LTR_OUT / "production_metrics.csv", force=True)
    show = res[["module", "method", "n", "mrr", "r@1", "r@5", "r@10", "r@20"]]
    kio.log.info("== PRODUCTION metrics (mean±std over folds) ==\n%s", show.to_string(index=False))
    print("\n", show.to_string(index=False))


if __name__ == "__main__":
    main()
