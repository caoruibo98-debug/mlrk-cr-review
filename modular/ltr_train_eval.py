"""
ltr_train_eval.py — 多证据 learning-to-rank：训练 + 诚实评估

模型：XGBoost LambdaMART（rank:ndcg），每模块独立。
任务：产物检索——给底物 S，在模块产物词表里把真产物排前面。

诚实设计（关键）：
  - 证据特征(菌/酶/文献)是"答案泄漏"——只有已知反应有，新颖候选检索时为 0。
    所以做两套：chem-only(可部署，能泛化) vs chem+evidence(展示 in-sample 拟合)。
    **新文献召回的 headline 用 chem-only。** 证据当置信度/可解释 overlay。
评估：
  1) scaffold-CV 检索 recall@k / MRR（干净集 OOF）
  2) 时间留出：train<2023 → test≥2023 = 新文献召回
  3) 证据 overlay：top 候选里多少有已知 菌/酶/文献
  4) 特征组重要性

用法：python scripts/ssrf/ml_ranking_kernel/modular/ltr_train_eval.py [--force]
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import kio  # noqa: E402
from ltr_build import EV_COLS, LTR_OUT  # noqa: E402

TOPK = [1, 3, 5, 10, 20]


def encode_all(smiles: set[str]):
    """ChemBERTa 编码所有 SMILES（缓存复用 v1 mol_embeddings + 增量）。"""
    import importlib
    enc_mod = importlib.import_module("encode_molecules")
    cache = kio.FEATURES_DIR / "mol_embeddings.parquet"
    emb = {}
    dim = 384
    if cache.exists():
        c = pd.read_parquet(cache)
        emb = {s: np.asarray(v, np.float32) for s, v in zip(c["smiles"], c["embedding"])}
        dim = len(next(iter(emb.values())))
    need = sorted(s for s in smiles if s and s not in emb)
    if need:
        kio.log.info("encoding %d new SMILES ...", len(need))
        cfg = kio.load_config()
        fe = enc_mod.FrozenEncoder(cfg["encoder"]["primary"], device=cfg["train"].get("device", "cuda"))
        vecs = fe.encode(need)
        dim = vecs.shape[1]
        for s, v in zip(need, vecs):
            emb[s] = np.asarray(v, np.float32)
    return emb, dim


def _l2(v):
    return v / max(float(np.linalg.norm(v)), 1e-9)


def rxn_feats(subs, prods, emb, dim):
    X = np.zeros((len(subs), 4 * dim), np.float32)
    ok = np.zeros(len(subs), bool)
    for i, (s, p) in enumerate(zip(subs, prods)):
        zs, zp = emb.get(s), emb.get(p)
        if zs is None or zp is None:
            continue
        zs, zp = _l2(zs), _l2(zp)
        X[i] = np.concatenate([zs, zp, zs - zp, zs * zp])
        ok[i] = True
    return X, ok


def build_X(df, emb, dim, use_evidence):
    X, ok = rxn_feats(df["substrate_smiles"].tolist(), df["product_smiles"].tolist(), emb, dim)
    if use_evidence:
        ev = df[EV_COLS].to_numpy(np.float32)
        X = np.hstack([X, ev])
    return X, ok


def train_ranker(X, y, qid):
    import xgboost as xgb
    order = np.argsort(qid, kind="mergesort")
    m = xgb.XGBRanker(objective="rank:ndcg", n_estimators=160, max_depth=5,
                      learning_rate=0.12, subsample=0.8, colsample_bytree=0.5,
                      min_child_weight=2, n_jobs=12, tree_method="hist")
    m.fit(X[order], y[order], qid=qid[order])
    return m


def retrieval_eval(model, test_subs_df, vocab, true_prod, emb, dim, use_evidence):
    """每测试底物：在模块词表打分，算真产物 recall@k / MRR。"""
    vP = vocab["product_smiles"].tolist()
    vB = vocab["product_block1"].tolist()
    case_ranks = []
    for _, r in test_subs_df.iterrows():
        s = r["substrate_smiles"]; sb = r["sb"]
        truth = true_prod.get(sb, set())
        if not truth:
            continue
        df = pd.DataFrame({"substrate_smiles": [s] * len(vP), "product_smiles": vP})
        for c in EV_COLS:
            df[c] = 0.0  # 检索时新颖候选无证据（答案泄漏不可用）
        X, ok = build_X(df, emb, dim, use_evidence)
        scores = np.full(len(vP), -1e9)
        if ok.any():
            scores[ok] = model.predict(X[ok])
        order = np.argsort(-scores, kind="mergesort")
        ranks = np.empty(len(scores)); ranks[order] = np.arange(1, len(scores) + 1)
        best = min((ranks[i] for i, b in enumerate(vB) if b in truth), default=None)
        if best is not None:
            case_ranks.append(best)
    cr = np.array(case_ranks, float)
    out = {"n": len(cr), "mrr": float(np.mean(1 / cr)) if len(cr) else 0.0}
    for k in TOPK:
        out[f"r@{k}"] = float(np.mean(cr <= k)) if len(cr) else 0.0
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--modules", nargs="+", default=["A", "B", "C", "D"])
    args = ap.parse_args()
    kio.setup_logging()

    pairs = pd.read_parquet(LTR_OUT / "pairs.parquet")
    vocab = pd.read_parquet(LTR_OUT / "vocab.parquet")
    rxn = pd.read_parquet(LTR_OUT / "reactions.parquet")
    emb, dim = encode_all(set(pairs["substrate_smiles"]) | set(pairs["product_smiles"]) | set(vocab["product_smiles"]))
    kio.log.info("emb dim=%d; pairs=%d", dim, len(pairs))

    MAX_WEAK = 2500
    rng = np.random.default_rng(0)
    results = []
    for m in args.modules:
        P = pairs[pairs["module"] == m].copy()
        # 控规模：weak 底物(scaffold_fold==-1) 子采样到 MAX_WEAK，干净底物全留
        weak_subs = P.loc[(P["scaffold_fold"] == -1) & (P["y"] == 1), "sb"].unique()
        if len(weak_subs) > MAX_WEAK:
            keep = set(rng.choice(weak_subs, MAX_WEAK, replace=False))
            P = P[(P["scaffold_fold"] >= 0) | (P["sb"].isin(keep))].copy()
        V = vocab[vocab["module"] == m].copy()
        R = rxn[rxn["module"] == m].copy()
        true_prod = R.groupby("sb")["pb"].agg(set).to_dict()
        clean = P[P["scaffold_fold"] >= 0].drop_duplicates("sb")
        folds = sorted(clean["scaffold_fold"].unique())
        kio.log.info("== module %s: pairs=%d vocab=%d clean_subs=%d folds=%s ==",
                     m, len(P), len(V), clean["sb"].nunique(), folds)

        for use_ev in [False, True]:
            tag = "chem+ev" if use_ev else "chem-only"
            # scaffold-CV 检索
            oof = []
            for f in folds:
                tr = P[(P["scaffold_fold"] != f)]
                Xtr, ok = build_X(tr, emb, dim, use_ev)
                trk = tr[ok]
                if (trk["y"] == 1).sum() < 3:
                    continue
                qid = pd.factorize(trk["sb"])[0]
                model = train_ranker(Xtr[ok], trk["y"].to_numpy(np.float32), qid)
                te_subs = clean[clean["scaffold_fold"] == f][["sb", "substrate_smiles"]]
                oof.append(retrieval_eval(model, te_subs, V, true_prod, emb, dim, use_ev))
            if oof:
                agg = {k: float(np.mean([o[k] for o in oof])) for k in oof[0] if k != "n"}
                agg["n"] = int(sum(o["n"] for o in oof))
                results.append({"module": m, "feat": tag, "eval": "scaffold_cv", **agg})

        # 时间留出（chem-only headline）：train<2023 + 全weak → test 干净 year>=2023
        future = clean[pd.to_numeric(clean["year"], errors="coerce") >= 2023]
        if len(future):
            tr = P[~P["sb"].isin(set(future["sb"]))]
            Xtr, ok = build_X(tr, emb, dim, False)
            trk = tr[ok]; qid = pd.factorize(trk["sb"])[0]
            model = train_ranker(Xtr[ok], trk["y"].to_numpy(np.float32), qid)
            r = retrieval_eval(model, future[["sb", "substrate_smiles"]], V, true_prod, emb, dim, False)
            results.append({"module": m, "feat": "chem-only", "eval": "temporal>=2023", **r})

    res = pd.DataFrame(results)
    kio.write_table(res, LTR_OUT / "ltr_metrics.csv", force=True)
    show = res[["module", "feat", "eval", "n", "mrr", "r@1", "r@5", "r@10", "r@20"]].round(3)
    kio.log.info("== LTR retrieval metrics ==\n%s", show.to_string(index=False))
    print("\n", show.to_string(index=False))


if __name__ == "__main__":
    main()
