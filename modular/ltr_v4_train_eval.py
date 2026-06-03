"""
ltr_v4_train_eval.py — 去污染重测（A 真实decoy + B 产物不相交切分）

针对 v3 两个污染：
  A. RunReactants decoy 是"怪结构" → 模型靠"像不像真代谢物"取巧（非变换）。
     修：(1) QED 过滤掉低质量 decoy；(2) 加"真实代谢物近邻 decoy"——词表里离**真产物**最近的真代谢物
        → 真/假都是真实代谢物且结构相近 → "可信度/相似度"都失效 → 逼模型用底物→产物变换。
  B. scaffold split 只留出底物、产物会复现 → 产物记忆。
     修：训练正样本里**剔除所有测试折真产物**（product-disjoint 正样本）。

对比 {tanimoto, random, LTR_chem, LTR_full(chem+ec+tan)}；报 per-fold mean±std + 特征组重要性。
诚实预期：若去污染后 LTR 仍 > Tanimoto 且 EC 重要性上来 → 真有变换信号；若 LTR≈Tanimoto → 之前的胜利是 decoy 取巧。

用法：python scripts/ssrf/ml_ranking_kernel/modular/ltr_v4_train_eval.py [--force]
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs, QED
from rdkit import RDLogger
from sklearn.model_selection import GroupKFold

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import kio  # noqa: E402
from ltr_build import LTR_OUT  # noqa: E402
from ltr_train_eval import encode_all, rxn_feats, TOPK  # noqa: E402
from ltr_v3_train_eval import fp, build_feats, train, within_sub_eval, agg, BEST  # noqa: E402

RDLogger.DisableLog("rdApp.*")
QED_MIN = 0.30          # 过滤 RunReactants junk
N_VOCAB_NN = 12         # 每底物加多少"真产物近邻"真代谢物 decoy
MAX_DECOY = 40


def qed_of(smiles, cache):
    if smiles in cache:
        return cache[smiles]
    m = Chem.MolFromSmiles(smiles) if isinstance(smiles, str) else None
    try:
        q = QED.qed(m) if m else 0.0
    except Exception:
        q = 0.0
    cache[smiles] = q
    return q


def build_realistic_candidates(cand, vocab, fpc, rng):
    """A：QED 过滤 RunReactants decoy + 加真产物近邻真代谢物 decoy。"""
    qcache = {}
    out_rows = []
    for m in ["A", "B", "C", "D"]:
        C = cand[cand.module == m]
        Vp = vocab[vocab.module == m][["product_block1", "product_smiles"]].drop_duplicates()
        vb = Vp["product_block1"].tolist(); vs = Vp["product_smiles"].tolist()
        vfp = [fpc.get(s) for s in vs]
        true_by_sub = C[C.y == 1].groupby("sb")["pb"].agg(set).to_dict()
        for sb, g in C.groupby("sb"):
            s = g["substrate_smiles"].iloc[0]; scaf = g["substrate_scaffold"].iloc[0]
            is_clean = bool(g["is_clean"].iloc[0]); yr = g["year"].iloc[0]
            tier = g["tier"].iloc[0] if "tier" in g.columns else ("gold" if is_clean else "weak")
            truth = true_by_sub.get(sb, set())
            pos = g[g.y == 1].drop_duplicates("pb")
            # 真实 decoy 来源1：QED 过滤过的同底物 RunReactants decoy
            dd = g[(g.y == 0)].drop_duplicates("pb")
            dd = dd[dd["product_smiles"].map(lambda x: qed_of(x, qcache)) >= QED_MIN]
            decoys = [(r.pb, r.product_smiles, int(r.ec_class)) for r in dd.itertuples()]
            # 真实 decoy 来源2：词表里离"真产物"最近的真代谢物（near-miss，真/假都真实）
            for _, pr in pos.iterrows():
                tfp = fpc.get(pr["product_smiles"])
                if tfp is None:
                    continue
                sims = np.array([DataStructs.TanimotoSimilarity(tfp, f) if f is not None else -1 for f in vfp])
                for j in np.argsort(-sims)[:N_VOCAB_NN * 2]:
                    if vb[j] in truth:
                        continue
                    decoys.append((vb[j], vs[j], 0))
            # 去重 + 限量
            seen = set(); uniq = []
            for pb, psm, ecc in decoys:
                if pb in truth or pb in seen:
                    continue
                seen.add(pb); uniq.append((pb, psm, ecc))
            if len(uniq) > MAX_DECOY:
                idx = rng.choice(len(uniq), MAX_DECOY, replace=False); uniq = [uniq[i] for i in idx]
            if not uniq or pos.empty:
                continue
            for _, pr in pos.iterrows():
                out_rows.append({"module": m, "sb": sb, "substrate_smiles": s, "substrate_scaffold": scaf,
                                 "pb": pr["pb"], "product_smiles": pr["product_smiles"],
                                 "ec_class": int(pr["ec_class"]), "y": 1, "is_clean": is_clean,
                                 "tier": tier, "year": yr})
            for pb, psm, ecc in uniq:
                out_rows.append({"module": m, "sb": sb, "substrate_smiles": s, "substrate_scaffold": scaf,
                                 "pb": pb, "product_smiles": psm, "ec_class": ecc, "y": 0,
                                 "is_clean": is_clean, "tier": tier, "year": yr})
    return pd.DataFrame(out_rows)


def main() -> None:
    ap = argparse.ArgumentParser(); ap.add_argument("--force", action="store_true"); args = ap.parse_args()
    kio.setup_logging()
    rng = np.random.default_rng(0)
    cand = pd.read_parquet(LTR_OUT / "candidates.parquet")
    vocab = pd.read_parquet(LTR_OUT / "vocab.parquet")
    allsm = set(cand.substrate_smiles) | set(cand.product_smiles) | set(vocab.product_smiles)
    fpc = {s: fp(s) for s in allsm}
    emb, dim = encode_all(allsm)

    C = build_realistic_candidates(cand, vocab, fpc, rng)
    C["scaf"] = C["substrate_scaffold"].fillna("NS_" + C["sb"])
    kio.log.info("realistic candidates=%d subs=%d (median/sub=%.0f)",
                 len(C), C.sb.nunique(), C.groupby("sb").size().median())

    rows = []; imp_rows = []
    VARIANTS = [("LTR_full", True, True), ("LTR_chem", False, False)]
    for m in ["A", "B", "C", "D"]:
        Cm = C[C.module == m].copy()
        clean = Cm[Cm.is_clean].drop_duplicates("sb").reset_index(drop=True)
        k = min(5, clean["scaf"].nunique())
        if k < 2 or len(clean) < 5:
            kio.log.info("module %s: too few clean subs (%d) skip", m, len(clean)); continue
        clean["fold"] = -1
        for f, (_, te) in enumerate(GroupKFold(k).split(clean, groups=clean["scaf"])):
            clean.iloc[te, clean.columns.get_loc("fold")] = f
        fold_of = dict(zip(clean.sb, clean.fold)); Cm["fold"] = Cm.sb.map(fold_of).fillna(-1).astype(int)
        kio.log.info("== module %s: clean_subs=%d folds=%d median_cand=%.0f ==",
                     m, len(clean), k, Cm[Cm.is_clean].groupby("sb").size().median())
        tan_f, rnd_f = [], []; per_var = {v[0]: [] for v in VARIANTS}; imp_acc = None
        for f in range(k):
            test = Cm[(Cm.fold == f) & (Cm.is_clean)]
            test_true = set(test.loc[test.y == 1, "pb"])
            # B: 产物不相交 —— 训练里剔除所有测试真产物的正样本
            tr = Cm[Cm.fold != f]
            tr = tr[~((tr.y == 1) & (tr.pb.isin(test_true)))]
            if test.sb.nunique() == 0 or (tr.y == 1).sum() < 5:
                continue
            tan_f.append(within_sub_eval(
                lambda g: np.array([DataStructs.TanimotoSimilarity(fpc.get(g.substrate_smiles.iloc[0]), fpc.get(p))
                                    if (fpc.get(g.substrate_smiles.iloc[0]) and fpc.get(p)) else 0.0
                                    for p in g["product_smiles"]]), test))
            r2 = np.random.default_rng(f); rnd_f.append(within_sub_eval(lambda g: r2.random(len(g)), test))
            for name, ue, ut in VARIANTS:
                Xtr, ok = build_feats(tr, emb, dim, fpc, ue, ut); trk = tr[ok]; qid = pd.factorize(trk.sb)[0]
                model = train(Xtr[ok], trk.y.to_numpy(np.float32), qid)
                def sf(g, _m=model, _ue=ue, _ut=ut):
                    Xt, okt = build_feats(g, emb, dim, fpc, _ue, _ut); sc = np.full(len(g), -1e9)
                    if okt.any(): sc[okt] = _m.predict(Xt[okt])
                    return sc
                per_var[name].append(within_sub_eval(sf, test))
                if name == "LTR_full" and imp_acc is None:
                    imp_acc = model.feature_importances_
        if tan_f:
            rows.append({"module": m, "method": "tanimoto", **agg(tan_f)})
            rows.append({"module": m, "method": "random", **agg(rnd_f)})
            for name in per_var:
                if per_var[name]:
                    rows.append({"module": m, "method": name, **agg(per_var[name])})
        if imp_acc is not None:
            n = 4 * dim
            imp_rows.append({"module": m, "structure": round(float(imp_acc[:n].sum()), 3),
                             "ec": round(float(imp_acc[n:n+7].sum()), 3),
                             "tanimoto": round(float(imp_acc[n+7:].sum()), 3)})
    res = pd.DataFrame(rows); kio.write_table(res, LTR_OUT / "v4_metrics.csv", force=True)
    imp = pd.DataFrame(imp_rows)
    if len(imp): kio.write_table(imp, LTR_OUT / "v4_feature_importance.csv", force=True)
    print("\n=== v4 去污染（真实decoy + 产物不相交）scaffold-CV mean±std ===")
    print(res[["module", "method", "n", "mrr", "r@1", "r@5", "r@10"]].to_string(index=False))
    if len(imp):
        print("\n特征组重要性(LTR_full):\n", imp.to_string(index=False))


if __name__ == "__main__":
    main()
