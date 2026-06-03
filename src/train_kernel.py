"""
train_kernel.py — Phase 4

规则无关、RDKit-free 打分的神经反应排序内核。
特征 = 冻结化学 LM 的反应表示 [z_s,z_p,z_s-z_p,z_s*z_p]（rxnfeat），**不含任何规则字段**。

损失：
  - within-substrate pairwise margin（真产物 > 同底物 decoy）  ← 直接优化排序
  - BCE（pos=1 / decoy=0，样本权重）
  - nnPU（pos + unlabeled，非负风险校正）                      ← 利用大未标注池

CV：scaffold split + random split 各做 GroupKFold；输出 **out-of-fold** 打分的候选全集
（每个测试底物：真产物 ∪ 该底物全部规则候选），供 evaluate_kernel 算 TopK/MRR/AUROC。

用法：
  python scripts/ssrf/ml_ranking_kernel/src/train_kernel.py [--smoke] [--force]
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import kio
import rxnfeat


def set_seed(seed: int):
    import torch
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)


def make_mlp(in_dim: int, hidden: list[int], dropout: float):
    import torch.nn as nn
    layers = []
    d = in_dim
    for h in hidden:
        layers += [nn.Linear(d, h), nn.ReLU(), nn.Dropout(dropout)]
        d = h
    layers += [nn.Linear(d, 1)]
    return nn.Sequential(*layers)


def nnpu_risk(g_pos, g_unl, prior: float):
    import torch
    # sigmoid 替代损失 l(z)=sigmoid(-z)（把 z 判为正的代价）
    l_pos_pos = torch.sigmoid(-g_pos).mean()
    l_pos_neg = torch.sigmoid(g_pos).mean()
    l_unl_neg = torch.sigmoid(g_unl).mean()
    neg = l_unl_neg - prior * l_pos_neg
    return prior * l_pos_pos + torch.clamp(neg, min=0.0)   # 非负风险（nnPU）


def fit_fold(Xtr, ytr, wtr, grp_tr, X_unl, tr_cfg, in_dim, seed):
    import torch
    import torch.nn as nn
    set_seed(seed)
    dev = tr_cfg.get("device", "cuda")
    dev = dev if torch.cuda.is_available() else "cpu"
    model = make_mlp(in_dim, tr_cfg["hidden"], tr_cfg["dropout"]).to(dev)
    opt = torch.optim.Adam(model.parameters(), lr=tr_cfg["lr"], weight_decay=tr_cfg["weight_decay"])
    bce = nn.BCEWithLogitsLoss(reduction="none")

    Xt = torch.tensor(Xtr, device=dev)
    yt = torch.tensor(ytr, dtype=torch.float32, device=dev)
    wt = torch.tensor(wtr, dtype=torch.float32, device=dev)
    Xu = torch.tensor(X_unl, device=dev) if len(X_unl) else None

    # within-substrate pos/neg 对
    pairs_pos, pairs_neg = [], []
    grp = np.asarray(grp_tr)
    for g in np.unique(grp):
        idx = np.where(grp == g)[0]
        pos_i = idx[ytr[idx] == 1]
        neg_i = idx[ytr[idx] == 0]
        for pi in pos_i:
            for ni in neg_i:
                pairs_pos.append(pi); pairs_neg.append(ni)
    has_pairs = len(pairs_pos) > 0
    if has_pairs:
        pp = torch.tensor(np.array(pairs_pos), device=dev)
        pn = torch.tensor(np.array(pairs_neg), device=dev)
    pos_mask = torch.tensor(ytr == 1, device=dev)

    margin = tr_cfg["pairwise_margin"]
    pw = tr_cfg["pairwise_weight"]
    puw = tr_cfg["pu_weight"]
    prior = tr_cfg["pu_prior"]
    rng = np.random.default_rng(seed)

    best_state, best_loss, bad = None, float("inf"), 0
    for ep in range(tr_cfg["epochs"]):
        model.train()
        opt.zero_grad()
        g_all = model(Xt).squeeze(-1)
        loss = (bce(g_all, yt) * wt).sum() / wt.sum().clamp(min=1e-9)
        if has_pairs:
            diff = g_all[pp] - g_all[pn]
            loss = loss + pw * torch.clamp(margin - diff, min=0.0).mean()
        if Xu is not None and puw > 0 and pos_mask.any():
            bidx = rng.integers(0, Xu.shape[0], size=min(2048, Xu.shape[0]))
            g_unl = model(Xu[torch.tensor(bidx, device=dev)]).squeeze(-1)
            loss = loss + puw * nnpu_risk(g_all[pos_mask], g_unl, prior)
        loss.backward()
        opt.step()

        lv = float(loss.detach())
        if lv < best_loss - 1e-4:
            best_loss, bad = lv, 0
            best_state = {k: v.detach().cpu().clone() for k, v in model.state_dict().items()}
        else:
            bad += 1
            if bad >= tr_cfg["patience"]:
                break
    if best_state is not None:
        model.load_state_dict(best_state)
    model.eval()
    return model


def score_model(model, X, dev):
    import torch
    if len(X) == 0:
        return np.zeros(0, dtype=np.float32)
    with torch.no_grad():
        g = model(torch.tensor(X, device=dev)).squeeze(-1)
        return torch.sigmoid(g).float().cpu().numpy()


def build_eval_universe(pos: pd.DataFrame, cand: pd.DataFrame, sub_ids: set) -> pd.DataFrame:
    """某些测试底物的候选全集：真产物 ∪ 该底物全部规则候选；y_eval=是否真产物。"""
    rows = []
    pos_s = pos[pos["substrate_block1"].isin(sub_ids)]
    cand_s = cand[cand["substrate_block1"].isin(sub_ids)]
    true_prod = pos_s.groupby("substrate_block1")["product_block1"].agg(set).to_dict()
    sub_smiles = pos_s.drop_duplicates("substrate_block1").set_index("substrate_block1")["substrate_smiles"].to_dict()
    cat = pos_s.drop_duplicates("substrate_block1").set_index("substrate_block1")["reaction_category"].to_dict()
    leak = pos_s.drop_duplicates("substrate_block1").set_index("substrate_block1")["leakage_flag"].to_dict()

    rclss_lookup = cand_s.set_index(["substrate_block1", "product_block1"])["rclss_prior"].to_dict() \
        if "rclss_prior" in cand_s.columns else {}

    for sub in sub_ids:
        seen = set()
        # 候选（规则生成）
        for _, r in cand_s[cand_s["substrate_block1"].eq(sub)].iterrows():
            pb = r["product_block1"]
            if pb in seen:
                continue
            seen.add(pb)
            rows.append({"substrate_block1": sub, "substrate_smiles": sub_smiles.get(sub, r.get("substrate_smiles")),
                         "product_smiles": r["product_smiles"], "product_block1": pb,
                         "rclss_prior": r.get("rclss_prior", np.nan),
                         "y_eval": int(pb in true_prod.get(sub, set())),
                         "reaction_category": cat.get(sub, ""), "leakage_flag": leak.get(sub, "")})
        # 真产物若规则没生成（generation miss）也要加入 → ml 可打分，rclss=NaN 排末
        for _, r in pos_s[pos_s["substrate_block1"].eq(sub)].iterrows():
            pb = r["product_block1"]
            if pb in seen:
                continue
            seen.add(pb)
            rows.append({"substrate_block1": sub, "substrate_smiles": r["substrate_smiles"],
                         "product_smiles": r["product_smiles"], "product_block1": pb,
                         "rclss_prior": rclss_lookup.get((sub, pb), np.nan),
                         "y_eval": 1, "reaction_category": cat.get(sub, ""), "leakage_flag": leak.get(sub, "")})
    return pd.DataFrame(rows)


def run_split(split_col: str, split_name: str, edges: pd.DataFrame, cand: pd.DataFrame,
              emb, dim, tr_cfg, seed, mode: str = "full") -> pd.DataFrame:
    import torch
    dev = tr_cfg.get("device", "cuda")
    dev = dev if torch.cuda.is_available() else "cpu"
    in_dim = rxnfeat.mode_ndim(mode, dim)
    l2 = tr_cfg.get("l2_normalize", True)
    labeled = edges[edges["y"].isin([0, 1])].copy()
    pu = edges[edges["y"] == -1].copy()
    pos = labeled[labeled["y"] == 1].copy()

    folds = sorted(f for f in labeled[split_col].unique() if f >= 0)
    kio.log.info("[%s|mode=%s] folds=%s labeled=%d pu=%d", split_name, mode, folds, len(labeled), len(pu))

    # PU 特征（一次）+ 保留其 fold 标签（防泄漏：PU 行也属于 71 个正底物，带 fold）
    Xu_all, oku = rxnfeat.build_features(pu["substrate_smiles"], pu["product_smiles"], emb, dim, l2, mode)
    Xu_all = Xu_all[oku]
    pu_fold = pu[split_col].to_numpy()[oku]   # 关键：按 split 取 PU 的 fold

    oof = []
    for f in folds:
        tr = labeled[labeled[split_col] != f]
        te_subs = set(labeled.loc[labeled[split_col] == f, "substrate_block1"])
        Xtr, okt = rxnfeat.build_features(tr["substrate_smiles"], tr["product_smiles"], emb, dim, l2, mode)
        trk = tr[okt]
        # 泄漏修复：排除测试 fold 底物的 PU 边（否则模型在无标注里见过测试底物候选）
        Xu_f = Xu_all[pu_fold != f]
        model = fit_fold(Xtr[okt], trk["y"].to_numpy().astype(np.float32),
                         trk["weight"].to_numpy().astype(np.float32),
                         trk["substrate_block1"].to_numpy(), Xu_f, tr_cfg, in_dim, seed + f)
        # OOF：测试底物候选全集打分
        uni = build_eval_universe(pos, cand, te_subs)
        if uni.empty:
            continue
        Xe, oke = rxnfeat.build_features(uni["substrate_smiles"], uni["product_smiles"], emb, dim, l2, mode)
        uni = uni[oke].copy()
        uni["ml_kernel_score"] = score_model(model, Xe[oke], dev)
        uni["fold"] = f
        uni["split"] = split_name
        oof.append(uni)
    return pd.concat(oof, ignore_index=True) if oof else pd.DataFrame()


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--smoke", action="store_true", help="少 epoch 快跑")
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()
    kio.setup_logging()
    cfg = kio.load_config()
    tr_cfg = dict(cfg["train"])
    tr_cfg["l2_normalize"] = cfg["reaction_features"].get("l2_normalize", True)
    seed = int(cfg["splits"]["random_seed"])
    if args.smoke:
        tr_cfg["epochs"] = 30

    edges = pd.read_parquet(kio.SPLITS_DIR / "edges_with_splits.parquet")
    cand_p = kio.DATASETS_DIR / "substrate_candidates.parquet"
    cand = pd.read_parquet(cand_p).rename(columns={"rclss_production_score": "rclss_prior"}) \
        if cand_p.exists() else edges[edges["y"] == 0].assign(rclss_prior=edges.get("rclss_prior"))
    cand["substrate_block1"] = cand["substrate_block1"].astype(str).str.upper()
    cand["product_block1"] = cand["product_block1"].astype(str).str.upper()

    emb, dim = rxnfeat.load_embeddings()
    kio.log.info("embeddings: %d mols, dim=%d -> reaction feat dim=%d", len(emb), dim, 4 * dim)

    # 规则无关断言：特征列全是 emb 衍生
    feat_names = rxnfeat.feature_block_names(dim)
    banned = ("rule", "ec", "rclss", "category", "uid", "smarts")
    assert not any(any(b in fn.lower() for b in banned) for fn in feat_names), "rule-derived feature leaked!"
    kio.log.info("[assert] rule-agnostic feature schema OK (%d emb-derived dims)", len(feat_names))

    for split_col, split_name in [("scaffold_fold", "scaffold"), ("random_fold", "random")]:
        oof = run_split(split_col, split_name, edges, cand, emb, dim, tr_cfg, seed)
        kio.write_table(oof, kio.PREDICTIONS_DIR / f"oof_scored_edges.{split_name}.parquet", force=args.force)
        if len(oof):
            kio.log.info("[%s] OOF scored: %d edges over %d substrates (%d positives)",
                         split_name, len(oof), oof["substrate_block1"].nunique(), int((oof.y_eval == 1).sum()))

    # ===== 最终模型：在全部 labeled 上训练，供 predict_kernel 用 =====
    import torch
    labeled = edges[edges["y"].isin([0, 1])].copy()
    pu = edges[edges["y"] == -1]
    Xall, oka = rxnfeat.build_features(labeled["substrate_smiles"], labeled["product_smiles"], emb, dim,
                                       tr_cfg.get("l2_normalize", True))
    lab = labeled[oka]
    Xu, oku = rxnfeat.build_features(pu["substrate_smiles"], pu["product_smiles"], emb, dim,
                                     tr_cfg.get("l2_normalize", True))
    final = fit_fold(Xall[oka], lab["y"].to_numpy().astype(np.float32),
                     lab["weight"].to_numpy().astype(np.float32),
                     lab["substrate_block1"].to_numpy(), Xu[oku], tr_cfg, 4 * dim, seed)
    torch.save({"state_dict": final.state_dict(), "in_dim": 4 * dim, "hidden": tr_cfg["hidden"],
                "dropout": tr_cfg["dropout"], "dim": dim, "l2_normalize": tr_cfg.get("l2_normalize", True)},
               kio.safe_output_path("models/kernel_final.pt"))
    kio.log.info("saved final model (trained on %d labeled edges)", len(lab))

    prov = {"feature_source": "frozen_chem_LM reaction features [z_s,z_p,z_s-z_p,z_s*z_p]",
            "n_feature_dims": len(feat_names), "rule_fields_used": [],
            "encoder_meta": json.loads((kio.FEATURES_DIR / "encoder_meta.json").read_text())
            if (kio.FEATURES_DIR / "encoder_meta.json").exists() else {}}
    kio.safe_output_path("models/feature_provenance.json").write_text(
        json.dumps(prov, indent=2, ensure_ascii=False), encoding="utf-8")
    kio.log.info("wrote feature_provenance.json")


if __name__ == "__main__":
    main()
