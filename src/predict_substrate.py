"""
predict_substrate.py — 端到端 CLI：输入一个物质(SMILES) → 预测产物排序

完整链路（非破坏，全部读 production，只写本包 outputs/predictions/）：
  输入 SMILES
   → 规则引擎生成候选 (SMARTS 匹配 + RCLSS + RunReactants)   [慢：引擎加载~20s + 打分]
   → 冻结 ChemBERTa 即时编码 substrate/product               [缺的才编码，已缓存的复用]
   → L-RCLSS 内核打分 (ml_kernel_score)
   → 排序 + 置信度 + 适用域(AD)警告 + 溯源链 JSON

⚠️ 诚实边界（见 analyze_report）：ml_kernel_score 是"候选可信度/初筛"分，主要由产物像不像
真实代谢物驱动（decoy bias），不是湿实验概率，也不保证建模了具体变换。AD 低 = 训练域外，别信。

用法：
  python scripts/ssrf/ml_ranking_kernel/src/predict_substrate.py --smiles "O=c1c(O)c(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12" --name quercetin --topn 15
  python scripts/ssrf/ml_ranking_kernel/src/predict_substrate.py --smiles "<SMILES>" --no-cache   # 强制走引擎
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import kio
import rxnfeat
from predict_kernel import load_final_model, confidence_tier
from train_kernel import score_model


# ---------- 候选生成 ----------

def candidates_for(smiles: str, name: str, use_cache: bool) -> tuple[pd.DataFrame, str]:
    block1 = kio.inchikey_block1(kio.smiles_to_inchikey(smiles))
    cache_p = kio.DATASETS_DIR / "substrate_candidates.parquet"
    if use_cache and cache_p.exists():
        cache = pd.read_parquet(cache_p).rename(columns={"rclss_production_score": "rclss_prior"})
        hit = cache[cache["substrate_block1"].astype(str).str.upper() == block1]
        if len(hit):
            kio.log.info("substrate %s 命中候选缓存 (%d 候选) — 跳过引擎", block1, len(hit))
            return hit.copy(), "cached"
    # 走 production 规则引擎（慢）
    kio.log.info("substrate %s 不在缓存 → 加载规则引擎生成候选（约 1-2 分钟）...", block1)
    sys.path.insert(0, str(kio.PROJECT_ROOT / "scripts" / "ssrf"))
    from rclss_similarity import RCLSSEngine
    from ssrf_filter import filter_compound
    from generate_candidates import generate_for_substrate
    t0 = time.time()
    eng = RCLSSEngine.load_default()
    kio.log.info("  引擎加载 %.1fs，开始打分...", time.time() - t0)
    cand = generate_for_substrate(eng, filter_compound, smiles, name)
    return cand, "rule_engine"


# ---------- 即时编码（缺的才算）----------

def embeddings_for(smiles_list: list[str]) -> tuple[dict, int]:
    emb, dim = ({}, 0)
    cache_p = kio.FEATURES_DIR / "mol_embeddings.parquet"
    if cache_p.exists():
        emb, dim = rxnfeat.load_embeddings()
    need = [s for s in smiles_list if s and s not in emb]
    if need:
        kio.log.info("即时编码 %d 个新结构 (ChemBERTa)...", len(need))
        from encode_molecules import FrozenEncoder
        cfg = kio.load_config()
        enc = FrozenEncoder(cfg["encoder"]["primary"], device=cfg["train"].get("device", "cuda"))
        vecs = enc.encode(need)
        for s, v in zip(need, vecs):
            emb[s] = np.asarray(v, dtype=np.float32)
        dim = dim or vecs.shape[1]
    return emb, dim


# ---------- 适用域 (applicability domain) ----------

def ad_score(sub_smiles: str, emb: dict, dim: int) -> tuple[float, str]:
    """substrate embedding 到 71 个训练底物的最大余弦相似度。低 = 域外，别信。"""
    z = emb.get(sub_smiles)
    if z is None:
        return float("nan"), "unknown"
    edges_p = kio.DATASETS_DIR / "kernel_edges.parquet"
    train_subs = pd.read_parquet(edges_p).query("y==1")["substrate_smiles"].dropna().unique()
    zt = rxnfeat._l2(z)
    sims = [float(np.dot(zt, rxnfeat._l2(emb[s]))) for s in train_subs if s in emb and s != sub_smiles]
    if not sims:
        return float("nan"), "unknown"
    mx = max(sims)
    flag = "in-domain" if mx >= 0.80 else "borderline" if mx >= 0.60 else "OUT-OF-DOMAIN(低置信)"
    return mx, flag


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--name", default=None, help="物质名称/缩写（如 2'-FL），自动解析成 SMILES")
    ap.add_argument("--smiles", default=None, help="直接给 SMILES（给了就不走名称解析）")
    ap.add_argument("--topn", type=int, default=15)
    ap.add_argument("--no-cache", action="store_true", help="即使缓存命中也强制走规则引擎")
    ap.add_argument("--out", default=None)
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()
    kio.setup_logging()
    if not args.name and not args.smiles:
        raise SystemExit("请给 --name <名称> 或 --smiles <SMILES>")

    if args.smiles:
        smi = kio.canonical_smiles(args.smiles)
        if smi is None:
            raise SystemExit(f"无法解析 SMILES: {args.smiles}")
        label = args.name or "query"
    else:
        from resolve import name_to_smiles
        kio.log.info("解析名称 '%s' → SMILES ...", args.name)
        smi, src = name_to_smiles(args.name)
        if smi is None:
            raise SystemExit(f"无法解析名称 '{args.name}'（PubChem/ChEBI 都没命中）。请改用 --smiles。")
        kio.log.info("  解析成功 [%s]: %s", src, smi)
        label = args.name
    args.name = "".join(c if c.isalnum() or c in "-_" else "_" for c in label)  # 文件名安全
    kio.log.info("== 预测底物 %s ==\n  canonical: %s", label, smi)

    cand, src = candidates_for(smi, args.name, use_cache=not args.no_cache)
    if cand is not None and not cand.empty:
        cand = cand.rename(columns={"rclss_production_score": "rclss_prior"})
        if "rclss_prior" not in cand.columns:
            cand["rclss_prior"] = np.nan
    if cand is None or cand.empty:
        kio.log.warning("规则引擎对该底物未生成任何候选（generation miss）。无法排序。")
        print(json.dumps({"substrate": smi, "name": args.name, "candidate_source": src,
                          "n_candidates": 0, "note": "no rule-generated candidates"}, ensure_ascii=False, indent=2))
        return
    kio.log.info("候选数: %d (source=%s)", len(cand), src)

    # 编码 + 打分
    emb, dim = embeddings_for([smi] + cand["product_smiles"].dropna().tolist())
    model, ckpt = load_final_model()
    import torch
    dev = "cuda" if torch.cuda.is_available() else "cpu"
    model.to(dev)
    X, ok = rxnfeat.build_features([smi] * len(cand), cand["product_smiles"].tolist(), emb, dim,
                                   ckpt.get("l2_normalize", True))
    cand = cand.reset_index(drop=True).copy()
    cand["ml_kernel_score"] = np.nan
    cand.loc[ok, "ml_kernel_score"] = score_model(model, X[ok], dev)
    cand["ml_confidence_tier"] = cand["ml_kernel_score"].map(
        lambda s: confidence_tier(s) if pd.notna(s) else "unscored")

    # 适用域
    ad_sim, ad_flag = ad_score(smi, emb, dim)

    ranked = cand.sort_values("ml_kernel_score", ascending=False).reset_index(drop=True)
    ranked["ml_rank"] = np.arange(1, len(ranked) + 1)

    # 产物加名字：本地 block1→name 全量映射；top-N 未命中再 PubChem 在线兜底
    from resolve import block1_to_name, ik_to_name_pubchem
    if "product_block1" in ranked.columns:
        ranked["product_name"] = ranked["product_block1"].map(block1_to_name)
    else:
        ranked["product_name"] = None
    for i in ranked.head(args.topn).index:
        nm = ranked.at[i, "product_name"]
        if nm is None or (isinstance(nm, float) and pd.isna(nm)) or not str(nm).strip():
            ik = ranked.at[i, "product_inchikey"] if "product_inchikey" in ranked.columns else None
            ranked.at[i, "product_name"] = ik_to_name_pubchem(ik) or "unnamed/novel"
    ranked["product_name"] = ranked["product_name"].fillna("")

    out_path = args.out or f"predictions/predict_{args.name}.parquet"
    kio.write_table(ranked, out_path, force=args.force)
    # 溯源链 JSON
    prov = {"input": {"name": args.name, "smiles": smi,
                      "inchikey": kio.smiles_to_inchikey(smi)},
            "applicability_domain": {"max_cosine_to_training_substrates": None if np.isnan(ad_sim) else round(ad_sim, 3),
                                     "flag": ad_flag},
            "candidate_source": src, "n_candidates": int(len(ranked)),
            "honest_boundary": "ml_kernel_score = 候选可信度/初筛分(主要由产物代谢物可信度驱动, decoy-bias)，非湿实验概率；AD 低=域外别信。",
            "top": []}
    for _, r in ranked.head(args.topn).iterrows():
        prov["top"].append({
            "rank": int(r["ml_rank"]), "product_name": r.get("product_name", ""),
            "product_smiles": r["product_smiles"],
            "product_inchikey": r.get("product_inchikey"),
            "ml_kernel_score": round(float(r["ml_kernel_score"]), 4) if pd.notna(r["ml_kernel_score"]) else None,
            "ml_confidence_tier": r["ml_confidence_tier"],
            "rclss_prior": None if pd.isna(r.get("rclss_prior")) else round(float(r["rclss_prior"]), 4),
            "rule_uid": r.get("rule_uid", "")})
    kio.safe_output_path(f"predictions/predict_{args.name}.json").write_text(
        json.dumps(prov, ensure_ascii=False, indent=2), encoding="utf-8")

    # 控制台
    kio.log.info("适用域 AD: max_cos=%.3f → %s", ad_sim if not np.isnan(ad_sim) else -1, ad_flag)
    show = ranked.head(args.topn)[["ml_rank", "product_name", "product_smiles",
                                   "ml_kernel_score", "ml_confidence_tier", "rclss_prior"]]
    print(f"\n=== {args.name} top-{args.topn} 预测产物 (source={src}, AD={ad_flag}) ===")
    print(show.to_string(index=False, max_colwidth=38))
    print(f"\n完整结果: outputs/{out_path}\n溯源JSON: outputs/predictions/predict_{args.name}.json")


if __name__ == "__main__":
    main()
