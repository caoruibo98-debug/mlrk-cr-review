"""
predict_kernel.py — Phase 7

用训练好的 L-RCLSS 内核给候选 (substrate→product) 边打分并排序，输出**完整溯源链 JSON**
（design_decisions.md 格式）+ ml_kernel_score + confidence tier。

保持可解释（调和 CLAUDE.md 原则3）：规则仍是机制来源；ML 只重排。
非破坏：ml 分作为附加列写到本包 outputs/predictions/，不改 production ssrf_filter 产物。

用法：
  # 给某底物（用 substrate_candidates 的候选）排序：
  python scripts/ssrf/ml_ranking_kernel/src/predict_kernel.py --substrate-block1 <BLOCK1> [--topn 20]
  # 或给一张候选边表（含 substrate_smiles/product_smiles[/rule_uid/rclss_prior]）打分：
  python scripts/ssrf/ml_ranking_kernel/src/predict_kernel.py --edges path.parquet --out scored.parquet
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
from train_kernel import make_mlp, score_model


def load_final_model():
    import torch
    ckpt = torch.load(kio.MODELS_DIR / "kernel_final.pt", map_location="cpu", weights_only=False)
    model = make_mlp(ckpt["in_dim"], ckpt["hidden"], ckpt["dropout"])
    model.load_state_dict(ckpt["state_dict"])
    model.eval()
    return model, ckpt


def confidence_tier(score: float) -> str:
    return "high" if score >= 0.66 else "medium" if score >= 0.33 else "low"


def score_edges(edges: pd.DataFrame) -> pd.DataFrame:
    import torch
    model, ckpt = load_final_model()
    emb, dim = rxnfeat.load_embeddings()
    assert dim == ckpt["dim"], "embedding dim mismatch vs trained model"
    X, ok = rxnfeat.build_features(edges["substrate_smiles"], edges["product_smiles"], emb, dim,
                                   ckpt.get("l2_normalize", True))
    out = edges.copy()
    out["ml_kernel_score"] = np.nan
    dev = "cuda" if torch.cuda.is_available() else "cpu"
    model.to(dev)
    out.loc[ok, "ml_kernel_score"] = score_model(model, X[ok], dev)
    out["ml_confidence_tier"] = out["ml_kernel_score"].map(
        lambda s: confidence_tier(s) if pd.notna(s) else "unscored")
    return out


def provenance_row(r: pd.Series) -> dict:
    """design_decisions.md 溯源链格式（ML 重排版）。"""
    return {
        "input": {"smiles": r.get("substrate_smiles"), "inchikey": r.get("substrate_inchikey", "")},
        "reaction_applied": {"rule_uid": r.get("rule_uid", ""),
                             "rclss_prior_score": _f(r.get("rclss_prior"))},
        "product": {"smiles": r.get("product_smiles"), "inchikey": r.get("product_inchikey", "")},
        "ml_kernel": {"score": _f(r.get("ml_kernel_score")),
                      "confidence_tier": r.get("ml_confidence_tier"),
                      "model": "L-RCLSS frozen-ChemBERTa reaction-diff MLP (rule-agnostic)"},
        "rank": int(r["ml_rank"]) if pd.notna(r.get("ml_rank")) else None,
    }


def _f(v):
    try:
        return None if v is None or (isinstance(v, float) and np.isnan(v)) else float(v)
    except Exception:
        return None


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--substrate-block1", default=None)
    ap.add_argument("--edges", default=None, help="候选边 parquet/csv（含 substrate_smiles, product_smiles）")
    ap.add_argument("--topn", type=int, default=20)
    ap.add_argument("--out", default=None)
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()
    kio.setup_logging()

    if args.edges:
        edges = kio.read_table(Path(args.edges))
    else:
        cand = pd.read_parquet(kio.DATASETS_DIR / "substrate_candidates.parquet")
        cand = cand.rename(columns={"rclss_production_score": "rclss_prior"})
        if args.substrate_block1:
            cand = cand[cand["substrate_block1"].astype(str).str.upper() == args.substrate_block1.upper()]
        edges = cand
    if edges.empty:
        raise SystemExit("No candidate edges to score.")

    scored = score_edges(edges)
    scored = scored.sort_values(["substrate_block1", "ml_kernel_score"], ascending=[True, False])
    scored["ml_rank"] = scored.groupby("substrate_block1")["ml_kernel_score"].rank(ascending=False, method="min")

    out_path = args.out or "predictions/kernel_scored_candidates.parquet"
    kio.write_table(scored, out_path, force=args.force)

    # 溯源链 JSON（top-N per substrate）
    prov = []
    for sub, g in scored.groupby("substrate_block1"):
        for _, r in g.nsmallest(args.topn, "ml_rank").iterrows():
            prov.append(provenance_row(r))
    kio.safe_output_path("predictions/kernel_provenance.json").write_text(
        json.dumps(prov, indent=2, ensure_ascii=False), encoding="utf-8")
    kio.log.info("scored %d edges; wrote provenance for %d ranked candidates", len(scored), len(prov))
    if args.substrate_block1:
        top = scored.head(min(args.topn, len(scored)))[
            ["product_smiles", "ml_kernel_score", "ml_confidence_tier", "rclss_prior"]]
        kio.log.info("top candidates for %s:\n%s", args.substrate_block1, top.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
