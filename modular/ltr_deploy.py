"""
ltr_deploy.py — 落地 LTR_chem 部署模型

在全部干净候选(clean_candidates)上、每模块训练 LTR_chem(纯 ChemBERTa 结构特征)终模型并保存。
这是诚实可部署的排序器：超裸 Tanimoto、不吃 EC/数据库 confound(无 EC 特征)。

输出 modular/outputs/ltr/models_clean/ltr_chem_{A,B,C,D}.json + meta
用法：python scripts/ssrf/ml_ranking_kernel/modular/ltr_deploy.py [--candidates clean_candidates.parquet] [--force]
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
from ltr_build import LTR_OUT  # noqa: E402
from ltr_train_eval import encode_all  # noqa: E402
from ltr_v3_train_eval import build_feats, train  # noqa: E402


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--candidates", default="clean_candidates.parquet")
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()
    kio.setup_logging()
    C = pd.read_parquet(LTR_OUT / args.candidates)
    emb, dim = encode_all(set(C.substrate_smiles) | set(C.product_smiles))
    models_dir = kio.safe_output_path("modular/ltr/models_clean/_m").parent
    meta_all = {}
    for m in ["A", "B", "C", "D"]:
        Cm = C[C.module == m]
        if (Cm.y == 1).sum() < 10:
            kio.log.info("module %s: too few positives, skip", m); continue
        # LTR_chem：纯结构（ue=False=不用EC, ut=False=不用Tanimoto特征）
        X, ok = build_feats(Cm, emb, dim, {}, False, False)  # fpc 空：chem-only 不需要指纹
        pk = Cm[ok]; qid = pd.factorize(pk.sb)[0]
        model = train(X[ok], pk.y.to_numpy(np.float32), qid)
        model.save_model(str(models_dir / f"ltr_chem_{m}.json"))
        meta_all[m] = {"module": m, "feat": "chem_only_reaction_diff", "dim": dim,
                       "feat_dim": 4 * dim, "n_train_rows": int(len(pk)),
                       "n_pos": int((pk.y == 1).sum()), "n_substrates": int(pk.sb.nunique())}
        kio.log.info("module %s: trained LTR_chem on %d rows (%d pos, %d subs) -> saved",
                     m, len(pk), int((pk.y == 1).sum()), pk.sb.nunique())
    (models_dir / "deploy_meta.json").write_text(json.dumps(meta_all, indent=2, ensure_ascii=False), encoding="utf-8")
    kio.log.info("== deploy models saved: %s ==", list(meta_all.keys()))


if __name__ == "__main__":
    main()
