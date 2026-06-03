"""
encode_molecules.py — Phase 2

用**冻结**的化学语言模型（ChemBERTa-2 / MolFormer）把 SMILES 编码成固定向量，缓存到
outputs/features/mol_embeddings.parquet（key = canonical SMILES）。

设计：
  - 打分路径 RDKit-free：SMILES → HF tokenizer → 冻结 LM → mean-pool（mask-aware）→ emb。
  - 编码器冻结（no_grad）：低数据（~90 正样本）下避免过拟合，且可一次性缓存。
  - canonical SMILES 在 build_kernel_dataset 已做（RDKit 离线）；这里只读已 canonical 的 SMILES。

用法：
  python scripts/ssrf/ml_ranking_kernel/src/encode_molecules.py [--smoke] [--force] [--encoder primary|fallback]
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import kio


def collect_smiles(cfg: dict, smoke: bool) -> list[str]:
    smis: set[str] = set()
    edges_p = kio.DATASETS_DIR / "kernel_edges.parquet"
    if edges_p.exists():
        e = pd.read_parquet(edges_p, columns=["substrate_smiles", "product_smiles"])
        smis |= set(e["substrate_smiles"].dropna())
        smis |= set(e["product_smiles"].dropna())
    cand_p = kio.DATASETS_DIR / "substrate_candidates.parquet"
    if cand_p.exists():
        c = pd.read_parquet(cand_p, columns=["substrate_smiles", "product_smiles"])
        smis |= set(c["substrate_smiles"].dropna())
        smis |= set(c["product_smiles"].dropna())
    smis = sorted(s for s in smis if isinstance(s, str) and s)
    if smoke:
        smis = smis[:200]
    return smis


class FrozenEncoder:
    def __init__(self, enc_cfg: dict, device: str = "cuda"):
        import torch
        from transformers import AutoTokenizer, AutoModel
        self.torch = torch
        self.name = enc_cfg["name"]
        self.hf = enc_cfg["hf_model"]
        self.max_length = int(enc_cfg.get("max_length", 256))
        self.batch_size = int(enc_cfg.get("batch_size", 256))
        self.pooling = enc_cfg.get("pooling", "mean")
        self.device = device if torch.cuda.is_available() else "cpu"
        kio.log.info("Loading frozen encoder %s (%s) on %s", self.name, self.hf, self.device)
        self.tok = AutoTokenizer.from_pretrained(self.hf, trust_remote_code=True)
        self.model = AutoModel.from_pretrained(self.hf, trust_remote_code=True,
                                               deterministic_eval=True) \
            if "MoLFormer" in self.hf else AutoModel.from_pretrained(self.hf, trust_remote_code=True)
        self.model.to(self.device).eval()
        for p in self.model.parameters():
            p.requires_grad_(False)
        self.dim = int(self.model.config.hidden_size)

    def encode(self, smiles: list[str]) -> np.ndarray:
        torch = self.torch
        vecs = []
        for i in range(0, len(smiles), self.batch_size):
            batch = smiles[i:i + self.batch_size]
            enc = self.tok(batch, padding=True, truncation=True,
                           max_length=self.max_length, return_tensors="pt").to(self.device)
            with torch.no_grad():
                out = self.model(**enc)
            hidden = out.last_hidden_state                  # (B, T, H)
            mask = enc["attention_mask"].unsqueeze(-1).float()  # (B, T, 1)
            summed = (hidden * mask).sum(1)
            counts = mask.sum(1).clamp(min=1e-9)
            mean = (summed / counts).float().cpu().numpy()
            vecs.append(mean)
            if (i // self.batch_size) % 10 == 0:
                kio.log.info("  encoded %d/%d", min(i + self.batch_size, len(smiles)), len(smiles))
        return np.vstack(vecs) if vecs else np.zeros((0, self.dim), dtype=np.float32)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--encoder", choices=["primary", "fallback"], default="primary")
    ap.add_argument("--smoke", action="store_true")
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    kio.setup_logging()
    cfg = kio.load_config()
    enc_cfg = cfg["encoder"][args.encoder]

    smis = collect_smiles(cfg, args.smoke)
    kio.log.info("Encoding %d unique SMILES with %s", len(smis), enc_cfg["name"])
    if not smis:
        raise SystemExit("No SMILES found. Run build_kernel_dataset.py first.")

    enc = FrozenEncoder(enc_cfg, device=cfg["train"].get("device", "cuda"))
    embs = enc.encode(smis)
    kio.log.info("Embedding matrix: %s", embs.shape)

    out = pd.DataFrame({"smiles": smis})
    out["encoder"] = enc_cfg["name"]
    out["embedding"] = list(embs.astype(np.float32))
    suffix = "" if args.encoder == "primary" else f".{enc_cfg['name']}"
    kio.write_table(out, kio.FEATURES_DIR / f"mol_embeddings{suffix}.parquet", force=args.force)
    # 维度元数据
    import json
    meta = {"encoder": enc_cfg["name"], "hf_model": enc_cfg["hf_model"], "dim": int(embs.shape[1]),
            "n_smiles": len(smis), "pooling": enc_cfg.get("pooling")}
    (kio.FEATURES_DIR / "encoder_meta.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")
    kio.log.info("encoder dim=%d", embs.shape[1])


if __name__ == "__main__":
    main()
