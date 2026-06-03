# L-RCLSS — Learned, Rule-Agnostic Reaction Ranking Kernel

与反应规则无关、可迁移、打分路径独立于 RDKit 的**神经反应排序内核**。
把 SSRF/RCLSS 从"人工网格校准的规则相似性"升级为"学习到的可迁移排序内核"。

完整设计见 `C:\Users\crb66\.claude\plans\snug-launching-backus.md`。

## 一句话定位
- **阶段1（稳 / RDKit+SMARTS）**：production SSRF/RCLSS（`scripts/ssrf/`）仍是候选生成器 + SMARTS gate + 溯源来源。
- **阶段2（拓宽 / 学习 / RDKit-free 打分）**：冻结化学语言模型（ChemBERTa-2 / MolFormer）编码 substrate/product →
  反应表示 `[z_s, z_p, z_s−z_p, z_s⊙z_p]` → 轻量排序头（pairwise + nnPU）→ `ml_kernel_score`。
- 内核**不消费** rule_id/EC/RCLSS/category 作特征 → 可迁移到任意来源的候选。

## 成功判据（标准式）
同一 split 下 `ML > 规则先验RCLSS`；且 `ML(scaffold) ≈ ML(random)`（小泛化 gap = 学到可迁移化学）。

## 安全契约（沿用 trainable_ranker_experiment 模式）
- **不修改** production SSRF/RCLSS/V3/V4 任何代码或产物；只读取它们的输出 + label 文件。
- 所有写出落在 `scripts/ssrf/ml_ranking_kernel/outputs/` 内（`kio.safe_output_path` 强制）；
  已存在文件需 `--force`。

## 流水线
```
src/build_kernel_dataset.py   # Phase 1  正/负/未标注边表 → outputs/datasets/kernel_edges.parquet
src/encode_molecules.py       # Phase 2  冻结化学LM embedding 缓存 → outputs/features/mol_embeddings.parquet
src/make_splits.py            # Phase 3  scaffold + random split + scaffold k-fold CV
src/train_kernel.py           # Phase 4  pairwise + nnPU 排序头（规则无关特征闸）
src/fuse_hybrid.py            # Phase 5  ML-only / gate×ML / late-fusion 变体
src/evaluate_kernel.py        # Phase 6  {RCLSS,ML,hybrid}×{scaffold,random} 矩阵 + success_gate.json
src/predict_kernel.py         # Phase 7  打分 + 溯源 JSON + 非破坏接入 ssrf_filter
src/build_network_negatives.py# Phase 1b AGREDA 隐式负样本（阶段化，不阻塞首版）
```

## 运行
```bash
cd "D:\CRB\Food models"
python scripts/ssrf/ml_ranking_kernel/src/build_kernel_dataset.py
python scripts/ssrf/ml_ranking_kernel/src/encode_molecules.py
python scripts/ssrf/ml_ranking_kernel/src/make_splits.py
python scripts/ssrf/ml_ranking_kernel/src/train_kernel.py
python scripts/ssrf/ml_ranking_kernel/src/evaluate_kernel.py
```
所有脚本支持 `--smoke`（小样本快跑）与 `--force`（覆盖本包 outputs）。
