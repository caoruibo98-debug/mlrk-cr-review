# Model Card — L-RCLSS v2 Modular Metabolite-Retrieval Ranker

生成：2026-05-30。位置 `scripts/ssrf/ml_ranking_kernel/modular/`（v1 保留作对照，production SSRF 未动）。

## 1. 这个模型做什么
输入一个食品组分物质（名称/SMILES）→ 结构路由到模块 A 小分子 / B 碳水 / C 蛋白 / D 脂肪 →
用该模块的 **LambdaMART(XGBoost) 排序器**在**该模块已知产物词表**里检索排序，输出 top-N 候选代谢物，
附：产物名 + 排序分 + 适用域(AD) + 证据 overlay(已知菌/酶/文献)。

**任务 = 检索（在已知代谢物词表里排），不是生成全新结构。** 排序分非湿实验概率。

## 2. 数据
- 来源：`FoodGut/.../foodgut_modular_positive_candidate_pool_v2.csv`（Codex 策展）。
- **10,550 条唯一单步反应**；训练对 168,800（正 10,550 / 负 158,250，负=同模块产物词表采样的"错配真代谢物"）。
- **干净文献正样本 80 条**（manual_literature_curated）= 唯一诚实测试集；A41/B13/C10/D16。
- 弱标签（MicrobeRX/VMH/DB 派生）= 预训练，每模块封顶 2500 底物。
- **泄漏控制**：干净反应的弱重复行被去掉；scaffold-CV 按底物骨架分组；干净集从未进调参。

## 3. 特征
- 结构/反应：冻结 ChemBERTa(384) 反应表示 `[z_s,z_p,z_s−z_p,z_s⊙z_p]`（1536）。
- 证据（菌/酶/文献）：**实测为"答案泄漏"**——加进排序使 r@10 腰斩（A 0.75→0.37 等），因检索时新颖候选无证据。
  → **证据不进排序，只作置信/解释 overlay。**

## 4. 评估与指标（干净集 scaffold 5-fold，per-fold mean±std）
| 模块 | LTR r@10 | **Tanimoto r@10** | random r@10 | LTR MRR | Tan MRR |
|---|---|---|---|---|---|
| A | 0.75±0.21 | **0.86±0.13** | 0.03 | 0.46 | 0.43 |
| B | 0.83±0.17 | 0.75±0.26 | 0.04 | 0.45 | 0.42 |
| C | 0.83±0.21 | **0.98±0.05** | 0.29 | 0.38 | 0.36 |
| D | 0.86±0.08 | 0.81±0.14 | 0.29 | 0.47 | 0.30 |

**关键结论（诚实）**：LTR 与 **Tanimoto 结构相似度基线打平**——A/C 上 Tanimoto 更好，B/D 上 LTR 略好但 CI 重叠。
两者都远超 random（任务可学），但**学习模型尚未证明比一行 Tanimoto 有增量**。原因：代谢物多是底物小改动，
相似度已近最优；当前负样本是跨底物词表采样（偏易），没逼模型学"具体变换"。

时间留出（新文献召回）严重欠样本（总 n≈5），不能下结论。

## 5. 适用边界（不能 claim）
- 不是湿实验准确率；不测"全新结构生成"（只在已知词表检索）。
- **未证明优于 Tanimoto**；n=18–47，CI 宽，单次 run。
- D 脂肪 / C 蛋白 干净测试仅 10–18 条。
- 证据 overlay 只覆盖"已知反应"。

## 6. 部署
- 终模型：`outputs/modular/ltr/models/ltr_{A,B,C,D}.json` + meta。
- 推理：`predict_modular.py --name <名称>|--smiles <SMILES>` → top-N + 证据 + AD + 名字 + provenance JSON。

## 7. 要真正"超过基线"的下一步
1. **同底物 hard 负样本**（同一底物经其它规则 RunReactants 的错产物）替代跨底物词表负 → 逼模型学变换而非相似度。
2. 反应类型/规则条件化特征。
3. Tanimoto 作显式特征/集成（至少不输基线）。
4. 扩近年干净反应，把时间留出做实。
5. 嵌套 CV 调参已做（weak 上选参），下一步多 seed 方差。
