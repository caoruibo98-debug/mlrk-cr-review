# CODE REVIEW 指南 — 模块化食品组分代谢物排序器(给你审我的逻辑)

目的:让你能逐条核对"模型是不是在学结构,而不是作弊"。每条都给**文件:函数**,你可以直接打开看。
全部代码在 `scripts/ssrf/ml_ranking_kernel/modular/`(共享工具在 `../src/`)。production SSRF 未动。

---

## 0. 任务定义(先对齐"我们在排什么")
给一个底物 S,在一组**候选产物**里把"真产物"排到前面。指标 = 真产物的 recall@k / MRR(per-substrate)。
- "候选产物"从哪来 = **作弊高发区**(见 §3)。
- 排序器 = XGBoost LambdaMART(`ltr_v3_train_eval.py:train`)。Tanimoto/ec_only/random 是**基线**(不学习)。

## 1. 数据流(按执行顺序)
```
build_food_test.py       # 单步文献正样本(gold) + 路由模块         → food_clean.parquet
ltr_build.py             # Codex 池 → 唯一反应 + tier(gold/silver/weak) → reactions.parquet
ltr_clean_candidates.py  # 对每底物规则 RunReactants 生成候选(部署忠实)  → clean_candidates.parquet  ★核心
ltr_clean_eval.py        # scaffold+产物不相交 多seed 评测 {6 方法}        → clean_metrics.csv        ★核心
ltr_deploy.py            # 全量训练 LTR_chem 终模型                       → models_clean/ltr_chem_*.json
predict_substrate_clean  # 部署推理:query→RunReactants→LTR_chem 排序+证据overlay (待接)
```
共享:`src/module_router.py`(结构路由 A/B/C/D)、`src/encode_molecules.py`(冻结 ChemBERTa)、`ltr_v3_train_eval.py`(build_feats/train/within_sub_eval)。

## 2. 特征 — 审"有没有混进会泄漏的东西"
`ltr_v3_train_eval.py:build_feats(df, emb, dim, fpc, use_ec, use_tan)`:
- 结构(永远在):`rxn_feats` = `[z_s, z_p, z_s−z_p, z_s⊙z_p]`,z 来自冻结 ChemBERTa(`encode_molecules.FrozenEncoder`)。**纯分子结构,无规则身份。**
- `use_ec`:候选生成规则的 EC class onehot(7 维)。**部署可用**(候选带其生成规则),但见 §3 第4个 confound。
- `use_tan`:Tanimoto(S,P) 标量。部署可用。
- **LTR_chem = use_ec=False, use_tan=False = 只结构**。这是我推荐部署的版本(不碰 EC/相似度元数据)。
- 审点:确认没有 rule_id / 数据库 id / label 衍生列进特征。`feature_block_names` 全是 emb 维。

## 3. 反作弊控制(★ 重点审这里)
排序好不好,关键看**候选集是怎么造的**——我犯过/堵过的 4 个 confound:

| # | 作弊机制 | 在哪检测到 | 怎么堵的 | 文件:位置 |
|---|---|---|---|---|
| 1 | decoy 是 RunReactants 怪结构→模型学"真代谢物 vs 怪结构" | 特征重要性 + 消融 | 部署忠实:候选就是规则真实产出(含怪的),模型必须在其中排;不再人为造"干净 vs 怪"对照 | `ltr_clean_candidates.py:main`(无 QED 偏置、无 vocab 注入) |
| 2 | scaffold 切分只留出底物,**产物跨训练/测试复现**→记忆产物 | 测试真产物 49–63% 也在训练 | **产物不相交切分**:训练里剔除所有测试真产物的正样本 | `ltr_clean_eval.py`: `tr = tr[~((tr.y==1)&(tr.pb.isin(test_true)))]` |
| 3 | 给 vocab decoy 打 EC=0、真产物有真 EC→学"有没有 EC 标签" | `ec_only` 基线 r@1=0.81 | **候选全用规则生成,EC 用生成规则的真实值(真/假一视同仁),不注入** | `ltr_clean_candidates.py`(cand 字典存规则 EC,真假同源) |
| 4 | 真产物的生成规则带 EC、decoy 规则常 EC=0(≈数据库成员偏置) | `ec_only` 仍偏高(r@1 .48–.80) | 部署用 **LTR_chem(无 EC 特征)**;并用"只 EC 规则"中和run(`--require-ec`)确认 LTR_chem 仍超 Tanimoto | `ltr_clean_eval.py:ec_score`(常驻探针) + `ltr_clean_candidates.py:--require-ec` |

另外两个永久守卫:
- **底物 scaffold 切分**(`ltr_clean_eval.py` GroupKFold on `scaf`):测试底物的骨架不在训练里 → 不能记忆底物。
- **generation miss 诚实剔除**:规则没生成出真产物的底物直接丢(`ltr_clean_candidates.py`: `if not gen_truth: continue`),并单独报 generation recall(=硬天花板,A25/B68/C45/D59%)。不偷偷注入真产物来刷分。

## 4. 作弊检测器(你可以一直拿它审我)
`ltr_clean_eval.py` 里 **`ec_only`**(score = `ec_class != 0`)和 **`random`** 是常驻基线。
- 任何"靠元数据虚高"的 artifact,`ec_only`/`random` 会跟着高 → 一眼看穿。
- 判据写死:**LTR 要同时 > tanimoto 且 > ec_only,r@5 上才算真信号**。
- `ltr_clean_candidates.py` 末尾打印 `[anti-cheat check] EC==0 rate pos vs decoy`,两者应接近(泄漏时会拉开)。

## 5. 结果怎么读(clean_metrics.csv)
- **r@1**:ec_only(confound)偏高,别用它下结论。
- **r@5 / MRR**:LTR_chem 全模块最高(r@5 0.81–1.0,超 tanimoto 和 ec_only)→ 这是真正的可部署信号。
- 特征重要性:structure≈0.99,EC≈0,Tanimoto≈0 → LTR 确实在用结构,没靠 EC。

## 6. 你审计时重点盯的地方(我自己列的薄弱点)
1. `ltr_clean_candidates.py` 的 EC 赋值:真产物和 decoy 的 ec_class 是否**真的同源**(都来自 `cand` 字典的生成规则)。这是堵 confound#3/#4 的关键。
2. `ltr_clean_eval.py` 产物不相交那一行是否对**所有** LTR 变体生效(是,tr 在变体循环外过滤)。
3. scaffold 列 `scaf` 是否真按 Bemis-Murcko(`kio.murcko_scaffold`),空值是否退化成"每底物自成一组"(是,`fillna("NS_"+sb)`)。
4. `within_sub_eval` 的 tie 处理:argsort 稳定排序,真产物**不会**因行序占便宜?→ ⚠️ 注意:clean_candidates 里真产物在前,若分数并列会被排前。LTR 分数连续(无并列)不受影响;但 `ec_only`/`random` 有并列 → 它们的 r@1 可能被行序高估。**这是 ec_only 偏高的部分原因之一,审计时要知道。**(修法:eval 里候选先 shuffle。)

## 7. 一句话给审稿人/你
排序器只吃分子结构(冻结 ChemBERTa),候选来自规则真实产出,切分同时切底物骨架和产物,且有 `ec_only`/`random` 常驻探针 + 4 个已修 confound 的可复现记录。**LTR_chem 在 top-5 召回上稳超相似度且不碰元数据 → 这是诚实可部署的增量。**
