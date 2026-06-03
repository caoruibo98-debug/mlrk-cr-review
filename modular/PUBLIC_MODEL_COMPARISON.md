# 与公开模型对比 — 定位 + 协议(诚实,不吹)

前提:公开模型也不完美,我们也不可能完美。要证明"我们的工作有意义",得在**性能**或**技术**上有可量化的差异。

## 1. 公开模型是什么(对比对象)
| 工具 | 干什么 | 输出 | 排序方式 |
|---|---|---|---|
| BioTransformer 3.0 | 代谢预测(EC/CYP450/gut/phase II) | 生成候选代谢物 | 反应似然/规则置信(启发式) |
| MicrobeRX | 肠道菌反应预测(**我们的规则源**) | 生成候选 | 规则/酶证据 |
| RetroRules | 反应模板(retro) | 模板匹配 | 模板 score |
| GLORYx / SyGMa / enviPath | 通用/CYP/环境代谢 | 生成候选 | 优先级规则/概率 |

**共性**:它们主要是**生成器**,排序靠**启发式**(规则置信/EC/概率),且**普遍不做泄漏审计、不分食品组分、不报 artifact**。

## 2. 我们 vs 它们(可量化)
| 维度 | 公开 | 我们 | 怎么证明 |
|---|---|---|---|
| 候选生成召回 | 规则 25–70% | 同源规则 25–68% | **打平**(诚实:用了相似规则,不声称多生成) |
| **排序质量** | 启发式(≈我们基线 ec_only/rclss/tanimoto) | **LTR_chem r@5 0.81–1.0** | 同一候选集上 LTR_chem > 这些启发式基线(已测) |
| 泄漏审计 | 罕见 | scaffold+产物不相交+4 confound 可复现 | CODE_REVIEW.md |
| 食品组分模块化 | 无 | A/B/C/D 各自模型 | module_router + per-module |
| 多证据可解释 | 部分(EC) | 结构+反应+菌+酶+文献 overlay + 适用域 | predict overlay |

**意义结论**:
- **性能**:不赢在"多生成",赢在**"同候选集上的学习排序 > 公开工具用的启发式排序"**(top-5 召回)。
- **技术**:赢在**泄漏审计严格 + 食品模块化 + 多证据 + 能自抓 artifact** 的方法学框架。

## 3. 真·头对头协议(待执行,需装 BioTransformer)
目标:在**同一批食品底物**上比 (a) 生成召回、(b) 排序 top-k 命中。
1. 取我们的 gold+silver 测试底物(structure 已 canonical)。
2. 跑 BioTransformer(`-b allHuman`/`superbio`,gut microbial)→ 得每底物候选代谢物。
3. 生成召回 = BioTransformer 是否产出真产物;对比我们规则的命中率。
4. 排序:BioTransformer 自带优先级 vs 我们 LTR_chem(在各自候选集上 recall@k)。
5. **公平点**:都在"各自生成命中"的子集上比排序;生成不同则分开报。
- 环境:BioTransformer 是 Java(`BioTransformer3.0.jar`),需 JRE + 下载。建议单列任务。
- 备选(无需装):MicrobeRX 是我们的规则源,等价于"我们生成 = MicrobeRX 生成";所以 vs MicrobeRX 的排序对比 = LTR_chem vs rclss/ec_only(已有)。

## 4. 现在已能下的诚实结论(无需再跑)
在同一规则候选集上:**LTR_chem 的 top-5 召回稳超公开工具所代表的启发式排序(EC/相似度/规则置信)**;生成召回与公开规则工具同档(诚实)。差异化主要是**排序 + 方法学严谨**,不是"生成更全"。
