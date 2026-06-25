# CAAM-scMAE BDD implementation plan

本目录整理 CAAM-scMAE 最新 BDD 与补充说明，作为后续实现、实验和论文决策的阶段化执行参考。

## 当前研究决策状态

CAAM-v1 已经完成工程 smoke 与 formal 接入 smoke，但早期实验暴露出两类问题：

```text
1. formal benchmark 曾强制 n_top_genes=0，导致 full-gene 输入、参数量和 mask budget 结论偏悲观；这不应作为主 benchmark 协议。
2. matched donor + strict effective budget 把 zero-to-zero corruption 当作硬失败；这比 scMAE/scNAME 的实际 masked modeling 协议更严格。
```

因此，本 BDD 现在进入 **Correction Pipeline**：先修正输入协议和 corruption 机制，再判断 AdvMask/Axial 是否有研究价值。不要继续把 `Axial + AdvMask synergy` 当作默认论文主线。

## 适用范围

目标包：

```text
methods/DeepLearning/CAAM_scMAE/
```

对外正式方法名仍为：

```text
caam_scmae
```

正式 benchmark 仍只允许注册最终方法名；内部消融和 correction variants 不得进入正式主方法列表。

## 当前方法身份

历史 CAAM-v1 的内部模型为：

```text
Model 0: Controlled-scMAE
Model A: Axial-scMAE
Model B: AdvMask-scMAE
Model C: CAAM-scMAE
```

但 Correction Pipeline 中，研究优先级临时调整为：

```text
1. 先比较 corruption 机制：scMAE-style shuffle / matched donor / nonzero-aware donor
2. 再判断 AdvMask 是否相对 random mask 有稳定增益
3. 最后才允许 Axial 重新进入 2x2 factorial 消融
```

## 文档优先级

执行时按以下优先级理解需求：

```text
correction_bdd_index.md
> literature_and_problem_reframing.md
> addendum_formal_benchmark_interface.md
> benchmark_contract.md
> phases/12_protocol_correction.md 到 phases/16_publication_decision.md
> 补充说明
> CAAM-scMAE 完整 BDD v1.0
> natural_language_model_overview.md
> 早期 BDD
> 早期模型设想说明
```

若早期 BDD 与 Correction Pipeline 冲突，以 Correction Pipeline 为准。

## 目录内容

```text
correction_bdd_index.md
  当前总控修订 BDD：冻结 CAAM-v1 主张，启动输入协议、corruption 和论文路线修正。

literature_and_problem_reframing.md
  汇总 scMAE、scNAME、scCluBench、sciLaMA 等相关工作对本项目的启发和边界。

addendum_formal_benchmark_interface.md
  正式 benchmark 接口补丁；仍约束 run.py 的 CLI/artifact 设计。

natural_language_model_overview.md
  CAAM-scMAE 的任务、目标、模型框架、理论支撑与重点风险的自然语言说明；用于帮助 Codex 理解研究思想，不替代工程约束。

execution_order.md
  阶段总顺序与验收边界，现已增加 Phase 12-16 correction pipeline。

benchmark_contract.md
  正式 benchmark、内部 ablation、GPU、unknown-K、artifact、HVG/full-gene 协议契约。

test_matrix.md
  BDD 要求的测试矩阵与阶段归属。

risk_and_stop_criteria.md
  高风险点、人工核查项、fail-fast 与止损条件；现在将 budget deficit 默认降级为诊断指标。

experiment_split_decision.md
  smoke/development/validation/sealed-test 数据集边界。

phases/
  每个阶段的具体操作步骤、产物和验收条件。
```

## Correction Pipeline 总原则

1. 主 benchmark 默认使用 label-free HVG feature space，而不是 full-gene input。
2. full-gene input 只作为 scalability stress test，不作为主表协议。
3. zero-to-zero corruption 不再默认 fail-fast；必须记录为 `zero_to_zero_rate` 和 `effective_corruption_rate`。
4. strict effective budget 只能作为显式开关，不得默认阻断训练。
5. matched donor 不得被默认视为正确；必须与 scMAE-style gene-wise shuffle 和 nonzero-aware donor 比较。
6. AdvMask 必须先在 MLP encoder 下证明相对 random mask 有稳定增益，才允许进入 full 模型。
7. Axial 必须在参数量、runtime 和 parameter-matched MLP 对照下证明价值，才允许作为论文主模块。
8. 不让真实标签进入训练、donor、context、gene module、mask selector、early stopping 或 resolution selection。
9. 不把 known-K 结果称作 unknown-K 或完全无监督。
10. 不把 correction variants 注册到 formal benchmark。
11. 遇到影响论文主张的机制变化，先写入 BDD 和 report，再执行代码。
12. sealed test datasets 不得用于调参、选 corruption、选 mask ratio 或决定是否保留 Axial/AdvMask。
