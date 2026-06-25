# Phase 16: Publication decision

## 目标

根据 corrected protocol 的真实结果，决定 CAAM 是否继续作为正向方法论文推进，还是转为协议分析论文或停止。

本阶段不写新模型，不调参，只做研究决策。

## previous_phase_check

执行前必须检查：

```text
1. Phase 12 protocol correction 是否完成。
2. Phase 13 corruption triad 是否完成。
3. Phase 14 AdvMask triage 是否完成。
4. 如果 Phase 15 被执行，是否满足 Axial re-entry gate。
5. 是否没有使用 sealed test 做机制选择。
6. 是否所有 reports 都记录了 seeds、datasets、n_top_genes、corruption_type。
```

如果 Phase 14 已经失败，则 Phase 15 可以没有执行，但必须明确记录原因。

## 决策输入

必须读取：

```text
PHASE12_PROTOCOL_CORRECTION_REPORT.md
PHASE13_CORRUPTION_TRIAD_REPORT.md
PHASE14_ADVMASK_TRIAGE_REPORT.md
PHASE15_AXIAL_REENTRY_REPORT.md if exists
experiment_split_decision.md
```

## 决策路线 A：正向方法论文

只有满足以下条件，才进入方法论文路线：

```text
1. 最佳 corrected method 在 development datasets 上稳定优于 scMAE-style baseline。
2. validation datasets 上保持趋势。
3. 提升超过 seed 波动。
4. biological interpretation 不下降。
5. runtime 和参数量可接受。
6. 与 scMAE/scNAME/DESC/scDeepCluster/scDCC/Leiden/Louvain 等基线比较公平。
```

可写贡献：

```text
sparsity-aware or effectiveness-aware masked corruption
robust mask selection under corrected protocol
scalable masked autoencoding for scRNA clustering
```

## 决策路线 B：协议分析论文

如果方法提升不稳定，但以下现象稳定：

```text
1. full-gene vs HVG 显著影响 masked AE。
2. zero_to_zero_rate/effective_corruption_rate 与表示质量有关。
3. matched donor、scMAE-style shuffle、nonzero-aware donor 在不同 sparsity 下表现不同。
4. AdvMask 或 Axial 在部分条件下失败可解释。
```

则转为 protocol analysis / diagnostic paper。

可写贡献：

```text
masked corruption protocol matters in scRNA-seq clustering
HVG feature space changes masked reconstruction behavior
effective corruption diagnostics predict training usefulness
```

## 决策路线 C：停止 CAAM

若满足任一条件，停止：

```text
1. corrected protocol 下最佳模型仍不优于 scMAE-style control。
2. AdvMask 无稳定增益。
3. Axial 不优于 parameter-matched MLP。
4. biological interpretation 下降。
5. 主要提升来自 known-K 或 oracle settings。
```

停止后不得再添加新模块挽救 CAAM；只能另立新项目。

## sealed test 使用规则

只有当路线 A 或 B 的 protocol 完全冻结后，才允许运行 sealed test：

```text
1. 不改变 corruption_type。
2. 不改变 n_top_genes。
3. 不改变 mask_ratio。
4. 不改变 model architecture。
5. 不改变 evaluation protocol。
```

sealed test 只能运行一次主协议；失败后不得回头调参。

## 输出

新增：

```text
BDD/CAAM_scMAE_BDD_plan/CAAM_PUBLICATION_DECISION.md
```

必须包含：

```text
chosen_route: method_paper | protocol_analysis | stop
main_evidence
negative_evidence
required_next_experiments
forbidden_claims
journal_target_tier
```

## 禁止行为

```text
不以工程 smoke 作为论文有效性证据
不以单个数据集结果作为主张
不把 no_positive_interaction 美化为 synergy
不把 validation/test 用作调参
不为了论文叙事删除负结果
```
