# Phase 15: Correction 4 — Axial re-entry

## 目标

只有 Phase 14 证明 AdvMask 有稳定正向价值后，才允许 Axial 重新进入实验。本阶段不是默认恢复 CAAM-v1，而是严格检查 Axial 是否值得保留。

## previous_phase_check

执行前必须检查：

```text
1. PHASE14_ADVMASK_TRIAGE_REPORT.md 是否存在。
2. AdvMask 是否满足进入下一阶段的最低条件。
3. 是否明确选定 corruption_type。
4. 是否仍未使用 validation/sealed test 做机制选择。
5. 是否没有把 advmask/control 注册进 formal benchmark。
6. 是否没有根据单个数据集结果调 loss 或 mask ratio。
```

若 AdvMask 未通过，不得执行本阶段。

## 实验设计 A：Axial vs parameter-matched MLP

先比较：

```text
MLP + best corruption + random/AdvMask
Axial + best corruption + random/AdvMask
parameter-matched MLP + best corruption + random/AdvMask
```

Feature space：

```text
input_mode = log1p
n_top_genes = 2000
scale_input = false
```

必须报告：

```text
student_params
encoder_params
decoder_params
runtime_seconds
peak_gpu_memory
ARI/NMI/ACC/F1
embedding collapse diagnostics
```

如果 Axial 没有优于 parameter-matched MLP，不得进入 2x2 factorial。

## 实验设计 B：恢复 2x2 factorial

只有实验 A 通过后，才运行：

```text
Y00 = control
Y10 = axial
Y01 = advmask
Y11 = full
```

使用同一：

```text
corruption_type
n_top_genes
mask_ratio
epochs
seeds
optimizer
feature space
```

计算：

```text
Delta_AB = Y11 - Y10 - Y01 + Y00
full_minus_axial
full_minus_advmask
full_minus_parammatched_mlp
```

## Synergy 条件

只有同时满足：

```text
full > axial
full > advmask
Delta_AB > 0
paired CI of Delta_AB excludes 0
```

才允许写：

```text
candidate synergy supported by internal ablation
```

在 validation/test 之前，不得写 `synergy_confirmed`。

## 输出

```text
axial_reentry_summary.csv
factorial_interaction_report.json
methods/DeepLearning/CAAM_scMAE/benchmark/PHASE15_AXIAL_REENTRY_REPORT.md
```

## 测试要求

新增或更新：

```text
test_axial_reentry.py
```

最低测试：

```text
1. Axial re-entry runner 必须读取 Phase 14 gate 文件或显式 --override_research_gate。
2. parameter-matched MLP 参数差距 <= 5%。
3. Delta_AB 计算正确。
4. 不把 internal variants 写入 method_manifest.yaml。
```

## 禁止行为

```text
不使用 sealed test
不跳过 parameter-matched MLP
不把参数量导致的提升写成结构贡献
不把 candidate_positive_interaction 写成 synergy_confirmed
不提交 results/
```
