# Phase 14: Correction 3 — AdvMask triage

## 目标

在 Phase 13 选出的最佳或前两名 corruption 上，判断 AdvMask 是否真的优于 random mask。本阶段不引入 Axial，也不运行 full 模型。

## previous_phase_check

执行前必须检查：

```text
1. corruption_triad_summary.csv 是否存在。
2. PHASE13_CORRUPTION_TRIAD_REPORT.md 是否存在。
3. 是否完成 development datasets 上的 control/random-mask corruption 对比，或明确记录失败原因。
4. 是否推荐了 1-2 个 corruption_type 进入本阶段。
5. 是否没有使用 validation/sealed test。
6. 是否没有注册 corruption variants。
```

若 1-4 不满足，不得继续。

## 实验设计

Datasets：

```text
data/processed/Quake_Smart-seq2_Lung.h5ad
data/其他/Mouse_Pancreas_1.h5ad
data/processed_scmae/Limb_Muscle.h5ad
```

Feature space：

```text
input_mode = log1p
n_top_genes = 2000
scale_input = false
```

Models：

```text
control = MLP + random mask
advmask = MLP + AdvMask selector
```

Corruption：使用 Phase 13 最优 corruption_type；如果前两名接近，则两个都跑。

Seeds：

```text
42, 2024, 3407
```

Epochs：

```text
3 epochs first
20 epochs only if the 3-epoch results show a positive trend
```

## 输出

```text
advmask_triage_summary.csv
advmask_triage_report.json
methods/DeepLearning/CAAM_scMAE/benchmark/PHASE14_ADVMASK_TRIAGE_REPORT.md
```

报告必须包含：

```text
1. 每个 dataset/corruption 的 control mean±std。
2. 每个 dataset/corruption 的 advmask mean±std。
3. advmask_minus_control for ARI/NMI/ACC/F1。
4. generator_grad_norm 是否真实大于 0。
5. mask entropy、mask Gini、top-gene concentration。
6. effective_corruption_rate 是否因 AdvMask 改变。
7. embedding collapse diagnostics。
8. 是否建议保留 AdvMask。
```

## 判断规则

AdvMask 进入下一阶段的最低条件：

```text
1. 至少 2/3 development datasets 上 ARI 的 advmask_minus_control > 0。
2. 平均提升大于 seed std 的 0.5 倍。
3. mask 不集中到少数 genes。
4. generator_grad_norm 真实非零。
5. embedding 没有 collapse。
6. biological interpretation 不明显下降。
```

若不满足，AdvMask 删除或降级为 supplementary negative result；不得运行 Axial/full 来替代这个判断。

## 测试要求

新增或更新：

```text
test_advmask_triage.py
```

最低测试：

```text
1. AdvMask triage runner 只运行 control/advmask。
2. summary 正确计算 advmask_minus_control。
3. generator grad norm 字段存在。
4. runner 不写入 formal benchmark manifest。
```

## 禁止行为

```text
不跑 axial/full
不改 manifest
不使用 validation/sealed test
不根据单个数据集调 loss
不把 3-epoch 趋势写成最终结论
不提交 results/
```
