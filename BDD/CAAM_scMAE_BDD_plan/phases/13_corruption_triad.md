# Phase 13: Correction 2 — Corruption triad

## 目标

一次性实现并比较三种 corruption（扰动）机制，避免继续默认坚持 matched donor。

三种 corruption：

```text
A. scmae_shuffle: scMAE-style gene-wise shuffle（基因列内随机打乱）
B. matched_donor: matched donor shuffle（匹配供体打乱）
C. nonzero_aware_donor: nonzero-aware donor shuffle（非零/变化感知供体打乱）
```

本阶段只评估 corruption 本身，因此第一轮只能使用：

```text
MLP encoder + random mask
```

不得运行 Axial，不得运行 full。

## previous_phase_check

执行本阶段前，Codex 必须检查 Phase 12：

```text
1. benchmark_mode 默认是否已改为 n_top_genes=2000。
2. strict_effective_budget 是否默认为 false。
3. zero_to_zero_rate/effective_corruption_rate/budget_deficit_rate 是否已写入 artifact。
4. validate_formal_smoke.py 是否不再硬编码 n_top_genes==0。
5. Phase 12 smoke 是否通过；若未通过，不得跑本阶段实验。
6. method_manifest.yaml 是否未加入 corruption variants。
7. 是否没有提交 results/ 或 data/smoke/*.h5ad。
```

若 1-5 未满足，不得继续。

## Corruption A: scMAE-style gene-wise shuffle

定义：

```text
对每个 gene column 在 cell 维度随机 permutation。
V_ij = X_{perm_j(i), j}
```

约束：

```text
不匹配 batch/library_size/zero_ratio
不使用 label
不保证 V_ij != X_ij
zero-to-zero 只记录为 diagnostic
```

必须保存：

```text
per_gene_permutation_seed
zero_to_zero_rate
effective_corruption_rate
mean_abs_delta
```

## Corruption B: matched donor shuffle

定义：

```text
V_ij = X_{r_ij, j}
r_ij 来自与 cell i 匹配的 donor pool
匹配条件包含 batch_code、library_size bin、zero_ratio bin
```

约束：

```text
r_ij != i
replacement value 来自同一 gene
不同 gene 可以独立采样 donor
fallback levels 必须记录
不使用 label
```

必须保存：

```text
donor_fallback_matched
donor_fallback_batch
donor_fallback_global
zero_to_zero_rate
effective_corruption_rate
mean_abs_delta
```

## Corruption C: nonzero-aware donor shuffle

定义：

```text
优先从 donor candidates 中选择能让 replacement value 与 original value 发生有效变化的位置。
```

最低实现：

```text
change-aware mode:
candidate r 满足 abs(X_rj - X_ij) > tolerance
```

可选增强：

```text
zero-to-nonzero mode:
若 X_ij == 0，优先选择 X_rj > 0
```

必须有 fallback：

```text
1. change-aware donor
2. matched donor
3. scMAE-style gene-wise shuffle
```

不得因为找不到 nonzero donor 而直接失败，除非 strict_effective_budget=true。

必须记录：

```text
nonzero_aware_success_rate
fallback_to_matched_rate
fallback_to_scmae_shuffle_rate
zero_to_zero_rate
effective_corruption_rate
mean_abs_delta
```

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
variant = control
encoder = MLP
mask_selector = random
```

Corruption：

```text
scmae_shuffle
matched_donor
nonzero_aware_donor
```

Seeds：

```text
42, 2024, 3407
```

Epochs：

```text
3 first
10 optional if all 3 datasets run cleanly
```

## 输出

每个 run 输出：

```text
metrics.json
corruption_stats.json
artifact_manifest.json
resolved_config.yaml
embedding_final.npy
```

汇总输出：

```text
results/CAAM_scMAE_correction/corruption_triad_summary.csv
results/CAAM_scMAE_correction/corruption_triad_report.json
```

报告写入：

```text
methods/DeepLearning/CAAM_scMAE/benchmark/PHASE13_CORRUPTION_TRIAD_REPORT.md
```

报告必须包含：

```text
1. 每个 dataset/corruption 的 ARI/NMI/ACC/F1 mean±std
2. zero_to_zero_rate mean±std
3. effective_corruption_rate mean±std
4. mean_abs_delta mean±std
5. scmae_shuffle_minus_matched_donor
6. nonzero_aware_minus_scmae_shuffle
7. nonzero_aware 是否只提高 mask prediction 但不提高 clustering
8. 推荐进入 Phase 14 的 corruption_type
```

## 判断规则

```text
1. 若 matched_donor 在 2/3 datasets 上弱于 scmae_shuffle，则 matched_donor 不得作为主 corruption。
2. 若 nonzero_aware_donor 提升 clustering 且没有明显 shortcut/collapse，则进入 Phase 14。
3. 若三种 corruption 均差异小于 seed 波动，则选择最简单的 scmae_shuffle 作为后续 baseline。
4. 若某 corruption 只提高 mask prediction，不提高 clustering，则不得写成主贡献。
```

## 测试要求

新增或更新：

```text
test_corruption_types.py
```

最少测试：

```text
1. scmae_shuffle 保持每个 gene 的边缘分布。
2. matched_donor 不选择 self donor。
3. nonzero_aware 优先选择 changed donor。
4. fallback 不使用 label。
5. 三种 corruption 都输出 corruption_stats.json 所需字段。
```

## 禁止行为

```text
不跑 axial/full
不改 manifest
不注册 corruption variants
不使用 validation/sealed test
不根据单个 dataset 的 ARI 自动改代码
不提交 results/
```
