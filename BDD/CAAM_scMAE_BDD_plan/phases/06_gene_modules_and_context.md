# Phase 6: Gene modules 与 fixed context

## 目标

为 Axial encoder 准备固定、label-free、跨 seed 可复用的 gene modules 与 context set。

## Gene module 操作步骤

1. 实现 `data/gene_modules.py`。
2. 对 canonical expression 副本按 gene 中心化和标准化。
3. 对 `X.T` 执行 TruncatedSVD。
4. 得到 gene embedding。
5. 使用 KMeans：

```text
n_clusters = n_gene_modules
n_init = 20
random_state = module_seed
```

6. 保存：

```text
gene_module_ids.npy
gene_module_assignment.npz
gene_module_builder.json
```

7. Model A/C 在同一数据集所有 seed 中使用完全相同的 gene modules。

## Context 操作步骤

1. 实现 `data/context_selection.py`。
2. 对 cells 执行 PCA。
3. 用 MiniBatchKMeans 生成 `context_size` 个中心。
4. 每个中心选择最近真实细胞。
5. 去重；不足时用 label-free farthest-point 补齐。
6. 使用固定 `context_seed`。
7. 保存：

```text
context_indices.npy
context_selection.json
```

## 验收条件

```text
gene_module_ids 可复现
context_indices 可复现
不使用 label
同一数据集不同训练 seed 共享 gene modules/context set
```

