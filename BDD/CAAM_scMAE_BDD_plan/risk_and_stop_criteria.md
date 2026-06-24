# CAAM-scMAE risks and stop criteria

## 最高风险点

1. Generator 梯度被截断。
2. Student 与 generator 串梯度。
3. Donor replacement 退化成整行 donor 或跨 gene value。
4. Context set 产生 self-copy 泄漏。
5. 标签通过 preprocessing、subsampling、context、gene modules、donor、early stopping 或 resolution selection 泄漏。
6. Corruption 暴露 single-value shortcut。
7. Mask generator 坍缩到少数 genes 或 HVG。
8. Attention 学到 batch/library size/zero ratio 等技术因素。
9. Axial 提升来自参数量而非结构。
10. 旧 artifact 污染 benchmark skip 逻辑。

## Generator gradient 硬约束

Generator step 必须使用：

```python
x_tilde = x * (1 - mask_st) + value.detach() * mask_st
```

不得退回：

```python
torch.where(mask_hard.bool(), value, x)
```

Generator step 中 student 参数必须 freeze，但 student forward 不得放入 `torch.no_grad()`，否则 loss 无法回传到 mask selector。

必须检查真实参数梯度：

```text
generator_grad_norm > 0
student_grad_norm == 0
```

不能只检查：

```python
mask.requires_grad
```

## Context cache 与复现性硬约束

必须满足：

```text
context_indices 在训练开始前固定并保存
每个 epoch 只刷新 context tokens，不重新选择 context cells
refresh context cache 时使用 model.eval() + torch.no_grad()
context key/value detach
query cell 如果在 context_indices 中，必须 self-exclusion
DataLoader 使用固定 torch.Generator
第一版 num_workers=0
embedding extraction 按原始 index 写回
checkpoint 保存 RNG states、context_indices、gene_module_ids
```

## Label leakage 核查

全仓库和 CAAM 包内重点搜索：

```text
label
labels
cell_type
n_clusters
true_cluster
```

这些字段不得进入：

```text
training Dataset
mask generator
corruption
encoder
loss
trainer
prototype/context selection
gene module builder
donor pool
early stopping
resolution selection
```

## Fail-fast 条件

训练时可使用的 label-free fail-fast：

```text
loss NaN/Inf
generator gradient 为 0
student gradient 为 0
effective mask deficit cells > 1%
normalized mask entropy < 0.2 持续 5 次检查
top 10% genes 获得 >80% mask 持续 5 次检查
mask Gini > 0.85 持续 5 次检查
effective rank < 0.3 * latent_dim
mean pairwise cosine > 0.95
per-dimension mean variance < 1e-4
generator grad norm > 10 * student grad norm 持续 3 次
```

Fail-fast 不得使用：

```text
ARI
NMI
真实类别
rare-cell label
```

## 止损条件

出现以下情况时停止继续加模块：

```text
A/B 都没有单独效果
C 没有正交互项
generator 持续坍缩
shortcut probe 很高
attention 强依赖 batch
参数匹配 MLP 达到同样结果
提升小于 seed 波动
```

## Synergy 声称条件

定义：

```text
Y00 = Model 0
Y10 = Model A
Y01 = Model B
Y11 = Model C
```

交互项：

```text
Delta_AB = Y11 - Y10 - Y01 + Y00
```

只有同时满足：

```text
Model C > Model A
Model C > Model B
Delta_AB > 0
paired confidence interval of Delta_AB excludes 0
```

才允许声称 axial encoder 与 adversarial masking 存在 synergy。

