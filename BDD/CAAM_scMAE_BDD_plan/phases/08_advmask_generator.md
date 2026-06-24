# Phase 8: AdvMask generator 与 relaxed top-k

## 目标

实现受约束的 adversarial mask selector。Generator 只选择 mask 位置，不生成 replacement value。

## Generator 操作步骤

1. 实现 `mask_generator/adversarial_mask.py`。
2. 输入 clean expression：

```text
X.shape = [B,G]
```

3. 结构：

```text
cell_context = CellMLP(X)             -> [B,d_g]
gene_embedding                       -> [G,d_g]
value_feature = ValueProjection(X)   -> [B,G,d_g]
mask_logits = ScoreMLP(tanh(c_i + e_j + v_ij))
```

4. Generator 不得接收：

```text
labels
n_clusters
batch_code
library_size
zero_ratio
context labels
cluster predictions
```

metadata 只允许 donor corruption 使用。

5. Eligibility-aware logits：

```python
masked_logits = logits.masked_fill(~eligibility, -large_value)
```

不得用随机重采样绕开 eligibility。

## Relaxed top-k 操作步骤

1. 实现 `mask_generator/relaxed_topk.py`。
2. 返回：

```text
mask_hard
mask_soft
mask_st
```

3. 满足：

```text
mask_hard in {0,1}
mask_hard.sum(dim=1) = k_i
mask_soft in [0,1]
mask_soft.sum(dim=1) approx k_i
mask_st forward value = mask_hard
mask_st backward gradient = mask_soft gradient
```

4. 即使启用 AMP，Gumbel noise 与 relaxed top-k 必须用 float32。

## Regularizer 操作步骤

实现 `mask_generator/regularizers.py`：

```text
coverage KL vs eligibility target
distortion target window
entropy regularizer
mask entropy
mask Gini
per-gene mask rate
```

## 验收条件

```text
generator 不输出 replacement value
eligibility 位置约束生效
relaxed top-k gradient test 通过
regularizer shape 与数值稳定
```

