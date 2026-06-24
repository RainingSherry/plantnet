# Phase 7: Model A Axial-scMAE

## 目标

验证 Axial encoder 是否优于 MLP encoder。Model A 只更换 encoder，不加入 adversarial mask。

## 模型定义

```text
Model A = Axial encoder + fixed-budget random mask + matched gene-wise donor
```

禁止：

```text
AdversarialMaskGenerator
DifficultyMaskGenerator
NeighborMix
clustering loss
contrastive loss
ZINB/NB loss
free value generator
```

## 操作步骤

1. 实现 `models/module_tokenizer.py`。
2. 使用 module pooling：

```text
U = X A_bar
T_im = ValueMLP(U_im) + ModuleEmbedding_m
```

不得构造完整 `[B,G,d]` token。

3. 实现 `models/gene_axis.py`：

```text
Pre-LN -> MHA -> Residual -> Pre-LN -> FFN -> Residual
```

gene-axis 只在同一细胞内部 module tokens 间 attention。

4. 实现 `models/cell_axis.py`：

```text
Query: 当前 batch cell tokens
Key/Value: fixed context cells 的同 module tokens
```

5. 实现 self-exclusion：

```text
query index == context index 时 attention weight 必须为 0
```

6. 实现 context cache：

```text
每 epoch 刷新一次
clean context expression
model.eval() + torch.no_grad()
context key/value detach
epoch 内固定
```

7. 实现 `models/axial_encoder.py`，输出：

```python
{
    "z": Tensor[B, latent_dim],
    "module_tokens": Tensor[B, M, d],
    "gene_attn": Tensor,
    "cell_attn": Tensor,
}
```

8. 实现 `trainers/axial_trainer.py`。

## 验收条件

```text
token shape 正确
gene attention shape 正确
cell attention shape 正确
self-exclusion test 通过
context_cache_checksum 可复现
Model A 不导入 adversarial/difficulty generator
```

