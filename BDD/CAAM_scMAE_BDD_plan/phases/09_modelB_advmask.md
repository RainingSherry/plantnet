# Phase 9: Model B AdvMask-scMAE

## 目标

验证 adversarial mask 是否优于 random mask。Model B 只替换 mask selector，不更换 encoder。

## 模型定义

```text
Model B = MLP encoder + adversarial mask selector + matched gene-wise donor
```

禁止：

```text
AxialEncoder
gene-axis attention
cell-axis attention
NeighborMix
clustering loss
contrastive loss
ZINB/NB loss
free value generator
```

## 操作步骤

1. 实现 `trainers/advmask_trainer.py`。
2. Stage 1 warm-up：

```text
student_warmup_epochs 内使用 random fixed-budget mask
只更新 student
不更新 generator
```

3. Stage 2 alternating training：

```text
每个 batch 执行 student step
每 generator_update_interval 个 batch 执行 generator step
```

4. Student step：

```text
generator 不更新
student 更新
mask 可 no_grad 或 detach
replacement V detach
```

5. Generator step：

```text
student parameters requires_grad=False
student eval mode
generator train mode
使用 mask_st
student forward 不得放入 torch.no_grad()
```

6. Generator step corruption 公式：

```python
x_tilde = x * (1 - mask_st) + value.detach() * mask_st
```

不得使用：

```python
torch.where(mask_hard.bool(), value, x)
```

7. Generator loss：

```text
L_phi =
  - L_rec_mask
  - beta * L_mask
  + lambda_coverage * R_coverage
  + lambda_distortion * R_distortion
  + lambda_entropy * R_entropy
```

代码注释必须说明 negative reconstruction/mask loss means generator maximizes student difficulty。

8. Generator step 后恢复 student 原始 train/eval 状态与 requires_grad。

## 验收条件

```text
student step: student_grad_norm > 0
student step: generator_grad_norm == 0
generator step: generator_grad_norm > 0
generator step: student_grad_norm == 0
optimizer parameter id 无交集
真实 generator 参数梯度非零
```

