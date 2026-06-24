# Phase 5: Model 0 Controlled-scMAE

## 目标

实现 CAAM 的受控基线，作为 Model A/B/C 因果对照。

## 模型定义

```text
Model 0 = MLP encoder + fixed-budget random mask + matched gene-wise donor
```

共享：

```text
Decoder
MaskHead
Loss
Mask conditioning
Student optimizer
Preprocessing
Mask budget
Training epochs
Clustering protocol
```

## 操作步骤

1. 实现 `models/caam_model.py` 中 control variant 构建。
2. 实现 `trainers/control_trainer.py`。
3. 每个 batch 执行：

```text
sample donor V and eligibility E
random fixed-budget mask M_hard from eligible positions
X_tilde = (1-M_hard)X + M_hard V
MLP encoder -> Z
MaskHead/Decoder -> mask_logits/X_hat
student loss backward
```

4. 实现 clean inference：

```text
clean X -> Encoder -> Z
```

5. 按原始 index 写回 embedding。
6. 保存：

```text
embedding_final.npy
embeddings_base.npy
training_history.json
model_checkpoint_last.pt
artifact_manifest.json
mask_stats.json
gradient_stats.json
embedding_stats.json
```

## 验收条件

```text
Model 0 可短训练
Model 0 不实例化 AxialEncoder
Model 0 不实例化 AdversarialMaskGenerator
shape/no-label/donor/budget/gradient tests 通过
```

