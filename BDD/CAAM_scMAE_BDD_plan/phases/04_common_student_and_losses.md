# Phase 4: 共同 student、random mask 与 loss

## 目标

实现四个模型共享的 student 组件，避免后续模型差异混入 decoder、mask head 或 loss。

## 操作步骤

1. 实现 `mask_generator/random_mask.py`。
2. fixed budget 下每个 cell 的 mask 数：

```text
k_i = min(round(mask_ratio * G), eligible_count_i)
```

3. 实现 `models/mlp_encoder.py`。
4. 实现共同：

```text
models/decoder.py
models/mask_head.py
losses/reconstruction.py
losses/mask_prediction.py
losses/loss_bundle.py
```

5. 默认 decoder mask conditioning：

```text
pred_detached
```

公式：

```text
P_M = sigmoid(mask_logits)
X_hat = Decoder([Z, stopgrad(P_M)])
```

6. reconstruction loss 必须分 masked 与 visible：

```text
L_rec = L_rec_mask + lambda_visible * L_rec_visible
```

7. mask loss 使用：

```python
binary_cross_entropy_with_logits(mask_logits, mask_hard)
```

## 验收条件

```text
四个模型共享 Decoder/MaskHead/Loss/Mask conditioning
random mask 满足每细胞预算
loss 不对全矩阵直接 mean
BCE 不提前 sigmoid
```

