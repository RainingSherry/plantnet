# Phase 10: Model C 与内部 ablation

## 目标

组合 Model A 与 Model B，验证是否存在 factorial interaction。Model C 不允许加入 A/B 之外的新机制。

## 模型定义

```text
Model C = Axial encoder + adversarial mask selector + matched gene-wise donor
```

禁止额外加入：

```text
NeighborMix
contrastive loss
clustering loss
ZINB/NB loss
second generator
free value generator
new decoder type
```

## 操作步骤

1. 实现 `trainers/full_trainer.py`。
2. 复用 Model A 的：

```text
gene modules
fixed context set
context cache
self-exclusion
AxialEncoder
```

3. 复用 Model B 的：

```text
AdversarialMaskGenerator
relaxed top-k
generator regularizers
alternating training
gradient isolation
```

4. 实现 `benchmark/run_ablation.py`。
5. 内部消融运行：

```text
Model 0: control
Model A: axial
Model B: advmask
Model C: full
parameter-matched MLP
```

6. 输出到：

```text
results/CAAM_scMAE_ablation/
```

7. 统计：

```text
Y00 = Model 0
Y10 = Model A
Y01 = Model B
Y11 = Model C
Delta_AB = Y11 - Y10 - Y01 + Y00
```

8. 对 paired dataset-seed observations 计算：

```text
mean
std
median
win/tie/loss
paired bootstrap 95% CI
interaction paired bootstrap 95% CI
```

## 验收条件

```text
Model C 仅为 A+B
内部 ablation 不污染正式 benchmark 主表
计算 C-A, C-B, C-baseline, Delta_AB
只有满足 BDD 条件才称 synergy
```

