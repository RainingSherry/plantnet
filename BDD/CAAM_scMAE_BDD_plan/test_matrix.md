# CAAM-scMAE test matrix

## Phase 1: package and config

必须测试：

```text
配置默认值完整
CLI 覆盖 YAML
resolved_config.yaml 可复现解析结果
variant 只能是 control/axial/advmask/full
```

## Phase 2: preprocessing and leakage

必须测试：

```text
benchmark_mode 不重复 preprocessing
standalone raw count 执行 normalize_total -> log1p -> HVG
scale_input 默认 false
training Dataset 不返回 label/cell_type/true_cluster/n_clusters
label 不进入训练调用链
```

## Phase 3: donor and corruption

必须测试：

```text
donor r != i
replacement 来自同一 gene
不同 gene 独立 donor
donor selection 不使用 label
ineligible 位置不被 mask
budget deficit 被记录
deficit 超阈值 fail-fast
```

## Phase 4: random mask and losses

必须测试：

```text
fixed budget 下每个 cell 的 mask 数等于可行预算
masked reconstruction loss 按 mask.sum 归一化
visible reconstruction loss 单独归一化
mask prediction 使用 BCEWithLogits
```

## Phase 5: Model 0

必须测试：

```text
Model 0 shape test
Model 0 不实例化 AxialEncoder
Model 0 不实例化 AdversarialMaskGenerator
student step 有梯度
embedding_final.npy shape 正确
```

## Phase 6: gene modules and context selection

必须测试：

```text
gene_module_ids 同 seed 可复现
同一数据集不同训练 seed 共享 gene_module_ids
context_indices 同 seed 可复现
context selection 不使用 label
```

## Phase 7: Model A

必须测试：

```text
tokenizer: [B,G] -> [B,M,d]
gene attention: [B,H,M,M]
cell attention: [B,M,H,C]
self-exclusion 后 query i -> context i attention weight = 0
context_cache_checksum 同 seed 一致
Model A 不实例化 adversarial/difficulty generator
```

## Phase 8: relaxed top-k

必须测试：

```text
mask_hard 为 0/1
mask_hard.sum(dim=1) = k_i
mask_soft 在 [0,1]
mask_st forward value 等于 mask_hard
mask_st backward gradient 来自 mask_soft
float32 下计算 Gumbel 和 relaxed top-k
```

## Phase 9: Model B gradient flow

必须测试：

```text
student step: student_grad_norm > 0
student step: generator_grad_norm == 0
generator step: generator_grad_norm > 0
generator step: student_grad_norm == 0
generator step 不用 torch.where hard mask
generator step student forward 不在 torch.no_grad 中
两个 optimizer parameter id 无交集
```

## Phase 10: Model C and ablation

必须测试：

```text
Model C 包含 AxialEncoder 与 AdversarialMaskGenerator
Model C 不包含新 loss 或 BDD 禁止模块
2x2 factorial 输出 Y00/Y10/Y01/Y11
interaction delta = Y11 - Y10 - Y01 + Y00
```

## Phase 11: formal benchmark

必须测试：

```text
methods/method_manifest.yaml 只注册正式 caam_scmae
caam_scmae 调用 --variant full
artifact_manifest 完整
embedding_final.npy 存在且 shape 正确
正式 eval 输出可生成
```

## Reproducibility hard test

同 seed 两次短训练，必须一致或在数值容差内一致：

```text
context_indices
gene_module_ids
first_batch_indices
first_mask_hard
first_donor_indices
first_loss
context_cache_checksum
```

