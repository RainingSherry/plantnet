# CAAM-scMAE execution order

本文件只描述阶段顺序。每个阶段的详细步骤见 `phases/`。

## Phase 1: 独立包与 CLI

创建 `methods/DeepLearning/CAAM_scMAE/` 独立包，建立配置、registry、`run.py`、目录骨架与基础 artifact 规范。

验收重点：

```text
run.py 支持 --variant control/axial/advmask/full
配置优先级为 默认值 < YAML < CLI
resolved_config.yaml 可保存
不创建仓库顶层重复 configs/data/tests/scripts
```

## Phase 2: 数据与 preprocessing

实现 canonical benchmark input 与 standalone preprocessing，并明确训练 Dataset 不返回 label。

验收重点：

```text
benchmark_mode 不重新 normalize/log1p/HVG/scale
standalone raw count 执行 normalize_total -> log1p -> HVG
scale_input 默认 false
training DataLoader 不含 label/cell_type/true_cluster/n_clusters
```

## Phase 3: Donor 与 corruption 基础

实现 DonorCandidateProvider、eligibility、matched gene-wise donor replacement、corruption 输出契约。

验收重点：

```text
donor r != i
replacement value 来自同一 gene
不同 masked gene 独立采样 donor
selector 只能选择 eligible 位置
budget deficit 可记录并可 fail-fast
```

## Phase 4: 共同 student、random mask 与 loss

实现随机固定预算 mask、MLP encoder、共同 Decoder、MaskHead、mask conditioning 与 student loss。

验收重点：

```text
masked/visible reconstruction 按位置数归一化
mask prediction 使用 BCEWithLogits
decoder/mask head/loss 在四个模型中共享
```

## Phase 5: Model 0 Controlled-scMAE

完成控制基线训练、推理、embedding 输出、基础 artifacts 和 Model 0 测试。

验收重点：

```text
MLP encoder + random fixed-budget mask + matched gene-wise donor
无 AxialEncoder
无 AdversarialMaskGenerator
embedding_final.npy shape = [N, latent_dim]
```

## Phase 6: Gene modules 与 fixed context

实现 label-free gene module builder 与 fixed context selection。

验收重点：

```text
gene_module_ids 在同一数据集不同 seed 共享
context_indices 在训练开始前固定并保存
context cells 不基于标签选择
```

## Phase 7: Model A Axial-scMAE

实现 module tokenizer、gene-axis attention、cell-axis context attention、context cache、self-exclusion 与 Axial encoder。

验收重点：

```text
gene-axis 不跨 cell 混合
cell-axis 使用 fixed context set
query cell 不能读取自身 clean context
Model A 不导入 adversarial/difficulty generator
```

## Phase 8: AdvMask generator 与 relaxed top-k

实现 adversarial mask generator、eligibility-aware logits、relaxed top-k、temperature schedule 与 generator regularizers。

验收重点：

```text
generator 只输出 mask logits/mask
不得输出 replacement value
relaxed top-k hard forward / soft backward
Gumbel 和 relaxed top-k 用 float32
```

## Phase 9: Model B AdvMask-scMAE

实现 alternating training 中的 generator step，并验证真实 generator 参数梯度。

验收重点：

```text
student step: student_grad_norm > 0, generator_grad_norm == 0
generator step: generator_grad_norm > 0, student_grad_norm == 0
generator step 不使用 torch.no_grad 包住 student forward
generator step 使用 x * (1-mask_st) + value.detach() * mask_st
```

## Phase 10: Model C 与内部消融

组合 Axial encoder 与 AdvMask generator，不引入新机制，完成 2x2 factorial 内部消融 runner。

验收重点：

```text
Model C 仅为 A+B
内部 ablation 输出到 results/CAAM_scMAE_ablation/
计算 Y11 - Y10 - Y01 + Y00
只有满足交互条件才称 synergy
```

## Phase 11: 正式 benchmark 接入

只注册 `caam_scmae`，对应 `--variant full`。

验收重点：

```text
methods/method_manifest.yaml 只注册 caam_scmae
envs/runtime_registry.yaml 注册或复用 CAAM runtime
scripts/run_formal_benchmark.py 能调用 CAAM run.py
正式输出目录只出现 caam_scmae
```

