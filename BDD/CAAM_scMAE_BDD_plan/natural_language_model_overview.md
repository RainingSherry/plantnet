# CAAM-scMAE natural language model overview

本文档用自然语言说明 CAAM-scMAE 的任务、目标、模型框架、核心公式、理论支撑、关键风险，以及必须由研究者本人判定的问题。它不替代 BDD 中的工程约束，而是帮助开发者理解为什么这些约束必须存在。

执行优先级说明：如果本文档与 `addendum_formal_benchmark_interface.md`、`benchmark_contract.md`、`test_matrix.md` 或 `risk_and_stop_criteria.md` 发生冲突，以那些工程约束文件为准。本文档只负责解释研究思想和判断边界。

---

## 1. 模型任务：CAAM-scMAE 在解决什么问题

CAAM-scMAE 面向单细胞 RNA-seq 聚类中的自监督表示学习。输入是一张细胞 × 基因表达矩阵：

```text
X ∈ R^{N × G}
N = cells 数量
G = genes 数量
X_ij = 第 i 个细胞中第 j 个基因的表达值
```

训练阶段不使用真实细胞类型标签。模型最终输出每个细胞的 embedding：

```text
Z ∈ R^{N × d}
```

然后将 `Z` 交给统一 benchmark 的聚类评测流程，例如 known-K KMeans 或 fixed-resolution Leiden。

CAAM-scMAE 不是分类模型，不直接预测真实 cell type；也不是端到端监督聚类模型。它的核心任务是：通过 masked autoencoding 学到更适合聚类的 cell embedding。

一句话任务定义：

```text
在不使用真实标签的条件下，通过更难、更受控、更不容易被 shortcut 解决的 mask-reconstruction 任务，学习更有聚类结构的单细胞表示。
```

---

## 2. 研究目标：为什么不是简单改 scMAE

原始 scMAE 的核心思路是：随机扰动一部分表达值，让 autoencoder 重建原始表达，并用 encoder 表示进行聚类。这个方向合理，但存在两个关键不足。

第一，普通 MLP encoder 基本把每个细胞的一整行表达向量直接压成 embedding，没有显式区分 gene axis 与 cell axis。scRNA-seq 数据天然是二维表格结构：行是细胞，列是基因。只用 MLP 会弱化“基因模块之间关系”和“细胞群体上下文”的建模。

第二，普通随机 mask 未必有足够训练价值。随机选中的位置可能太容易恢复，也可能太无意义；模型可能学到的是扰动痕迹、零值模式、表达量边际分布，而不是真正的上下文结构。

CAAM-scMAE 的目标是同时解决这两个问题：

```text
Axial encoder：显式建模 gene-axis 与 cell-axis 上下文。
AdvMask selector：选择更难但仍然合法的 mask 位置。
Matched gene-wise donor corruption：让 replacement value 来自真实数据分布，避免人造噪声和单点 shortcut。
```

因此，CAAM-scMAE 不是“把 GAN、attention、scMAE 拼起来”。它更准确的研究定位是：

```text
Constrained adversarial masking + axial context encoder for masked single-cell representation learning.
```

---

## 3. 四模型关系：为什么必须有 Model 0/A/B/C

内部必须实现四个模型，它们构成严格的 2×2 因子设计：

```text
Model 0: Controlled-scMAE = MLP encoder + random mask
Model A: Axial-scMAE      = Axial encoder + random mask
Model B: AdvMask-scMAE    = MLP encoder + adversarial mask
Model C: CAAM-scMAE       = Axial encoder + adversarial mask
```

这四个模型分别回答不同问题。

Model 0 是受控基线。它不是原始 scMAE 的简单复刻，而是在 CAAM 的共同 corruption、共同 decoder、共同 mask head、共同 loss 下构造的公平起点。

Model A 只回答 Axial encoder 是否有贡献。

Model B 只回答 AdvMask selector 是否有贡献。

Model C 是最终方法，只允许等于 Model A 与 Model B 的组合，不得引入额外 loss、额外 clustering head、额外 contrastive loss 或其他未声明机制。

正式 benchmark 主表只注册：

```text
caam_scmae = Model C
```

Model 0/A/B 只能用于内部 ablation，不能进入正式主方法列表。否则论文叙事会从“一个方法与 baselines 比较”混乱成“四个方法同时参赛”。

---

## 4. 总体模型框架

CAAM-scMAE 的训练流程可以自然语言理解为以下步骤。

第一步，读取表达矩阵 `X`。在 benchmark 模式下，输入应被视为 canonical benchmark input，不能重复 normalize/log1p/HVG/scale。standalone 模式下，如果输入是 raw count，才执行 normalize_total、log1p、HVG 等预处理。

第二步，基于 label-free 信息构建 donor pool 和 eligibility。donor pool 可以使用 batch、library size、zero ratio 等技术变量做匹配，但不能使用 cell type label。

第三步，mask selector 选择要扰动的位置。Model 0/A 使用 random fixed-budget mask；Model B/C 使用 adversarial mask selector。

第四步，从 matched gene-wise donor 中取 replacement value，构造 corrupted input。

第五步，student encoder 将 corrupted input 编码成 cell embedding。

第六步，mask head 预测哪些位置被 mask；decoder 重建原始表达。

第七步，训练结束后，推理阶段只输入 clean X，使用 encoder 输出 embedding_final.npy。正式用于聚类评测的是 embedding，不是 decoder 输出，也不是 mask head 输出。

整体数据流：

```text
X
  -> eligibility E
  -> donor replacement V
  -> mask M
  -> corrupted input X_tilde
  -> encoder
  -> cell embedding Z
  -> decoder / mask head during training
  -> embedding_final.npy during inference
```

---

## 5. 核心公式：数据、donor、mask、corruption

表达矩阵记为：

```text
X ∈ R^{N × G}
```

其中 `N` 是细胞数，`G` 是基因数。

mask 记为：

```text
M ∈ {0,1}^{N × G}
```

`M_ij = 1` 表示第 i 个细胞的第 j 个基因被 mask 或替换。

每个 cell 的固定预算为：

```text
k = round(mask_ratio × G)
```

但如果某个 cell 的 eligible 位置不足，则使用：

```text
k_i = min(k, number_of_eligible_positions_in_cell_i)
```

eligibility 的含义是：只有 replacement value 与原值确实不同、且可合法替换的位置才允许被 mask。

matched gene-wise donor replacement 定义为：

```text
V_ij = X_{r_ij, j},  where r_ij ≠ i
```

含义是：第 i 个细胞第 j 个基因的替换值，来自另一个 donor cell 的同一个基因 j。

注意三条硬规则：

```text
1. donor cell 不能等于自己：r_ij ≠ i
2. replacement 必须来自同一个 gene：不能跨 gene 取值
3. 每个 (cell, gene) 独立采 donor：不能整行 donor
```

student step 中可以使用 hard mask 构造 corrupted input：

```text
X_tilde = (1 - M_hard) ⊙ X + M_hard ⊙ V
```

其中 `⊙` 表示逐元素乘法。

Generator step 中必须使用 straight-through mask：

```text
M_st = M_hard + M_soft - stopgrad(M_soft)
```

前向传播时 `M_st` 的数值等于 `M_hard`，反向传播时梯度从 `M_soft` 回到 generator。

Generator step 的 corrupted input 必须写成连续形式：

```text
X_tilde = X ⊙ (1 - M_st) + stopgrad(V) ⊙ M_st
```

不能写成 hard `torch.where(mask_hard.bool(), V, X)`，因为这样会截断 generator 的梯度。

---

## 6. 核心公式：student loss

student 的目标是：根据 corrupted input 重建原始表达，并预测哪些位置被 mask。

masked reconstruction loss：

```text
L_rec_mask = sum_{i,j} M_ij (X_ij - X_hat_ij)^2 / (sum_{i,j} M_ij + eps)
```

visible reconstruction loss：

```text
L_rec_visible = sum_{i,j} (1 - M_ij) (X_ij - X_hat_ij)^2 / (sum_{i,j} (1 - M_ij) + eps)
```

总 reconstruction loss：

```text
L_rec = L_rec_mask + lambda_visible × L_rec_visible
```

mask prediction loss 使用 BCEWithLogitsLoss：

```text
L_mask = BCEWithLogits(mask_logits, M_hard)
```

student 总损失：

```text
L_student = L_rec + lambda_mask × L_mask
```

这里有两个关键点。

第一，masked 和 visible reconstruction 必须分别归一化，不能简单对全矩阵求平均。否则 mask ratio 改变时 loss 尺度会变化。

第二，mask head 输出必须是 logits，不能提前 sigmoid 后再传给 BCEWithLogitsLoss。

---

## 7. 核心公式：generator loss

generator 的任务不是重建表达，而是选择更难但合法的 mask 位置。因此它的优化方向与 student 部分相反。

Generator loss 可以理解为：

```text
L_generator
  = - L_rec_mask
    - beta_mask × L_mask
    + lambda_coverage × R_coverage
    + lambda_distortion × R_distortion
    + lambda_entropy × R_entropy
```

其中：

```text
- L_rec_mask
```

鼓励 generator 选择让 student 更难重建的位置。

```text
- beta_mask × L_mask
```

鼓励 generator 选择让 mask head 更难识别的位置。

coverage regularizer 防止 generator 总是 mask 少数几个基因。

```text
R_coverage = KL(normalized_mask_frequency || normalized_eligibility_frequency)
```

直观含义：mask 的基因分布不能严重偏离 eligible 位置的整体分布。

distortion regularizer 防止 replacement 过弱或过强。

```text
D_ij = |V_ij - X_ij| / (std_gene_j + eps)
```

如果平均 distortion 太小，mask 没训练价值；如果太大，模型可能只学异常值痕迹。

entropy regularizer 防止 mask selector 坍缩到少数位置。

```text
R_entropy = - normalized_entropy(mask_distribution)
```

这里的负号含义是：优化时鼓励更高 entropy，避免过早坍缩。

---

## 8. 核心公式：Axial encoder

Axial encoder 的输入不是完整的 `[B, G, d]` gene token 矩阵，因为那会导致显存过高。CAAM 使用 gene modules 压缩基因维度。

设 gene module assignment 为：

```text
A ∈ {0,1}^{G × M}
```

其中 `M` 是 gene module 数量。归一化后的 assignment 为：

```text
A_bar_{jm} = A_{jm} / (sum_j A_{jm} + eps)
```

module-level expression 为：

```text
U = X × A_bar
U ∈ R^{N × M}
```

也就是把 `G` 个基因压缩成 `M` 个 gene modules。

每个 cell 的 module token 为：

```text
T_{i,m} = ValueMLP(U_{i,m}) + ModuleEmbedding_m
```

Gene-axis attention 在同一个细胞内部的 modules 之间做 attention：

```text
T_gene = Attention_over_modules(T)
```

Cell-axis attention 让当前 query cell 访问 fixed context cells：

```text
T_cell = Attention(query = current cell tokens,
                   key/value = fixed context cell tokens)
```

最终 cell embedding 可以通过 module pooling 得到：

```text
Z_i = Projection(mean_m T_cell_{i,m})
```

关键限制：

```text
1. gene-axis attention 不能跨 cell 混合。
2. cell-axis attention 只能访问 fixed context set。
3. context set 必须 label-free。
4. query cell 如果在 context set 中，必须 self-exclusion，不能读取自己的 clean context。
```

---

## 9. 核心公式：2×2 因子设计与交互项

定义四个模型在某个指标上的表现为：

```text
Y00 = Model 0 = MLP + random mask
Y10 = Model A = Axial + random mask
Y01 = Model B = MLP + adversarial mask
Y11 = Model C = Axial + adversarial mask
```

Axial 的单独贡献：

```text
Effect_A = Y10 - Y00
```

AdvMask 的单独贡献：

```text
Effect_B = Y01 - Y00
```

二者组合后的交互项：

```text
Delta_AB = Y11 - Y10 - Y01 + Y00
```

只有同时满足：

```text
Y11 > Y10
Y11 > Y01
Delta_AB > 0
Delta_AB 的 paired confidence interval 不跨 0
```

才允许声称 Axial encoder 与 AdvMask selector 存在 synergy。

如果 Model C 只是比 Model 0 好，但不比 Model A 或 Model B 好，就不能说组合机制成功。那最多说明某一个模块有效，或者提升来自参数量、训练噪声、数据偶然性。

---

## 10. 理论支撑：为什么这个方法有道理

### 10.1 Masked autoencoding 的理论直觉

Masked autoencoding 的核心思想是：遮住输入的一部分，让模型根据上下文恢复它。若 mask 设计合理，模型不能只记住输入，而必须学习变量之间的依赖关系。

在 scRNA-seq 中，这种依赖关系可能体现为：

```text
同一细胞内基因程序之间的关系
同一细胞群体中表达模式相似性的关系
细胞状态连续变化中的局部结构
```

### 10.2 Matched donor corruption 的理论直觉

如果 replacement value 是自由生成的随机噪声，模型可能很容易发现异常值。这样训练任务会退化成“找噪声”，而不是“理解上下文”。

Matched gene-wise donor corruption 让替换值来自真实数据中同一个 gene 的表达分布，因此更接近数据流形。模型若想识别或修复扰动，就必须利用上下文，而不是只靠单点数值异常。

### 10.3 AdvMask 的理论直觉

随机 mask 可能太容易，也可能太无效。AdvMask 试图选择更难的位置，类似 hard example mining 或 adversarial data selection。

但 CAAM 的 adversarial 是受约束的：generator 只能选择位置，不能生成值。这一点非常重要。否则 generator 会变成自由攻击器，制造非生物学、非数据分布的扰动。

### 10.4 Axial encoder 的理论直觉

表达矩阵是二维结构。Axial encoder 分别建模 gene axis 和 cell axis，可以比单纯 MLP 更明确地利用表格结构。

Gene-axis attention 学习同一个细胞内 gene modules 的依赖。

Cell-axis attention 学习当前细胞相对于固定 context cells 的群体位置。

### 10.5 Fixed context 与 self-exclusion 的理论直觉

固定 context set 让模型拥有跨细胞参照系，但如果 query cell 可以读取自己的 clean context，就会偷看答案。因此 self-exclusion 是防止数据泄漏的核心约束，不是工程优化。

### 10.6 因子实验的理论直觉

Model 0/A/B/C 的 2×2 设计用于回答“到底是哪一部分有效”。如果没有这个设计，审稿人无法判断提升来自 Axial、AdvMask、参数量、训练步数，还是其他 confounder。

---

## 11. Generator gradient 是硬风险

AdvMask 是否真的成立，取决于 generator 是否能从 loss 收到真实梯度。

不能只检查：

```python
mask.requires_grad
```

必须检查真实 generator 参数：

```text
generator_grad_norm > 0
```

Generator step 中 student 参数要 freeze，但 student forward 不能放入 `torch.no_grad()`。原因是 loss 仍然需要通过 student 计算图回传到 `X_tilde`，再回到 `mask_st` 和 generator。

错误逻辑：

```text
student 不更新，所以 student forward 应该 no_grad
```

正确逻辑：

```text
student 参数不更新，但 loss 对 X_tilde 的梯度必须保留
```

否则 Model B/C 的 AdvMask 会表面训练、实际无梯度。

---

## 12. Context cache 与复现性是硬风险

Axial 模型最容易出现两个问题：不可复现和 self-copy 泄漏。

必须保证：

```text
context_indices 训练前固定并保存
每个 epoch 只刷新 context tokens，不重新选择 context cells
refresh context cache 时使用 eval mode + no_grad
context key/value detach
query cell 如果在 context set 中，必须 self-exclusion
DataLoader 使用固定 torch.Generator
第一版 num_workers=0
embedding extraction 按原始 index 写回
checkpoint 保存 RNG states、context_indices、gene_module_ids
```

同 seed 两次短训练必须检查：

```text
context_indices
gene_module_ids
first_batch_indices
first_mask_hard
first_donor_indices
first_loss
context_cache_checksum
```

如果这个 reproducibility test 不通过，不允许进入大规模 benchmark。

---

## 13. Benchmark 与 ablation 的边界

正式 benchmark 只回答：

```text
CAAM-scMAE 作为最终方法是否优于已有方法？
```

因此正式 benchmark 只注册：

```text
caam_scmae -> --variant full
```

内部 ablation 回答：

```text
Axial 是否有用？
AdvMask 是否有用？
二者是否协同？
参数量是否解释了提升？
```

因此内部 ablation 必须跑：

```text
control
axial
advmask
full
parameter-matched MLP
```

但这些不进入正式主表。

---

## 14. 必须由研究者本人判定的地方

以下问题不能由 Codex 自行决定。如果实现过程中遇到这些问题，Codex 必须写 TODO 并停止，等待研究者判断。

### 14.1 论文最终主张强度

需要研究者判断：

```text
最终论文是主张 CAAM-scMAE 是稳定有效的新方法，
还是主张某些机制只在特定数据集/规模上有效，
或者将结果收束为负结果/机制分析？
```

Codex 不能根据局部实验自动改写论文主张。

### 14.2 主指标与主表叙事

需要研究者判断：

```text
主指标以 ARI、NMI、ACC 还是综合排名为主？
known-K KMeans 是主结果还是 oracle-K 补充？
leiden_fixed 是否作为 unknown-K 主协议呈现？
```

Codex 不能把 known-K 结果称作 fully unsupervised，也不能自己决定主表叙事。

### 14.3 development / validation / test 数据集划分

需要研究者判断：

```text
哪些数据集用于开发调参？
哪些数据集用于验证设计？
哪些数据集必须封存为最终 test？
```

Codex 不得在 test 数据集上根据 ARI/NMI 反复修改模型。

### 14.4 是否继续推进 Model C

需要研究者判断：

```text
如果 Model A 和 Model B 单独效果都弱，是否仍继续 Model C？
如果 C 没有正交互项，是否停止继续堆模块？
如果 parameter-matched MLP 达到同样效果，是否放弃 Axial 叙事？
```

Codex 只能报告证据，不能自动扩大模型。

### 14.5 是否允许调整研究机制

以下变化会改变研究机制，必须由研究者批准：

```text
新增 contrastive loss
新增 clustering loss
新增 ZINB loss
新增 supervised label path
让 generator 生成 replacement value
把 donor 改为整行 donor
让 context set 依赖标签或 pseudo-label
用 ARI/NMI 做 early stopping
把 oracle sweep 放入主表
```

Codex 不得为了提升结果自行添加这些内容。

### 14.6 关键超参数是否改变默认研究设定

以下超参数可以在 development set 上调，但不能由 Codex 根据 test 表现自行改：

```text
mask_ratio
context_size
gene_module_count
lambda_visible
lambda_mask
lambda_coverage
lambda_distortion
lambda_entropy
generator_update_interval
temperature schedule
latent_dim
encoder depth
```

如果需要修改这些参数，必须说明修改原因、预期影响、是否只基于 development/validation，而不是 test。

### 14.7 失败条件与止损

如果出现以下情况，需要研究者决定是否止损：

```text
generator 长期梯度异常
mask 长期坍缩到少数 genes
shortcut probe 显示模型主要识别扰动痕迹
attention 主要学习 batch/library size/zero ratio
不同 seed 波动大于平均提升
大数据集显存成本不可接受
```

Codex 不能通过不断加模块来掩盖这些问题。

### 14.8 正式 benchmark 何时开始

正式 benchmark 只能在以下条件满足后开始：

```text
shape tests 通过
label leakage tests 通过
donor/corruption tests 通过
generator gradient-flow tests 通过
context self-exclusion tests 通过
reproducibility hard test 通过
artifact contract 通过
```

是否进入正式 benchmark，由研究者确认。

---

## 15. Codex 必须停止并询问的情况

Codex 遇到以下情况必须停止，不得自行发明解决方案：

```text
BDD 未定义但会改变模型机制的问题
为了提高指标需要新增 loss 或模块
测试失败但可以通过放宽测试绕过
formal runner 接口与 BDD 冲突
需要使用真实 label 参与训练路径
GPU/环境不确定导致结果不可复现
metrics 格式与 benchmark summary 不兼容
artifact 缺失但想让 run 标记成功
```

---

## 16. 最容易偏离的地方

实现中不得：

```text
让 generator 生成 replacement value
把 donor replacement 写成整行 donor
让 label 进入 training、donor、context、gene module 或 early stopping
把 Model 0/A/B/C 混成一个不可拆模型
generator step 使用 hard torch.where
generator step 用 no_grad 包住 student forward
让 query attend 自己的 clean context
把 known-K 结果写成 fully unsupervised
把 control、axial、advmask 注册进正式 benchmark
根据 test ARI/NMI 改模型
为了补效果新增 BDD 未声明机制
```

---

## 17. 一句话总结

CAAM-scMAE 的核心是：用受约束的 adversarial mask selection 构造更难但合法的自监督任务，用 axial context encoder 显式建模 gene-axis 和 cell-axis 上下文，从而学习更难被 shortcut 解决、更适合单细胞聚类的 cell embedding。
