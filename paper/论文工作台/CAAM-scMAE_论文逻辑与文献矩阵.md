# CAAM-scMAE 论文逻辑与文献矩阵

更新时间：2026-06-26

当前写作主线已经转向 protocol-analysis / diagnostic paper。本文档保留历史 CAAM 构想和路线调整记录；新的 source-to-claim 文献矩阵见：

```text
paper/论文工作台/Protocol_Analysis_Literature_Matrix_2026-06-26.md
```

## 0. 当前论文目标

目标不是简单把 TabPFN、GAN 和 scMAE 拼接在一起，而是围绕单细胞聚类中的一个明确问题构造方法：

> 高稀疏、高维 scRNA-seq 表达矩阵中，如何设计更有信息量且不泄漏标签的 masked reconstruction 任务，使模型学到更适合 clustering 的细胞表示？

### 2026-06-26 Phase14 证据更新

Phase14 在 3 个 development datasets 上比较了 `control=random mask` 与 `advmask=constrained adversarial mask selector`：

```text
corruption_type = scmae_shuffle
seeds = 42, 2024, 3407
epochs = 3
complete runs = 18
```

结论：

```text
AdvMask generator 梯度真实非零；
未观察到 mask concentration 或 embedding collapse；
但 mean ARI delta = 0.002403，小于 0.5 * seed std reference = 0.006518；
Phase14 gate_result = fail；
recommendation = drop_or_downgrade_advmask。
```

因此当前论文不能继续默认写成 “adversarial masking 是主贡献”，也不能启动 “Axial + AdvMask synergy” 叙事。更稳妥的当前路线是 masked autoencoding protocol analysis / diagnostic paper：研究 feature space、corruption semantics、effective corruption diagnostics 和 learned masking failure 对 scRNA-seq clustering 表示学习的影响。

### 2026-06-26 Attention/Phase16 路线决策更新

AdvMask 降级后，又做了一个 development-only attention/context smoke：

```text
datasets = Limb_Muscle, Mouse_Pancreas_1, Quake_Smart-seq2_Lung
seed = 42
epochs = 3
corruption_type = scmae_shuffle
mask_selector = random
variants = standard MLP, current Axial, parameter-matched MLP
```

结果：

```text
mean ARI control = 0.561970
mean ARI axial = 0.191264
mean ARI mlp_parammatched = 0.584531
```

因此当前 Axial 实现也不支持作为主贡献。Phase16 BDD 决策已记录为：

```text
chosen_route = protocol_analysis
positive method paper = fail
current CAAM Axial/AdvMask route = stop as positive method
```

这不是说所有 attention / TabPFN-like 思路都不可能，而是说当前 CAAM 的 Axial 实现不能继续作为证据链。若未来重启 attention，应作为新的机制设计问题，重新设 gate，并与 parameter-matched MLP 比较。

历史拟定方法名（当前证据不支持作为主标题）：

```text
CAAM-scMAE: Context-aware Adversarial Axial Masked Autoencoder for single-cell RNA-seq clustering
```

当前关键词应调整为：

```text
Protocol-aware: 纠正 HVG/full-gene、strict effective budget 和 zero-to-zero 语义。
Corruption-aware: 比较 scMAE-style shuffle、matched donor、nonzero-aware donor。
Diagnostic: 记录 zero_to_zero_rate、effective_corruption_rate、mean_abs_delta、mask collapse 等。
Negative evidence: AdvMask 作为主贡献不被当前 development 证据支持。
Phase16 decision: 正向 CAAM 方法论文失败，protocol-analysis / diagnostic paper 路线保留。
```

## 1. 核心文献矩阵

| 方向 | 代表文献 | 已知结论 | 对本论文的意义 |
|---|---|---|---|
| masked AE for scRNA-seq | scMAE, Bioinformatics 2024 | gene-wise shuffle corruption + mask predictor + weighted MSE 在 15 个真实数据集上表现强，并能识别 rare cell types | 直接基底；必须先复现 scMAE-compatible baseline |
| benchmark | scCluBench, 2025/2026 | 36 个数据集、统一流程；scMAE 在深度模型中鲁棒；scCDCG 等 graph 方法更强；大规模和高稀疏会导致 OOM/NaN 等问题 | 论文评估协议应靠近 scCluBench，不能只在一个数据集讲故事 |
| tabular foundation model | TabPFN, Nature 2025 | 对表格 cell 使用二维结构建模：行内 feature attention 和列内 sample attention；小中型表格上强 | 给 CAAM 的 bi-axial encoder 提供结构启发 |
| single-cell foundation models | Geneformer, scGPT, scFoundation, CellFM | 大模型可学习通用 cell embedding，但 clustering 任务上未必优于任务专用模型 | 本文可以定位为 task-specific clustering representation learner，而不是通用 foundation model |
| graph clustering | scCDCG, scDSC, scGNN, AttentionAE-sc | 图方法能建模细胞关系，但可能过平滑、embedding collapse 或对图构建敏感 | 本文可避免显式图构建，用 context attention 捕获 population signal |
| adversarial / generative masking | GAN/WGAN-GP, adversarial training, DOLORIS sparsity mask | 生成式或对抗式机制必须受约束，否则会产生 shortcut、mode collapse 或不稳定 | adversarial mask 不能只靠 discriminator，必须定义合法扰动集合 |

## 2. 论文主线逻辑

### 2.1 现有问题

scMAE 的强点是 masked reconstruction，但它仍有两个不足：

1. Encoder 多为 MLP，将一个 cell 的整行表达直接压缩为 embedding，未显式区分：

```text
cell 内部 gene-gene relation
population 中同一 gene across cells 的分布关系
```

2. Random mask 不区分位置价值。在高稀疏表达矩阵中，大量 mask 可能落在零值或容易恢复的位置，训练信号弱。

scCluBench 指出真实 scRNA-seq 数据常见高稀疏、高噪声、大规模和多 tissue 异质性，因此方法创新必须同时回答：

```text
表示是否更适合 clustering？
mask 是否真的提供有价值训练信号？
方法是否能在高稀疏数据上稳定运行？
```

### 2.2 本文核心假设

本文假设：

> 对 scRNA-seq clustering 来说，一个好的 masked reconstruction 任务不应只随机破坏表达值，而应优先选择对细胞类型区分、基因相关性和群体结构更敏感的位置；同时 encoder 应显式利用表达矩阵的二维上下文。

### 2.3 方法贡献

建议论文贡献写成三条：

1. 提出 bi-axial context encoder，将 scMAE 的 MLP 编码器替换为轻量级双轴上下文编码器，同时捕获 gene-axis 和 cell-axis 信息。
2. 提出 constrained adversarial mask selector，在固定 mask budget、无标签泄漏、覆盖度受控的条件下学习更有价值的 mask 位置。
3. 在 scCluBench-style 多数据集协议下验证 clustering performance、rare cell discovery、embedding distinguishability、mask diagnostics 和 scalability。

## 3. 方法定义：自然语言

给定表达矩阵 `X`，训练时首先生成 mask `M`。被 mask 的表达值不是简单置零，而是使用 gene-wise shuffle 或 matched donor replacement 生成扰动值。模型输入扰动后的表达矩阵，并输出：

```text
1. cell embedding z_i，用于聚类；
2. reconstructed expression x_hat_i，用于重建损失；
3. predicted mask m_hat_i，用于迫使模型识别哪些位置被破坏。
```

第一阶段先做 scMAE-compatible baseline：

```text
gene-wise shuffle corruption
random mask
MLP encoder
weighted MSE + mask prediction loss
```

第二阶段加入 bi-axial encoder：

```text
gene-axis: 在同一 cell 内建模 gene module 之间关系
cell-axis: 对 context cells 中相同 gene/module 的状态做 attention
```

第三阶段加入 adversarial mask selector：

```text
generator 学习选择更难、更有信息量的位置
student 学习在这些位置被扰动后仍恢复原始表达和稳定聚类表示
```

## 4. 方法定义：LaTeX 公式

表达矩阵：

$$
X\in\mathbb{R}^{N\times G},
$$

其中 \(N\) 为细胞数，\(G\) 为基因数。

gene-wise shuffle corruption：

$$
X'_{ij}=X_{\pi_j(i),j},
$$

其中 \(\pi_j\) 是第 \(j\) 个基因在细胞维度上的随机置换。

mask 后输入：

$$
\widetilde{X}=(1-M)\odot X + M\odot X'.
$$

scMAE-style encoder：

$$
z_i=f_\theta(\widetilde{x}_i).
$$

本文希望从只依赖 cell 内上下文：

$$
p(x_{ij}\mid X_{i,-j})
$$

扩展到同时使用 cell 内上下文和 population 上下文：

$$
p(x_{ij}\mid X_{i,-j},X_{-i,j}).
$$

mask predictor：

$$
\widehat{M}=h_\theta(Z).
$$

weighted reconstruction loss：

$$
\mathcal{L}_{rec}
=
\frac{1}{NG}
\sum_{i=1}^{N}\sum_{j=1}^{G}
\Omega_{ij}
(\widehat{X}_{ij}-X_{ij})^2,
$$

其中：

$$
\Omega_{ij}=1+\lambda M_{ij}.
$$

mask prediction loss：

$$
\mathcal{L}_{mask}
=
-
\frac{1}{NG}
\sum_{i=1}^{N}\sum_{j=1}^{G}
\left[
M_{ij}\log\widehat{M}_{ij}
+(1-M_{ij})\log(1-\widehat{M}_{ij})
\right].
$$

总损失：

$$
\mathcal{L}
=
\mathcal{L}_{rec}
+\gamma\mathcal{L}_{mask}.
$$

受约束 adversarial mask：

$$
\min_\theta \max_{\phi\in\Phi}
\mathcal{L}_{rec}(\theta,\phi)
+\gamma\mathcal{L}_{mask}(\theta,\phi).
$$

合法 mask 空间：

$$
\Phi=
\left\{
M:
\sum_j M_{ij}=k_i,\ 
\text{no label leakage},\ 
\text{gene coverage controlled},\ 
\text{sparsity-aware}
\right\}.
$$

## 5. 方法定义：代码雏形

```python
import torch
import torch.nn as nn
import torch.nn.functional as F


def gene_shuffle_corruption(x: torch.Tensor, mask: torch.Tensor) -> torch.Tensor:
    """Shuffle each gene across cells, then mix original and shuffled values by mask."""
    batch_size, n_genes = x.shape
    shuffled_columns = []
    for gene_idx in range(n_genes):
        perm = torch.randperm(batch_size, device=x.device)
        shuffled_columns.append(x[perm, gene_idx])
    x_prime = torch.stack(shuffled_columns, dim=1)
    return x * (1.0 - mask) + x_prime * mask


class BiAxialEncoder(nn.Module):
    """A lightweight TabPFN-inspired encoder for scRNA-seq matrices."""

    def __init__(self, n_genes: int, n_modules: int = 64, d_model: int = 128, n_heads: int = 4):
        super().__init__()
        self.gene_embed = nn.Embedding(n_genes, d_model)
        self.value_proj = nn.Linear(1, d_model)
        self.module_assign = nn.Parameter(torch.randn(n_genes, n_modules))
        self.gene_attn = nn.MultiheadAttention(d_model, n_heads, batch_first=True)
        self.cell_attn = nn.MultiheadAttention(d_model, n_heads, batch_first=True)
        self.out = nn.Linear(d_model, d_model)

    def forward(self, x: torch.Tensor, context_tokens: torch.Tensor | None = None) -> torch.Tensor:
        batch_size, n_genes = x.shape
        gene_ids = torch.arange(n_genes, device=x.device)

        h = self.value_proj(x[..., None]) + self.gene_embed(gene_ids)[None, :, :]

        assignment = torch.softmax(self.module_assign, dim=1)
        tokens = torch.einsum("bgd,gm->bmd", h, assignment)

        tokens, _ = self.gene_attn(tokens, tokens, tokens)

        if context_tokens is not None:
            tokens, _ = self.cell_attn(tokens, context_tokens, context_tokens)

        return self.out(tokens.mean(dim=1))


class MaskGenerator(nn.Module):
    """Constrained top-k mask selector."""

    def __init__(self, n_genes: int, d_model: int):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(n_genes + d_model, 512),
            nn.ReLU(),
            nn.Linear(512, n_genes),
        )

    def forward(self, x: torch.Tensor, z: torch.Tensor, mask_ratio: float = 0.3):
        logits = self.net(torch.cat([x, z], dim=1))
        k = max(1, int(round(mask_ratio * x.shape[1])))
        selected = logits.topk(k, dim=1).indices
        mask = torch.zeros_like(x)
        mask.scatter_(1, selected, 1.0)
        return mask, logits


def caam_scmae_loss(x, encoder, generator, decoder, mask_head, gamma: float = 1.0):
    with torch.no_grad():
        z_clean = encoder(x)

    mask, _ = generator(x, z_clean)
    x_tilde = gene_shuffle_corruption(x, mask)

    z = encoder(x_tilde)
    x_hat = decoder(z)
    mask_logits = mask_head(z)

    rec_weight = 1.0 + 4.0 * mask
    loss_rec = (rec_weight * (x_hat - x).pow(2)).mean()
    loss_mask = F.binary_cross_entropy_with_logits(mask_logits, mask)

    return loss_rec + gamma * loss_mask
```

## 6. 当前最重要的实验顺序

不要直接跑正式 36 数据集 benchmark。建议顺序：

```text
Stage 1: scMAE-compatible reproduction
  - random mask
  - gene-wise shuffle
  - MLP encoder
  - weighted MSE + mask predictor

Stage 2: encoder ablation
  - MLP encoder
  - gene-axis module encoder
  - bi-axial context encoder

Stage 3: mask strategy ablation
  - random mask
  - observed-only mask
  - entropy-regularized adversarial mask
  - coverage-controlled adversarial mask

Stage 4: formal benchmark
  - selected variants only
  - multiple datasets
  - multiple seeds
  - ACC/NMI/ARI
  - embedding distinguishability
  - marker-overlap annotation
  - runtime and memory
```

## 7. 论文暂定结构

### Title

```text
Context-aware Adversarial Masked Autoencoding for Single-cell RNA-seq Clustering
```

### Abstract 逻辑

1. 背景：scRNA-seq clustering 是细胞类型发现的核心任务。
2. 痛点：高维、高稀疏、噪声强；现有 masked AE 没有充分利用二维上下文，random mask 信息量不稳定。
3. 方法：提出 CAAM-scMAE，结合 bi-axial context encoder 和 constrained adversarial mask selector。
4. 结果：在 scCluBench-style 多数据集上提升 ACC/NMI/ARI，并改善 rare cell type discovery 与 embedding distinguishability。
5. 意义：为单细胞聚类提供 task-specific、可扩展、自监督表示学习框架。

### Introduction 逻辑

```text
P1: scRNA-seq 聚类的重要性。
P2: 传统、AE、图模型、foundation model 的局限。
P3: scMAE 证明 masked reconstruction 有效，但仍有上下文和 mask 信息量不足的问题。
P4: TabPFN 提醒我们表格数据应使用 row/column 二维结构。
P5: 本文提出 CAAM-scMAE。
P6: contributions。
```

### Methods

```text
1. Problem formulation
2. scMAE-compatible masked reconstruction
3. Bi-axial context encoder
4. Constrained adversarial mask selector
5. Optimization
6. Clustering and evaluation
```

### Experiments

```text
1. Datasets and protocols
2. Baselines
3. Overall clustering performance
4. Ablation study
5. Rare cell type discovery
6. Embedding distinguishability
7. Mask diagnostics
8. Scalability and robustness
```

## 8. 需要继续补的内容

1. 为每篇核心文献补完整 BibTeX。
2. 明确目标期刊：Nature Methods / Genome Biology / Briefings in Bioinformatics / Bioinformatics / NAR Genomics and Bioinformatics。
3. 下载或确认 scCluBench 原始代码和数据协议。
4. 明确本文最终是否使用 plant single-cell 数据作为特色实验。
5. 写出第一版 Introduction 中文草稿，再翻译成英文。
6. 等实验结果出来后重写 Abstract，不要提前夸大。

## 9. 参考链接

- scMAE: Fang Z, Zheng R, Li M. [scMAE: a masked autoencoder for single-cell RNA-seq clustering](https://academic.oup.com/bioinformatics/article/40/1/btae020/7564641). Bioinformatics, 2024.
- scCluBench: Xu P et al. [scCluBench: Comprehensive Benchmarking of Clustering Algorithms for Single-Cell RNA Sequencing](https://arxiv.org/abs/2512.02471). arXiv, 2025; AAAI 2026 version available from [AAAI proceedings](https://ojs.aaai.org/index.php/AAAI/article/view/37110/41072).
- TabPFN: Hollmann N et al. [Accurate predictions on small data with a tabular foundation model](https://www.nature.com/articles/s41586-024-08328-6). Nature, 2025.
- Geneformer: Theodoris CV et al. [Transfer learning enables predictions in network biology](https://www.nature.com/articles/s41586-023-06139-9). Nature, 2023.
- scGPT: Cui H et al. [scGPT: toward building a foundation model for single-cell multi-omics using generative AI](https://www.nature.com/articles/s41592-024-02201-0). Nature Methods, 2024.
- scFoundation: Hao M et al. [Large-scale foundation model on single-cell transcriptomics](https://www.nature.com/articles/s41592-024-02305-7). Nature Methods, 2024.
- Foundation model evaluation: Liu T et al. [Evaluating the Utilities of Foundation Models in Single-cell Data Analysis](https://advanced.onlinelibrary.wiley.com/doi/10.1002/advs.202514490). Advanced Science, 2026.
