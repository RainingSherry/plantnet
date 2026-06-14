# 2024-2026 Single-cell Representation Methods Scan

本轮检索聚焦近年顶会/顶刊和公开源码，结论是：单细胞表征方法正在从单纯 AE/GNN 转向 foundation model、知识图、contrastive/triplet 以及可解释的图先验。

## Foundation Models

| 方法 | 来源 | 关键思想 | 源码 |
| --- | --- | --- | --- |
| scGPT | Nature Methods 2024 | gene token + expression token transformer，支持 zero-shot embedding/reference mapping | https://github.com/bowang-lab/scGPT |
| scFoundation | Nature Methods 2024 | 100M 参数，xTrimoGene 架构，50M+ 人类单细胞预训练 | https://github.com/biomap-research/scFoundation |
| scPRINT | Nature Communications 2025 | 50M+ cells 预训练，zero-shot denoising、cell embedding、label prediction、GRN inference | https://github.com/cantinilab/scPRINT |
| TranscriptFormer | Science 2025 / CZI | 跨物种 generative cell atlas，112M cells、12 species | https://github.com/czi-ai/transcriptformer |
| CellFM | Nature Communications 2025 | 100M human cells，RetNet-style large model | paper: https://www.nature.com/articles/s41467-025-59926-5 |

迁移判断：

- 这些模型对人/鼠数据最有价值；当前数据大量是植物基因，human vocab/checkpoint 直接迁移会出现基因不匹配。
- 因此本轮没有强行安装/下载大模型，而是新增安全 scGPT runner，并明确记录 fallback。

## Graph / Prior-knowledge Models

| 方法 | 来源 | 关键思想 | 源码 |
| --- | --- | --- | --- |
| scNET | Nature Methods 2025 | 同时建模 gene-gene PPI 和 cell-cell similarity，交替传播学习 context-specific gene/cell embedding | https://github.com/madilabcode/scNET |
| Cell-GraphCompass | National Science Review 2025 | 把每个 cell 建成 gene graph，用图结构 foundation model 表征单细胞 | https://github.com/epang-ucas/Cell-Graph-Compass |
| scHNTL | Bioinformatics 2025 | high-order neighbors + triplet loss，先图注意力 AE，再增强 separability | https://github.com/SWJTU-ML/scHNTL-code |

迁移到 NeighborMix 的部分：

- scNET / graph-attention：加入 per-edge gate，用 anchor/neighbor latent 和图特征学习边可靠性。
- scHNTL：保留“不能只拉近邻居，也要增强 separability”的判断，但没有直接复刻 TF1/Keras 实现。
- Cell-GraphCompass：当前没有植物 gene graph 先验，暂不复刻大模型；保留“图先验应显式进入结构”的思路。

## Benchmarks / Caution

| 来源 | 发现 |
| --- | --- |
| Biology-driven insights into single-cell foundation models, Genome Biology 2025 | 六个 scFM 在不同任务上没有一个方法稳定全胜，简单模型在资源受限/特定数据上仍常更有效。 |
| Zero-shot evaluations of Geneformer/scGPT, 2025 | zero-shot foundation model 不一定超过简单 HVG/PCA baseline，尤其有 batch/species/domain mismatch 时。 |
| 2026 zero-shot embedding benchmark | Geneformer 和 scGPT 相对较强，但很多场景仍不稳定，scFoundation 等可能低于 HVG baseline。 |

这与本轮 scGPT fallback 和 NeighborMix 结果一致：在植物/跨域数据上，直接套 human foundation model 不应被默认视为强基线。

## 本轮可执行迁移

新增的 `canm_gated_cut_mix` 和 `canm_gated_cut_warm` 是本轮对 graph-attention/gating 思路的迁移：

```text
edge_features = [base_weight, cosine_similarity, distance, mutual_knn, snn]
gate = sigmoid(MLP([z_i, z_j, |z_i-z_j|, z_i*z_j, edge_features]))
mixed_x_i = alpha * x_i + (1-alpha) * sum_j normalized(base_w_ij * gate_ij) * x_j
```

实验表明：门控没有改善当前最好结果，说明单纯让 pseudo reconstruction 训练 gate 会把它重新推向“邻域平滑器”。下一步更合理的是把 gate 的监督信号改成 edge-cut 目标，例如用稳定 pseudo partition 或 contrastive/triplet 直接训练 cross-edge probability。
