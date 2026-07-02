# scMAE 改良论文逐篇结合方案与实现约束

> 生成目的：把 `参考文献/02_整理索引` 中的每一篇论文，转化为面向当前 `methods/DeepLearning/scMAE` 原型的**可实现结合方案**，避免继续进行无目标拼接。
>
> 适用仓库/分支：`RainingSherry/plantnet@wide-research`
>
> 主要输入资料：
> - `methods/DeepLearning/scMAEs/参考文献/00_scMAE改良方法整理总报告.md`
> - `methods/DeepLearning/scMAEs/参考文献/02_整理索引.csv`
> - `methods/DeepLearning/scMAEs/01_PDF论文_按推荐程度排序/`
> - 当前代码：`methods/DeepLearning/scMAE` 与 `methods/DeepLearning/scMAEs`
>
> 本文档不是最终论文方法，而是后续实现和筛选的“设计约束 + 逐篇论文接入说明”。

## 0. 总判断：后续不再做“论文名驱动拼接”，而改为“scMAE 缺口驱动改良”

当前已有实验说明，scMAE baseline 已经稳定，普通增加 loss、prototype、token、graph regularizer 很容易只产生微小扰动甚至退化。NeighborMix 的价值不在于它一定是最终答案，而在于它暴露出 scMAE 的关键缺口：**scMAE 缺少可靠的 cell-cell manifold modeling**。

因此所有论文结合必须遵守以下原则：

1. **scMAE 主体不可丢**：保留表达扰动、mask prediction、masked expression reconstruction、latent embedding 聚类这条主线。
2. **NeighborMix 是主线探针，不是信仰**：只有能提升邻居图可靠性、边界保护、稀有细胞保护、mix policy 的方法，才适合与 NeighborMix 直接结合。
3. **外来论文方法必须落到明确接入点**：不得只写“加入某某 loss”；必须说明它进入 scMAE 的哪一层：输入、mask、encoder、decoder、graph controller、semantic target、clustering head、evaluation diagnostic。
4. **输入结构必须适配 scRNA**：raw counts、log-normalized expression、scaled expression 应分开使用；ZINB/NB/token/diffusion 不能直接使用 scaled expression 当 count。
5. **组合必须分阶段**：先 scMAE warmup，再 teacher/graph/retrieval，再 NeighborMix，再 clustering/semantic target；禁止一开始堆满所有 loss。
6. **任何组合都要证明不冲突**：新增模块必须回答它是否破坏 scMAE reconstruction、是否破坏 NeighborMix 的局部流形信号、是否吞并 rare cell、是否造成过平滑。
7. **结果必须按证据等级管理**：
   - legacy-lightweight：只做灵感筛选；
   - independent-full：可作为候选证据；
   - formal-verified：多 seed、多数据集、可进入论文。
8. **所有候选方法必须至少输出机制诊断**：neighbor purity、edge survival、mixed-cell fraction、boundary entropy、rare-cell recall、embedding variance、cluster mass max/min。
9.请你仅需输出融合方案并在Melanoma_5K上进行一轮烟测即可。
---

## 2. 当前代码应先修正的基础设施问题

### 2.1 区分 legacy lightweight 与 independent full

`common/model.py + common/variant_defs.py` 只能用于灵感筛选，不应作为正式论文结果。每个候选模型若要进入论文，应建立独立文件夹，并明确：
- 来源论文；
- 采用的结构；
- 修改点；
- 与原论文差异；
- 与 scMAE 的接入点；
- 禁止声称的内容。

### 2.2 输入预处理必须按方法类型分流

建议统一数据对象中保存三路输入：

```text
raw_counts:      NB/ZINB/dropout likelihood 使用
log_expr:        scMAE reconstruction / NeighborMix 使用
scaled_expr:     MLP/Transformer encoder 使用
rank_or_token:   BEiT/MaskGIT/ScDiVa/VQ 等离散 token 使用
```

不得让所有方法默认使用 `scale_input=True` 后的矩阵。

### 2.3 Graph / Retrieval 不能只做静态 PCA-KNN

当前静态 PCA-KNN + mean context 不足以代表 TabR/Graphormer/GraphMAE/GNN 类方法。后续图模块必须至少包含：
- EMA teacher embedding 建图；
- mutual-KNN；
- 多轮 refresh edge persistence；
- edge confidence；
- boundary / rare-cell veto；
- 推理阶段是否使用 graph/context 的明确设计。

---

## 3. 逐篇论文结合方案

说明：
- “采用部分”表示该论文中最适合借鉴的结构。
- “接入形式”表示它在 scMAE 中应接入的位置。
- “建议损失”给出可实现公式或损失组合。
- “注意事项”列出最容易犯的错误。
- “代码”来自索引中的 `github_url`，若无则标记未提供。

### 1. ScDiVa: Masked Discrete Diffusion for Joint Modeling of Single-Cell Identity and Expression

- **索引推荐**：rank=1，level=强烈推荐，score=97
- **来源/年份**：2026 / arXiv / ICML-style preprint
- **PDF**：`001_强烈推荐_ScDiVa_Masked_Discrete_Diffusion_for_Joint_Modeling_of_Single-Cell_Identity_and_Expression.pdf`
- **代码**：https://github.com/SindiLab/ScDiVa
- **方法类型**：Mask/Target + Dropout-aware diffusion
- **采用部分**：采用吸收式离散 mask、时间步难度调度和表达 token 预测；不要训练完整重型 diffusion 作为主干。
- **与 scMAE 的最佳结合形式**：在 scMAE 的 corruption 层加入 dropout-like absorbing mask；token 必须按 gene-specific quantile/rank 生成，同时保留 log-expression 重构分支。
- **建议损失/训练方式**：L = L_scMAE_rec + λ_mask CE(mask) + λ_tok CE(token_on_mask) + λ_huber Huber(log_expr_on_mask)。
- **额外注意问题**：不能在 scaled input 上做 count token；先区分 raw counts/log-normalized/scaled 三路输入。
- **重复/备注**：未提供

### 2. MaskFeat: Masked Feature Prediction for Self-Supervised Visual Pre-Training

- **索引推荐**：rank=2，level=强烈推荐，score=95
- **来源/年份**：2022 / CVPR
- **PDF**：`002_强烈推荐_MaskFeat_Masked_Feature_Prediction_for_Self-Supervised_Visual_Pre-Training.pdf`
- **代码**：https://github.com/facebookresearch/SlowFast/tree/main/projects/maskfeat
- **方法类型**：Semantic target
- **采用部分**：采用 masked feature prediction 的思想，而不是照搬图像 patch；目标应从低维统计扩展到 gene module/pathway/rank feature。
- **与 scMAE 的最佳结合形式**：在 scMAE decoder 外加 semantic-target head，预测 masked genes 所属模块的均值、zero-rate、rank profile、pathway score。
- **建议损失/训练方式**：L = L_scMAE + λ_feat SmoothL1(module_feature)；feature target 只在 masked/module-block 区域计算。
- **额外注意问题**：当前只预测 mean/std/zero-rate 太弱；必须改成 gene-set/module 级 target。
- **重复/备注**：未提供

### 3. Improved Techniques for Training Consistency Models

- **索引推荐**：rank=3，level=强烈推荐，score=93
- **来源/年份**：2024 / ICLR
- **PDF**：`003_强烈推荐_Improved_Techniques_for_Training_Consistency_Models.pdf`
- **代码**：https://github.com/openai/consistency_models
- **方法类型**：Robust denoising + EMA teacher
- **采用部分**：采用 Pseudo-Huber、噪声强度调度和相邻噪声级一致性；不建议训练完整生成式 consistency model。
- **与 scMAE 的最佳结合形式**：把 teacher 用于 latent target 和 graph construction；student 输入 corrupted expression，teacher 输入 clean/log-normalized expression。
- **建议损失/训练方式**：L = L_scMAE + λ_cons ||z_s(noisy,t)-sg(z_t(clean,t-1))|| + λ_huber Huber(expr)。
- **额外注意问题**：teacher 不能只做同一 batch latent MSE；应同时输出稳定邻居图。
- **重复/备注**：未提供

### 4. JOAO: Automated Data Augmentations for Graph Contrastive Learning

- **索引推荐**：rank=4，level=强烈推荐，score=92
- **来源/年份**：2021 / ICML
- **PDF**：`004_强烈推荐_JOAO_Automated_Data_Augmentations_for_Graph_Contrastive_Learning.pdf`
- **代码**：https://github.com/Shen-Lab/GraphCL_Automated
- **方法类型**：Policy search for mask/mix
- **采用部分**：采用自适应策略搜索思想，搜索 mask type、mask ratio、mix strength、neighbor refresh interval。
- **与 scMAE 的最佳结合形式**：作为 policy controller 包在 scMAE/NeighborMix 外层；搜索空间必须小：random/dropout/module mask × mix_alpha × graph_refresh。
- **建议损失/训练方式**：L_policy 以验证 ARI proxy、loss plateau、graph purity 或 teacher-student agreement 为 reward；训练主损失仍是 scMAE。
- **额外注意问题**：不要让 RL/AutoAugment 成为主贡献；否则工程复杂且不可解释。
- **重复/备注**：未提供

### 5. TabR: Tabular Deep Learning Meets Nearest Neighbors

- **索引推荐**：rank=5，level=强烈推荐，score=91
- **来源/年份**：2024 / ICLR
- **PDF**：`005_强烈推荐_TabR_Tabular_Deep_Learning_Meets_Nearest_Neighbors.pdf`
- **代码**：https://github.com/yandex-research/tabular-dl-tabr
- **方法类型**：Retrieval / reliable neighbor context
- **采用部分**：采用 nearest-neighbor retrieval，但不简单均值邻居；用检索来产生可靠邻居候选和上下文。
- **与 scMAE 的最佳结合形式**：用 EMA teacher embedding 建 ANN/mutual-KNN；检索到的 top-k 进入 attention context 或控制 NeighborMix 边。
- **建议损失/训练方式**：L = L_scMAE + λ_ret ||z_i - sg(Aggregate(z_N(i)))||；NeighborMix 只用高置信 first reliable neighbor。
- **额外注意问题**：当前静态 PCA-KNN + mean context 太弱，且推理时 context 消失；需保持训练/推理一致或明确只作正则。
- **重复/备注**：未提供

### 6. scVGAE: ZINB-Based Variational Graph Autoencoder for Single-Cell RNA-Seq Imputation

- **索引推荐**：rank=6，level=强烈推荐，score=90
- **来源/年份**：2024 / arXiv
- **PDF**：`006_强烈推荐_scVGAE_ZINB-Based_Variational_Graph_Autoencoder_for_Single-Cell_RNA-Seq_Imputation.pdf`
- **代码**：https://github.com/STOmics/scVGAE
- **方法类型**：Distribution / VAE uncertainty
- **采用部分**：采用 NB/ZINB/variational uncertainty 作为重构和不确定性模块，不建议替换整个 scMAE backbone。
- **与 scMAE 的最佳结合形式**：保留 scMAE encoder；decoder 分出 μ/θ/π 或 latent μ/logσ，边界细胞可根据 posterior variance 降低 mix 强度。
- **建议损失/训练方式**：L = L_scMAE + λ_NB NLL_NB/ZINB(raw_counts|μ,θ,π) + β KL(q(z|x)||p(z))。
- **额外注意问题**：必须使用 raw counts + size factor 分支；在 scaled expression 上做 ZINB 是结构性错误。
- **重复/备注**：未提供

### 7. DinoBloom: A Foundation Model for Generalizable Cell Embeddings in Hematology

- **索引推荐**：rank=7，level=强烈推荐，score=88
- **来源/年份**：2024 / MICCAI / arXiv
- **PDF**：`007_强烈推荐_DinoBloom_A_Foundation_Model_for_Generalizable_Cell_Embeddings_in_Hematology.pdf`
- **代码**：https://github.com/marrlab/DinoBloom
- **方法类型**：Teacher-student self-distillation
- **采用部分**：采用 EMA teacher、center/sharpen 或 BYOL predictor，目标是稳定表示和图，而非外部 foundation weight。
- **与 scMAE 的最佳结合形式**：teacher 输入弱扰动/clean view，student 输入强 mask view；teacher embedding 同时用于 consensus graph。
- **建议损失/训练方式**：L = L_scMAE + λ_distill ||p(z_s)-sg(z_t)|| 或 CE(centered teacher probs)；可加 collapse diagnostics。
- **额外注意问题**：不能只加一个 teacher MSE；需要 warmup、EMA decay schedule、variance/covariance 防塌陷监控。
- **重复/备注**：未提供

### 8. scCello: Cell-ontology guided transcriptome foundation model

- **索引推荐**：rank=8，level=强烈推荐，score=87
- **来源/年份**：2024 / NeurIPS 2024 / arXiv
- **PDF**：`008_强烈推荐_scCello_Cell-ontology_guided_transcriptome_foundation_model.pdf`
- **代码**：https://github.com/DeepGraphLearning/scCello
- **方法类型**：Knowledge / interpretability only
- **采用部分**：采用 ontology/LLM/gene text 作为弱知识、解释和后验分析，不作为无监督聚类主监督。
- **与 scMAE 的最佳结合形式**：可把 cell ontology/marker/pathway 用于高置信簇命名、错误边解释、rare-cell validation；训练默认不用标签语义。
- **建议损失/训练方式**：训练损失不应包含 supervised label CE；可选 L_knowledge 仅用于 gene-module attention bias。
- **额外注意问题**：使用外部模型或指令数据会破坏公平性；普通计算机一区论文应保持无监督主体。
- **重复/备注**：未提供

### 9. LangCell: Language-Cell Pre-training for Cell Identity Understanding

- **索引推荐**：rank=9，level=强烈推荐，score=86
- **来源/年份**：2024 / arXiv
- **PDF**：`009_强烈推荐_LangCell_Language-Cell_Pre-training_for_Cell_Identity_Understanding.pdf`
- **代码**：https://github.com/PharMolix/LangCell
- **方法类型**：Knowledge / interpretability only
- **采用部分**：采用 ontology/LLM/gene text 作为弱知识、解释和后验分析，不作为无监督聚类主监督。
- **与 scMAE 的最佳结合形式**：可把 cell ontology/marker/pathway 用于高置信簇命名、错误边解释、rare-cell validation；训练默认不用标签语义。
- **建议损失/训练方式**：训练损失不应包含 supervised label CE；可选 L_knowledge 仅用于 gene-module attention bias。
- **额外注意问题**：使用外部模型或指令数据会破坏公平性；普通计算机一区论文应保持无监督主体。
- **重复/备注**：未提供

### 10. Celler: A Genomic Language Model for Long-Tailed Single-Cell Annotation

- **索引推荐**：rank=10，level=强烈推荐，score=85
- **来源/年份**：2025 / arXiv
- **PDF**：`010_强烈推荐_Celler_A_Genomic_Language_Model_for_Long-Tailed_Single-Cell_Annotation.pdf`
- **代码**：https://github.com/ckqqqq/PyCeller
- **方法类型**：Rare-cell / long-tail protection
- **采用部分**：采用长尾 reweighting、hard sample mining 和 rare-cluster protection，而不是直接做监督长尾分类。
- **与 scMAE 的最佳结合形式**：在 NeighborMix 中为小簇/低密度但稳定区域降低跨簇 mix；prototype 只对高置信核心细胞启用。
- **建议损失/训练方式**：L = L_scMAE + λ_tail weighted prototype/entropy + λ_boundary consistency；mix_weight_i 按 rare-risk 降低。
- **额外注意问题**：不要用当前 KMeans 伪簇直接定义稀有类；用多轮 teacher consensus + density 估计。
- **重复/备注**：未提供

### 11. scMamba: Scalable Foundation Model for Single-Cell Multi-Omics Integration

- **索引推荐**：rank=11，level=强烈推荐，score=84
- **来源/年份**：2025 / arXiv
- **PDF**：`011_强烈推荐_scMamba_Scalable_Foundation_Model_for_Single-Cell_Multi-Omics_Integration.pdf`
- **代码**：https://github.com/ZhengYu-AI-Lab/scMamba
- **方法类型**：Gene sequence / signal modeling
- **采用部分**：采用长程依赖或多尺度信号建模思想，但基因没有天然顺序，必须先构造 gene module/order。
- **与 scMAE 的最佳结合形式**：把 HVG 按共表达/pathway/PPI module 排列或分组，再做 module-block mask/SSM/temporal-style encoder。
- **建议损失/训练方式**：L = L_scMAE + λ_module SmoothL1(module target) + λ_cons multi-view consistency。
- **额外注意问题**：禁止直接按随机 HVG 顺序套序列模型；否则学到的是索引伪结构。
- **重复/备注**：未提供

### 12. A Survey on Foundation Language Models for Single-cell Biology

- **索引推荐**：rank=12，level=强烈推荐，score=83
- **来源/年份**：2025 / ACL 2025
- **PDF**：`012_强烈推荐_A_Survey_on_Foundation_Language_Models_for_Single-cell_Biology.pdf`
- **代码**：未提供
- **方法类型**：Knowledge / interpretability only
- **采用部分**：采用 ontology/LLM/gene text 作为弱知识、解释和后验分析，不作为无监督聚类主监督。
- **与 scMAE 的最佳结合形式**：可把 cell ontology/marker/pathway 用于高置信簇命名、错误边解释、rare-cell validation；训练默认不用标签语义。
- **建议损失/训练方式**：训练损失不应包含 supervised label CE；可选 L_knowledge 仅用于 gene-module attention bias。
- **额外注意问题**：使用外部模型或指令数据会破坏公平性；普通计算机一区论文应保持无监督主体。
- **重复/备注**：未提供

### 13. Masked Modeling for Single-cell Clustering of scRNA-seq Data

- **索引推荐**：rank=13，level=强烈推荐，score=82
- **来源/年份**：2023 / OpenReview
- **PDF**：`013_强烈推荐_Masked_Modeling_for_Single-cell_Clustering_of_scRNA-seq_Data.pdf`
- **代码**：未提供
- **方法类型**：Single-cell masked clustering objective
- **采用部分**：采用其 single-cell masked modeling 和 clustering-aware target，但不替换 scMAE 主任务。
- **与 scMAE 的最佳结合形式**：先 warmup scMAE，再加入高置信 cluster-aware target；目标是提升 latent cluster separability。
- **建议损失/训练方式**：L = L_scMAE + λ_cluster KL(P||Q) + λ_cons high-confidence pair consistency。
- **额外注意问题**：早期 pseudo label 不可靠；低置信/边界细胞不参与硬 KL。
- **重复/备注**：未提供

### 14. CICL: scRNA-seq Data Clustering by Cluster-aware Iterative Contrastive Learning

- **索引推荐**：rank=14，level=强烈推荐，score=81
- **来源/年份**：2023 / arXiv
- **PDF**：`014_强烈推荐_CICL_scRNA-seq_Data_Clustering_by_Cluster-aware_Iterative_Contrastive_Learning.pdf`
- **代码**：https://github.com/WHY-17/Circle
- **方法类型**：Contrastive / pseudo-positive mining
- **采用部分**：采用高置信正样本挖掘和无/弱负样本对比；避免普通 InfoNCE 把同类或连续状态当负样本。
- **与 scMAE 的最佳结合形式**：正样本来自 teacher consensus neighbors、同 pseudo-core、augmentation views；负样本只用远距离高置信不同簇。
- **建议损失/训练方式**：L = L_scMAE + λ_ins InfoNCE_filtered + λ_cluster swapped/cluster consistency。
- **额外注意问题**：必须 false-negative filtering；稀有细胞不能被大簇拉走。
- **重复/备注**：未提供

### 15. scAGC: Learning Adaptive Cell Graphs with Contrastive Guidance for Single-Cell Clustering

- **索引推荐**：rank=15，level=强烈推荐，score=80
- **来源/年份**：2025 / arXiv
- **PDF**：`015_强烈推荐_scAGC_Learning_Adaptive_Cell_Graphs_with_Contrastive_Guidance_for_Single-Cell_Clustering.pdf`
- **代码**：未提供
- **方法类型**：Graph / neighbor reliability
- **采用部分**：采用图方法中的建图、边置信、edge mask、attention、high-order neighbor；不要直接堆多层 GNN。
- **与 scMAE 的最佳结合形式**：图模块作为 NeighborMix controller：输出 edge_confidence、edge_persistence、boundary_score，控制是否 mix。
- **建议损失/训练方式**：L = L_scMAE + λ_edge BCE(masked edge reconstruction) + λ_graph neighbor consistency + λ_dropedge robustness。
- **额外注意问题**：核心风险是过平滑和错边放大；必须 shallow、residual、DropEdge、EMA graph refresh。
- **重复/备注**：未提供

### 16. Computational Methods for Single-Cell Multi-Omics Integration and Alignment

- **索引推荐**：rank=16，level=强烈推荐，score=79
- **来源/年份**：2022 / Genomics Proteomics Bioinformatics / arXiv
- **PDF**：`016_强烈推荐_Computational_Methods_for_Single-Cell_Multi-Omics_Integration_and_Alignment.pdf`
- **代码**：未提供
- **方法类型**：Batch/multi-view integration
- **采用部分**：采用 batch alignment、multi-view fusion、shared/private latent 的思想作为多数据集鲁棒性模块。
- **与 scMAE 的最佳结合形式**：scMAE 主干输出 shared latent；batch/platform view 只用 adversarial/OT alignment 或 view consistency。
- **建议损失/训练方式**：L = L_scMAE + λ_align MMD/OT/adversarial + λ_view multi-view consistency。
- **额外注意问题**：不要把当前单数据集聚类问题复杂化；作为跨批次扩展和附加实验。
- **重复/备注**：未提供

### 17. Interpretable Deep Learning in Single-Cell Omics

- **索引推荐**：rank=17，level=强烈推荐，score=78
- **来源/年份**：2024 / arXiv / review
- **PDF**：`017_强烈推荐_Interpretable_Deep_Learning_in_Single-Cell_Omics.pdf`
- **代码**：未提供
- **方法类型**：Background / cautious reference
- **采用部分**：仅作为背景或消融参考；不要直接作为主模块。
- **与 scMAE 的最佳结合形式**：先判断它补 scMAE 哪个缺口，再决定是否落到 mask、graph、loss、semantic target 或 boundary controller。
- **建议损失/训练方式**：默认使用 L_scMAE，新增 loss 必须有明确机制诊断。
- **额外注意问题**：没有明确结构适配前不进入 formal benchmark。
- **重复/备注**：未提供

### 18. BEiT: BERT Pre-Training of Image Transformers

- **索引推荐**：rank=18，level=强烈推荐，score=77
- **来源/年份**：2022 / ICLR
- **PDF**：`018_强烈推荐_BEiT_BERT_Pre-Training_of_Image_Transformers.pdf`
- **代码**：https://github.com/microsoft/unilm/tree/master/beit
- **方法类型**：Masked modeling / discriminator target
- **采用部分**：采用 token 化、replaced-expression detection、自引导 mask、条件多阶段 mask，但必须映射到基因表达结构。
- **与 scMAE 的最佳结合形式**：保留 scMAE corruption；加 gene-specific rank tokens、mask difficulty schedule、ELECTRA/VIME 判别头。
- **建议损失/训练方式**：L = L_scMAE_rec + λ_mask CE(replaced) + λ_tok CE(rank_token) + λ_curriculum adaptive mask。
- **额外注意问题**：避免图像/文本 token 机械迁移；gene token 应按 gene-specific distribution 或 module 生成。
- **重复/备注**：未提供

### 19. data2vec: A General Framework for Self-supervised Learning in Speech, Vision and Language

- **索引推荐**：rank=19，level=强烈推荐，score=76
- **来源/年份**：2022 / ICML
- **PDF**：`019_强烈推荐_data2vec_A_General_Framework_for_Self-supervised_Learning_in_Speech_Vision_and_Language.pdf`
- **代码**：https://github.com/facebookresearch/fairseq/tree/main/examples/data2vec
- **方法类型**：Teacher-student self-distillation
- **采用部分**：采用 EMA teacher、center/sharpen 或 BYOL predictor，目标是稳定表示和图，而非外部 foundation weight。
- **与 scMAE 的最佳结合形式**：teacher 输入弱扰动/clean view，student 输入强 mask view；teacher embedding 同时用于 consensus graph。
- **建议损失/训练方式**：L = L_scMAE + λ_distill ||p(z_s)-sg(z_t)|| 或 CE(centered teacher probs)；可加 collapse diagnostics。
- **额外注意问题**：不能只加一个 teacher MSE；需要 warmup、EMA decay schedule、variance/covariance 防塌陷监控。
- **重复/备注**：未提供

### 20. MultiMAE: Multi-modal Multi-task Masked Autoencoders

- **索引推荐**：rank=20，level=强烈推荐，score=75
- **来源/年份**：2022 / ECCV
- **PDF**：`020_强烈推荐_MultiMAE_Multi-modal_Multi-task_Masked_Autoencoders.pdf`
- **代码**：https://github.com/EPFL-VILAB/MultiMAE
- **方法类型**：Masked modeling / discriminator target
- **采用部分**：采用 token 化、replaced-expression detection、自引导 mask、条件多阶段 mask，但必须映射到基因表达结构。
- **与 scMAE 的最佳结合形式**：保留 scMAE corruption；加 gene-specific rank tokens、mask difficulty schedule、ELECTRA/VIME 判别头。
- **建议损失/训练方式**：L = L_scMAE_rec + λ_mask CE(replaced) + λ_tok CE(rank_token) + λ_curriculum adaptive mask。
- **额外注意问题**：避免图像/文本 token 机械迁移；gene token 应按 gene-specific distribution 或 module 生成。
- **重复/备注**：未提供

### 21. AudioMAE: Masked Autoencoders that Listen

- **索引推荐**：rank=21，level=推荐，score=74
- **来源/年份**：2022 / NeurIPS
- **PDF**：`021_推荐_AudioMAE_Masked_Autoencoders_that_Listen.pdf`
- **代码**：https://github.com/facebookresearch/AudioMAE
- **方法类型**：Gene sequence / signal modeling
- **采用部分**：采用长程依赖或多尺度信号建模思想，但基因没有天然顺序，必须先构造 gene module/order。
- **与 scMAE 的最佳结合形式**：把 HVG 按共表达/pathway/PPI module 排列或分组，再做 module-block mask/SSM/temporal-style encoder。
- **建议损失/训练方式**：L = L_scMAE + λ_module SmoothL1(module target) + λ_cons multi-view consistency。
- **额外注意问题**：禁止直接按随机 HVG 顺序套序列模型；否则学到的是索引伪结构。
- **重复/备注**：未提供

### 22. I-JEPA: Self-Supervised Learning from Images with Joint-Embedding Predictive Architecture

- **索引推荐**：rank=22，level=推荐，score=73
- **来源/年份**：2023 / CVPR
- **PDF**：`022_推荐_I-JEPA_Self-Supervised_Learning_from_Images_with_Joint-Embedding_Predictive_Architecture.pdf`
- **代码**：https://github.com/facebookresearch/ijepa
- **方法类型**：Teacher-student self-distillation
- **采用部分**：采用 EMA teacher、center/sharpen 或 BYOL predictor，目标是稳定表示和图，而非外部 foundation weight。
- **与 scMAE 的最佳结合形式**：teacher 输入弱扰动/clean view，student 输入强 mask view；teacher embedding 同时用于 consensus graph。
- **建议损失/训练方式**：L = L_scMAE + λ_distill ||p(z_s)-sg(z_t)|| 或 CE(centered teacher probs)；可加 collapse diagnostics。
- **额外注意问题**：不能只加一个 teacher MSE；需要 warmup、EMA decay schedule、variance/covariance 防塌陷监控。
- **重复/备注**：未提供

### 23. MaskGIT: Masked Generative Image Transformer

- **索引推荐**：rank=23，level=推荐，score=72
- **来源/年份**：2022 / CVPR
- **PDF**：`023_推荐_MaskGIT_Masked_Generative_Image_Transformer.pdf`
- **代码**：https://github.com/google-research/maskgit
- **方法类型**：Masked modeling / discriminator target
- **采用部分**：采用 token 化、replaced-expression detection、自引导 mask、条件多阶段 mask，但必须映射到基因表达结构。
- **与 scMAE 的最佳结合形式**：保留 scMAE corruption；加 gene-specific rank tokens、mask difficulty schedule、ELECTRA/VIME 判别头。
- **建议损失/训练方式**：L = L_scMAE_rec + λ_mask CE(replaced) + λ_tok CE(rank_token) + λ_curriculum adaptive mask。
- **额外注意问题**：避免图像/文本 token 机械迁移；gene token 应按 gene-specific distribution 或 module 生成。
- **重复/备注**：未提供

### 24. Flow Matching for Generative Modeling

- **索引推荐**：rank=24，level=推荐，score=71
- **来源/年份**：2023 / ICLR
- **PDF**：`024_推荐_Flow_Matching_for_Generative_Modeling.pdf`
- **代码**：未提供
- **方法类型**：Generative / denoising auxiliary
- **采用部分**：采用生成模型的噪声调度、denoising target、rare-cell augmentation 思想；不建议训练重型生成主干。
- **与 scMAE 的最佳结合形式**：在 scMAE latent 上做轻量 denoising/consistency，或生成难 corruption；生成样本只用于辅助，不进入正式真实样本评估。
- **建议损失/训练方式**：L = L_scMAE + λ_denoise ||z_clean-z_denoised|| + λ_score/flow latent regularization。
- **额外注意问题**：生成模型容易过重且难验证真实性；必须用 marker/pathway/DEG 保真作为附加评估。
- **重复/备注**：未提供

### 25. Anomaly Transformer

- **索引推荐**：rank=25，level=推荐，score=70
- **来源/年份**：2022 / ICLR
- **PDF**：`025_推荐_Anomaly_Transformer.pdf`
- **代码**：https://github.com/thuml/Anomaly-Transformer
- **方法类型**：Background / cautious reference
- **采用部分**：仅作为背景或消融参考；不要直接作为主模块。
- **与 scMAE 的最佳结合形式**：先判断它补 scMAE 哪个缺口，再决定是否落到 mask、graph、loss、semantic target 或 boundary controller。
- **建议损失/训练方式**：默认使用 L_scMAE，新增 loss 必须有明确机制诊断。
- **额外注意问题**：没有明确结构适配前不进入 formal benchmark。
- **重复/备注**：未提供

### 26. BGRL: Large-Scale Representation Learning on Graphs via Bootstrapping

- **索引推荐**：rank=26，level=推荐，score=69
- **来源/年份**：2022 / ICLR
- **PDF**：`026_推荐_BGRL_Large-Scale_Representation_Learning_on_Graphs_via_Bootstrapping.pdf`
- **代码**：https://github.com/nerdslab/bgrl
- **方法类型**：Graph / neighbor reliability
- **采用部分**：采用图方法中的建图、边置信、edge mask、attention、high-order neighbor；不要直接堆多层 GNN。
- **与 scMAE 的最佳结合形式**：图模块作为 NeighborMix controller：输出 edge_confidence、edge_persistence、boundary_score，控制是否 mix。
- **建议损失/训练方式**：L = L_scMAE + λ_edge BCE(masked edge reconstruction) + λ_graph neighbor consistency + λ_dropedge robustness。
- **额外注意问题**：核心风险是过平滑和错边放大；必须 shallow、residual、DropEdge、EMA graph refresh。
- **重复/备注**：未提供

### 27. Graph Barlow Twins

- **索引推荐**：rank=27，level=推荐，score=68
- **来源/年份**：2022 / Knowledge-Based Systems
- **PDF**：`027_推荐_Graph_Barlow_Twins.pdf`
- **代码**：https://github.com/pbielak/graph-barlow-twins
- **方法类型**：Graph / neighbor reliability
- **采用部分**：采用图方法中的建图、边置信、edge mask、attention、high-order neighbor；不要直接堆多层 GNN。
- **与 scMAE 的最佳结合形式**：图模块作为 NeighborMix controller：输出 edge_confidence、edge_persistence、boundary_score，控制是否 mix。
- **建议损失/训练方式**：L = L_scMAE + λ_edge BCE(masked edge reconstruction) + λ_graph neighbor consistency + λ_dropedge robustness。
- **额外注意问题**：核心风险是过平滑和错边放大；必须 shallow、residual、DropEdge、EMA graph refresh。
- **重复/备注**：未提供

### 28. Graphormer

- **索引推荐**：rank=28，level=推荐，score=67
- **来源/年份**：2021 / NeurIPS
- **PDF**：`028_推荐_Graphormer.pdf`
- **代码**：https://github.com/microsoft/Graphormer
- **方法类型**：Graph / neighbor reliability
- **采用部分**：采用图方法中的建图、边置信、edge mask、attention、high-order neighbor；不要直接堆多层 GNN。
- **与 scMAE 的最佳结合形式**：图模块作为 NeighborMix controller：输出 edge_confidence、edge_persistence、boundary_score，控制是否 mix。
- **建议损失/训练方式**：L = L_scMAE + λ_edge BCE(masked edge reconstruction) + λ_graph neighbor consistency + λ_dropedge robustness。
- **额外注意问题**：核心风险是过平滑和错边放大；必须 shallow、residual、DropEdge、EMA graph refresh。
- **重复/备注**：未提供

### 29. Deep Adaptive Fuzzy Clustering for Evolutionary Unsupervised Representation Learning

- **索引推荐**：rank=29，level=推荐，score=66
- **来源/年份**：2021 / arXiv / pattern recognition
- **PDF**：`029_推荐_Deep_Adaptive_Fuzzy_Clustering_for_Evolutionary_Unsupervised_Representation_Learning.pdf`
- **代码**：未提供
- **方法类型**：Boundary / clustering head
- **采用部分**：采用软聚类、模糊 membership、rough lower/upper approximation 或社区发现作为诊断和后期约束。
- **与 scMAE 的最佳结合形式**：先用 scMAE/teacher graph 得到核心细胞；只对核心细胞加强聚类，边界细胞允许高熵。
- **建议损失/训练方式**：L = L_scMAE + λ_core KL(P||Q) + λ_fuzzy entropy/balance + λ_sep center separation。
- **额外注意问题**：不要早期强行 DEC；known-k 假设需在论文中明确。
- **重复/备注**：未提供

### 30. Fuzzy Rough Sets Based on Fuzzy Quantification

- **索引推荐**：rank=30，level=推荐，score=65
- **来源/年份**：2023 / Fuzzy Sets and Systems
- **PDF**：`030_推荐_Fuzzy_Rough_Sets_Based_on_Fuzzy_Quantification.pdf`
- **代码**：未提供
- **方法类型**：Boundary / clustering head
- **采用部分**：采用软聚类、模糊 membership、rough lower/upper approximation 或社区发现作为诊断和后期约束。
- **与 scMAE 的最佳结合形式**：先用 scMAE/teacher graph 得到核心细胞；只对核心细胞加强聚类，边界细胞允许高熵。
- **建议损失/训练方式**：L = L_scMAE + λ_core KL(P||Q) + λ_fuzzy entropy/balance + λ_sep center separation。
- **额外注意问题**：不要早期强行 DEC；known-k 假设需在论文中明确。
- **重复/备注**：未提供

### 31. scPilot: Large Language Model Reasoning Toward Automated Single-Cell Analysis and Discovery

- **索引推荐**：rank=31，level=可用，score=64
- **来源/年份**：2026 / arXiv / OpenReview
- **PDF**：`031_可用_scPilot_Large_Language_Model_Reasoning_Toward_Automated_Single-Cell_Analysis_and_Discovery.pdf`
- **代码**：https://github.com/HelloWorldLTY/scEval/tree/master/scPilot
- **方法类型**：Knowledge / interpretability only
- **采用部分**：采用 ontology/LLM/gene text 作为弱知识、解释和后验分析，不作为无监督聚类主监督。
- **与 scMAE 的最佳结合形式**：可把 cell ontology/marker/pathway 用于高置信簇命名、错误边解释、rare-cell validation；训练默认不用标签语义。
- **建议损失/训练方式**：训练损失不应包含 supervised label CE；可选 L_knowledge 仅用于 gene-module attention bias。
- **额外注意问题**：使用外部模型或指令数据会破坏公平性；普通计算机一区论文应保持无监督主体。
- **重复/备注**：未提供

### 32. ChatCell: Natural Language Interface for Single-Cell Analysis

- **索引推荐**：rank=32，level=可用，score=63
- **来源/年份**：2024 / arXiv
- **PDF**：`032_可用_ChatCell_Natural_Language_Interface_for_Single-Cell_Analysis.pdf`
- **代码**：https://github.com/zjunlp/ChatCell
- **方法类型**：Knowledge / interpretability only
- **采用部分**：采用 ontology/LLM/gene text 作为弱知识、解释和后验分析，不作为无监督聚类主监督。
- **与 scMAE 的最佳结合形式**：可把 cell ontology/marker/pathway 用于高置信簇命名、错误边解释、rare-cell validation；训练默认不用标签语义。
- **建议损失/训练方式**：训练损失不应包含 supervised label CE；可选 L_knowledge 仅用于 gene-module attention bias。
- **额外注意问题**：使用外部模型或指令数据会破坏公平性；普通计算机一区论文应保持无监督主体。
- **重复/备注**：未提供

### 33. iBOT: Image BERT Pre-Training with Online Tokenizer

- **索引推荐**：rank=33，level=可用，score=62
- **来源/年份**：2022 / ICLR
- **PDF**：`033_可用_iBOT_Image_BERT_Pre-Training_with_Online_Tokenizer.pdf`
- **代码**：https://github.com/bytedance/ibot
- **方法类型**：Teacher-student self-distillation
- **采用部分**：采用 EMA teacher、center/sharpen 或 BYOL predictor，目标是稳定表示和图，而非外部 foundation weight。
- **与 scMAE 的最佳结合形式**：teacher 输入弱扰动/clean view，student 输入强 mask view；teacher embedding 同时用于 consensus graph。
- **建议损失/训练方式**：L = L_scMAE + λ_distill ||p(z_s)-sg(z_t)|| 或 CE(centered teacher probs)；可加 collapse diagnostics。
- **额外注意问题**：不能只加一个 teacher MSE；需要 warmup、EMA decay schedule、variance/covariance 防塌陷监控。
- **重复/备注**：未提供

### 34. Masked Siamese Networks for Label-Efficient Learning

- **索引推荐**：rank=34，level=可用，score=61
- **来源/年份**：2022 / ECCV
- **PDF**：`034_可用_Masked_Siamese_Networks_for_Label-Efficient_Learning.pdf`
- **代码**：https://github.com/facebookresearch/msn
- **方法类型**：Teacher-student self-distillation
- **采用部分**：采用 EMA teacher、center/sharpen 或 BYOL predictor，目标是稳定表示和图，而非外部 foundation weight。
- **与 scMAE 的最佳结合形式**：teacher 输入弱扰动/clean view，student 输入强 mask view；teacher embedding 同时用于 consensus graph。
- **建议损失/训练方式**：L = L_scMAE + λ_distill ||p(z_s)-sg(z_t)|| 或 CE(centered teacher probs)；可加 collapse diagnostics。
- **额外注意问题**：不能只加一个 teacher MSE；需要 warmup、EMA decay schedule、variance/covariance 防塌陷监控。
- **重复/备注**：未提供

### 35. VideoMAE: Masked Autoencoders are Data-Efficient Learners for Self-Supervised Video Pre-Training

- **索引推荐**：rank=35，level=可用，score=60
- **来源/年份**：2022 / NeurIPS
- **PDF**：`035_可用_VideoMAE_Masked_Autoencoders_are_Data-Efficient_Learners_for_Self-Supervised_Video_Pre-Training.pdf`
- **代码**：https://github.com/MCG-NJU/VideoMAE
- **方法类型**：Gene sequence / signal modeling
- **采用部分**：采用长程依赖或多尺度信号建模思想，但基因没有天然顺序，必须先构造 gene module/order。
- **与 scMAE 的最佳结合形式**：把 HVG 按共表达/pathway/PPI module 排列或分组，再做 module-block mask/SSM/temporal-style encoder。
- **建议损失/训练方式**：L = L_scMAE + λ_module SmoothL1(module target) + λ_cons multi-view consistency。
- **额外注意问题**：禁止直接按随机 HVG 顺序套序列模型；否则学到的是索引伪结构。
- **重复/备注**：未提供

### 36. DINOv2: Learning Robust Visual Features without Supervision

- **索引推荐**：rank=36，level=可用，score=59
- **来源/年份**：2023 / arXiv / CV foundation model
- **PDF**：`036_可用_DINOv2_Learning_Robust_Visual_Features_without_Supervision.pdf`
- **代码**：https://github.com/facebookresearch/dinov2
- **方法类型**：Teacher-student self-distillation
- **采用部分**：采用 EMA teacher、center/sharpen 或 BYOL predictor，目标是稳定表示和图，而非外部 foundation weight。
- **与 scMAE 的最佳结合形式**：teacher 输入弱扰动/clean view，student 输入强 mask view；teacher embedding 同时用于 consensus graph。
- **建议损失/训练方式**：L = L_scMAE + λ_distill ||p(z_s)-sg(z_t)|| 或 CE(centered teacher probs)；可加 collapse diagnostics。
- **额外注意问题**：不能只加一个 teacher MSE；需要 warmup、EMA decay schedule、variance/covariance 防塌陷监控。
- **重复/备注**：未提供

### 37. ImageBind: One Embedding Space To Bind Them All

- **索引推荐**：rank=37，level=可用，score=58
- **来源/年份**：2023 / CVPR
- **PDF**：`037_可用_ImageBind_One_Embedding_Space_To_Bind_Them_All.pdf`
- **代码**：https://github.com/facebookresearch/ImageBind
- **方法类型**：Background / cautious reference
- **采用部分**：仅作为背景或消融参考；不要直接作为主模块。
- **与 scMAE 的最佳结合形式**：先判断它补 scMAE 哪个缺口，再决定是否落到 mask、graph、loss、semantic target 或 boundary controller。
- **建议损失/训练方式**：默认使用 L_scMAE，新增 loss 必须有明确机制诊断。
- **额外注意问题**：没有明确结构适配前不进入 formal benchmark。
- **重复/备注**：未提供

### 38. Latent Diffusion Models

- **索引推荐**：rank=38，level=可用，score=57
- **来源/年份**：2022 / CVPR
- **PDF**：`038_可用_Latent_Diffusion_Models.pdf`
- **代码**：https://github.com/CompVis/latent-diffusion
- **方法类型**：Generative / denoising auxiliary
- **采用部分**：采用生成模型的噪声调度、denoising target、rare-cell augmentation 思想；不建议训练重型生成主干。
- **与 scMAE 的最佳结合形式**：在 scMAE latent 上做轻量 denoising/consistency，或生成难 corruption；生成样本只用于辅助，不进入正式真实样本评估。
- **建议损失/训练方式**：L = L_scMAE + λ_denoise ||z_clean-z_denoised|| + λ_score/flow latent regularization。
- **额外注意问题**：生成模型容易过重且难验证真实性；必须用 marker/pathway/DEG 保真作为附加评估。
- **重复/备注**：未提供

### 39. TS2Vec: Towards Universal Representation of Time Series

- **索引推荐**：rank=39，level=可用，score=56
- **来源/年份**：2022 / AAAI
- **PDF**：`039_可用_TS2Vec_Towards_Universal_Representation_of_Time_Series.pdf`
- **代码**：https://github.com/zhihanyue/ts2vec
- **方法类型**：Gene sequence / signal modeling
- **采用部分**：采用长程依赖或多尺度信号建模思想，但基因没有天然顺序，必须先构造 gene module/order。
- **与 scMAE 的最佳结合形式**：把 HVG 按共表达/pathway/PPI module 排列或分组，再做 module-block mask/SSM/temporal-style encoder。
- **建议损失/训练方式**：L = L_scMAE + λ_module SmoothL1(module target) + λ_cons multi-view consistency。
- **额外注意问题**：禁止直接按随机 HVG 顺序套序列模型；否则学到的是索引伪结构。
- **重复/备注**：未提供

### 40. TimesNet: Temporal 2D-Variation Modeling for General Time Series Analysis

- **索引推荐**：rank=40，level=可用，score=55
- **来源/年份**：2023 / ICLR
- **PDF**：`040_可用_TimesNet_Temporal_2D-Variation_Modeling_for_General_Time_Series_Analysis.pdf`
- **代码**：https://github.com/thuml/TimesNet
- **方法类型**：Gene sequence / signal modeling
- **采用部分**：采用长程依赖或多尺度信号建模思想，但基因没有天然顺序，必须先构造 gene module/order。
- **与 scMAE 的最佳结合形式**：把 HVG 按共表达/pathway/PPI module 排列或分组，再做 module-block mask/SSM/temporal-style encoder。
- **建议损失/训练方式**：L = L_scMAE + λ_module SmoothL1(module target) + λ_cons multi-view consistency。
- **额外注意问题**：禁止直接按随机 HVG 顺序套序列模型；否则学到的是索引伪结构。
- **重复/备注**：未提供

### 41. Non-stationary Transformers

- **索引推荐**：rank=41，level=可用，score=54
- **来源/年份**：2022 / NeurIPS
- **PDF**：`041_可用_Non-stationary_Transformers.pdf`
- **代码**：https://github.com/thuml/Nonstationary_Transformers
- **方法类型**：Gene sequence / signal modeling
- **采用部分**：采用长程依赖或多尺度信号建模思想，但基因没有天然顺序，必须先构造 gene module/order。
- **与 scMAE 的最佳结合形式**：把 HVG 按共表达/pathway/PPI module 排列或分组，再做 module-block mask/SSM/temporal-style encoder。
- **建议损失/训练方式**：L = L_scMAE + λ_module SmoothL1(module target) + λ_cons multi-view consistency。
- **额外注意问题**：禁止直接按随机 HVG 顺序套序列模型；否则学到的是索引伪结构。
- **重复/备注**：未提供

### 42. MaskGAE: Masked Graph Modeling Meets Graph Autoencoders

- **索引推荐**：rank=42，level=可用，score=53
- **来源/年份**：2023 / KDD
- **PDF**：`042_可用_MaskGAE_Masked_Graph_Modeling_Meets_Graph_Autoencoders.pdf`
- **代码**：https://github.com/EdisonLeeeee/MaskGAE
- **方法类型**：Graph / neighbor reliability
- **采用部分**：采用图方法中的建图、边置信、edge mask、attention、high-order neighbor；不要直接堆多层 GNN。
- **与 scMAE 的最佳结合形式**：图模块作为 NeighborMix controller：输出 edge_confidence、edge_persistence、boundary_score，控制是否 mix。
- **建议损失/训练方式**：L = L_scMAE + λ_edge BCE(masked edge reconstruction) + λ_graph neighbor consistency + λ_dropedge robustness。
- **额外注意问题**：核心风险是过平滑和错边放大；必须 shallow、residual、DropEdge、EMA graph refresh。
- **重复/备注**：未提供

### 43. GraphGPS: General Powerful Scalable Graph Transformers

- **索引推荐**：rank=43，level=可用，score=52
- **来源/年份**：2022 / NeurIPS
- **PDF**：`043_可用_GraphGPS_General_Powerful_Scalable_Graph_Transformers.pdf`
- **代码**：https://github.com/rampasek/GraphGPS
- **方法类型**：Graph / neighbor reliability
- **采用部分**：采用图方法中的建图、边置信、edge mask、attention、high-order neighbor；不要直接堆多层 GNN。
- **与 scMAE 的最佳结合形式**：图模块作为 NeighborMix controller：输出 edge_confidence、edge_persistence、boundary_score，控制是否 mix。
- **建议损失/训练方式**：L = L_scMAE + λ_edge BCE(masked edge reconstruction) + λ_graph neighbor consistency + λ_dropedge robustness。
- **额外注意问题**：核心风险是过平滑和错边放大；必须 shallow、residual、DropEdge、EMA graph refresh。
- **重复/备注**：未提供

### 44. Mole-BERT: Rethinking Pre-training Graph Neural Networks for Molecules

- **索引推荐**：rank=44，level=可用，score=51
- **来源/年份**：2023 / ICLR
- **PDF**：`044_可用_Mole-BERT_Rethinking_Pre-training_Graph_Neural_Networks_for_Molecules.pdf`
- **代码**：https://github.com/junxia97/Mole-BERT
- **方法类型**：Graph / neighbor reliability
- **采用部分**：采用图方法中的建图、边置信、edge mask、attention、high-order neighbor；不要直接堆多层 GNN。
- **与 scMAE 的最佳结合形式**：图模块作为 NeighborMix controller：输出 edge_confidence、edge_persistence、boundary_score，控制是否 mix。
- **建议损失/训练方式**：L = L_scMAE + λ_edge BCE(masked edge reconstruction) + λ_graph neighbor consistency + λ_dropedge robustness。
- **额外注意问题**：核心风险是过平滑和错边放大；必须 shallow、residual、DropEdge、EMA graph refresh。
- **重复/备注**：未提供

### 45. SAINT: Improved Neural Networks for Tabular Data via Row Attention and Contrastive Pre-training

- **索引推荐**：rank=45，level=可用，score=50
- **来源/年份**：2021 / NeurIPS Workshop
- **PDF**：`045_可用_SAINT_Improved_Neural_Networks_for_Tabular_Data_via_Row_Attention_and_Contrastive_Pre-training.pdf`
- **代码**：https://github.com/somepago/saint
- **方法类型**：Tabular architecture / row attention
- **采用部分**：采用表格学习中的 feature tokenizer、row attention、feature importance；不替换生物 mask 主线。
- **与 scMAE 的最佳结合形式**：在 scMAE encoder 前加入 gene-wise tokenizer/gate，或用 row-neighbor attention 替代静态 context。
- **建议损失/训练方式**：L = L_scMAE + λ_gate sparsity + λ_ret row-neighbor consistency。
- **额外注意问题**：表格方法通常面向监督分类；无监督聚类下需避免标签先验。
- **重复/备注**：未提供

### 46. TabPFN: A Transformer That Solves Small Tabular Classification Problems in a Second

- **索引推荐**：rank=46，level=可用，score=49
- **来源/年份**：2023 / ICLR
- **PDF**：`046_可用_TabPFN_A_Transformer_That_Solves_Small_Tabular_Classification_Problems_in_a_Second.pdf`
- **代码**：https://github.com/automl/TabPFN
- **方法类型**：Tabular architecture / row attention
- **采用部分**：采用表格学习中的 feature tokenizer、row attention、feature importance；不替换生物 mask 主线。
- **与 scMAE 的最佳结合形式**：在 scMAE encoder 前加入 gene-wise tokenizer/gate，或用 row-neighbor attention 替代静态 context。
- **建议损失/训练方式**：L = L_scMAE + λ_gate sparsity + λ_ret row-neighbor consistency。
- **额外注意问题**：表格方法通常面向监督分类；无监督聚类下需避免标签先验。
- **重复/备注**：未提供

### 47. Self-Guided Masked Autoencoder

- **索引推荐**：rank=47，level=可用，score=48
- **来源/年份**：2024 / NeurIPS
- **PDF**：`047_可用_Self-Guided_Masked_Autoencoder.pdf`
- **代码**：https://github.com/cvlab-kaist/Self-Guided-MAE
- **方法类型**：Masked modeling / discriminator target
- **采用部分**：采用 token 化、replaced-expression detection、自引导 mask、条件多阶段 mask，但必须映射到基因表达结构。
- **与 scMAE 的最佳结合形式**：保留 scMAE corruption；加 gene-specific rank tokens、mask difficulty schedule、ELECTRA/VIME 判别头。
- **建议损失/训练方式**：L = L_scMAE_rec + λ_mask CE(replaced) + λ_tok CE(rank_token) + λ_curriculum adaptive mask。
- **额外注意问题**：避免图像/文本 token 机械迁移；gene token 应按 gene-specific distribution 或 module 生成。
- **重复/备注**：未提供

### 48. Masked Autoencoders Are Scalable Vision Learners

- **索引推荐**：rank=48，level=可用，score=47
- **来源/年份**：2022 / CVPR
- **PDF**：`048_可用_Masked_Autoencoders_Are_Scalable_Vision_Learners.pdf`
- **代码**：https://github.com/facebookresearch/mae
- **方法类型**：Masked modeling / discriminator target
- **采用部分**：采用 token 化、replaced-expression detection、自引导 mask、条件多阶段 mask，但必须映射到基因表达结构。
- **与 scMAE 的最佳结合形式**：保留 scMAE corruption；加 gene-specific rank tokens、mask difficulty schedule、ELECTRA/VIME 判别头。
- **建议损失/训练方式**：L = L_scMAE_rec + λ_mask CE(replaced) + λ_tok CE(rank_token) + λ_curriculum adaptive mask。
- **额外注意问题**：避免图像/文本 token 机械迁移；gene token 应按 gene-specific distribution 或 module 生成。
- **重复/备注**：未提供

### 49. Masked Autoencoders Are Scalable Vision Learners

- **索引推荐**：rank=48，level=可用，score=47
- **来源/年份**：CVPR 2022
- **PDF**：`048_可用_Masked_Autoencoders_Are_Scalable_Vision_Learners.pdf`
- **代码**：https://github.com/facebookresearch/mae
- **方法类型**：Masked modeling / discriminator target
- **采用部分**：采用 token 化、replaced-expression detection、自引导 mask、条件多阶段 mask，但必须映射到基因表达结构。
- **与 scMAE 的最佳结合形式**：保留 scMAE corruption；加 gene-specific rank tokens、mask difficulty schedule、ELECTRA/VIME 判别头。
- **建议损失/训练方式**：L = L_scMAE_rec + λ_mask CE(replaced) + λ_tok CE(rank_token) + λ_curriculum adaptive mask。
- **额外注意问题**：避免图像/文本 token 机械迁移；gene token 应按 gene-specific distribution 或 module 生成。
- **重复/备注**：duplicate of rank 48

### 50. Adaptive-Masking Policy with Deep Reinforcement Learning for Self-Supervised Medical Image Segmentation

- **索引推荐**：rank=49，level=可用，score=46
- **来源/年份**：2023 / ICME
- **PDF**：`049_可用_Adaptive-Masking_Policy_with_Deep_Reinforcement_Learning_for_Self-Supervised_Medical_Image_Segmentation.pdf`
- **代码**：未提供
- **方法类型**：Policy search for mask/mix
- **采用部分**：采用自适应策略搜索思想，搜索 mask type、mask ratio、mix strength、neighbor refresh interval。
- **与 scMAE 的最佳结合形式**：作为 policy controller 包在 scMAE/NeighborMix 外层；搜索空间必须小：random/dropout/module mask × mix_alpha × graph_refresh。
- **建议损失/训练方式**：L_policy 以验证 ARI proxy、loss plateau、graph purity 或 teacher-student agreement 为 reward；训练主损失仍是 scMAE。
- **额外注意问题**：不要让 RL/AutoAugment 成为主贡献；否则工程复杂且不可解释。
- **重复/备注**：未提供

### 51. The Dynamic Duo of Collaborative Masking and Target for Advanced Masked Autoencoder Learning

- **索引推荐**：rank=50，level=可用，score=45
- **来源/年份**：2025 / AAAI
- **PDF**：`050_可用_The_Dynamic_Duo_of_Collaborative_Masking_and_Target_for_Advanced_Masked_Autoencoder_Learning.pdf`
- **代码**：未提供
- **方法类型**：Masked modeling / discriminator target
- **采用部分**：采用 token 化、replaced-expression detection、自引导 mask、条件多阶段 mask，但必须映射到基因表达结构。
- **与 scMAE 的最佳结合形式**：保留 scMAE corruption；加 gene-specific rank tokens、mask difficulty schedule、ELECTRA/VIME 判别头。
- **建议损失/训练方式**：L = L_scMAE_rec + λ_mask CE(replaced) + λ_tok CE(rank_token) + λ_curriculum adaptive mask。
- **额外注意问题**：避免图像/文本 token 机械迁移；gene token 应按 gene-specific distribution 或 module 生成。
- **重复/备注**：未提供

### 52. Reinforcement Learning meets Masked Video Modeling: Trajectory-Guided Adaptive Token Selection

- **索引推荐**：rank=51，level=可用，score=44
- **来源/年份**：2025 / arXiv
- **PDF**：`051_可用_Reinforcement_Learning_meets_Masked_Video_Modeling_Trajectory-Guided_Adaptive_Token_Selection.pdf`
- **代码**：未提供
- **方法类型**：Policy search for mask/mix
- **采用部分**：采用自适应策略搜索思想，搜索 mask type、mask ratio、mix strength、neighbor refresh interval。
- **与 scMAE 的最佳结合形式**：作为 policy controller 包在 scMAE/NeighborMix 外层；搜索空间必须小：random/dropout/module mask × mix_alpha × graph_refresh。
- **建议损失/训练方式**：L_policy 以验证 ARI proxy、loss plateau、graph purity 或 teacher-student agreement 为 reward；训练主损失仍是 scMAE。
- **额外注意问题**：不要让 RL/AutoAugment 成为主贡献；否则工程复杂且不可解释。
- **重复/备注**：未提供

### 53. Anatomically-guided Masked Autoencoder with Domain-Adaptive Prompting for multimodal cerebral aneurysm detection and segmentation

- **索引推荐**：rank=52，level=可用，score=43
- **来源/年份**：2026 / npj Digital Medicine
- **PDF**：`052_可用_Anatomically-guided_Masked_Autoencoder_with_Domain-Adaptive_Prompting_for_multimodal_cerebral_aneurysm_detection_and_segmentation.pdf`
- **代码**：未提供
- **方法类型**：Masked modeling / discriminator target
- **采用部分**：采用 token 化、replaced-expression detection、自引导 mask、条件多阶段 mask，但必须映射到基因表达结构。
- **与 scMAE 的最佳结合形式**：保留 scMAE corruption；加 gene-specific rank tokens、mask difficulty schedule、ELECTRA/VIME 判别头。
- **建议损失/训练方式**：L = L_scMAE_rec + λ_mask CE(replaced) + λ_tok CE(rank_token) + λ_curriculum adaptive mask。
- **额外注意问题**：避免图像/文本 token 机械迁移；gene token 应按 gene-specific distribution 或 module 生成。
- **重复/备注**：未提供

### 54. DAP-MAE: Domain-Adaptive Point Cloud Masked Autoencoder for Effective Cross-Domain Learning

- **索引推荐**：rank=53，level=可用，score=42
- **来源/年份**：2025 / arXiv
- **PDF**：`053_可用_DAP-MAE_Domain-Adaptive_Point_Cloud_Masked_Autoencoder_for_Effective_Cross-Domain_Learning.pdf`
- **代码**：未提供
- **方法类型**：Masked modeling / discriminator target
- **采用部分**：采用 token 化、replaced-expression detection、自引导 mask、条件多阶段 mask，但必须映射到基因表达结构。
- **与 scMAE 的最佳结合形式**：保留 scMAE corruption；加 gene-specific rank tokens、mask difficulty schedule、ELECTRA/VIME 判别头。
- **建议损失/训练方式**：L = L_scMAE_rec + λ_mask CE(replaced) + λ_tok CE(rank_token) + λ_curriculum adaptive mask。
- **额外注意问题**：避免图像/文本 token 机械迁移；gene token 应按 gene-specific distribution 或 module 生成。
- **重复/备注**：未提供

### 55. An Empirical Study of Multiple Masking in Masked Autoencoder / Conditional MAE

- **索引推荐**：rank=54，level=可用，score=41
- **来源/年份**：2025 / OpenReview/withdrawn
- **PDF**：`054_可用_An_Empirical_Study_of_Multiple_Masking_in_Masked_Autoencoder_Conditional_MAE.pdf`
- **代码**：未提供
- **方法类型**：Masked modeling / discriminator target
- **采用部分**：采用 token 化、replaced-expression detection、自引导 mask、条件多阶段 mask，但必须映射到基因表达结构。
- **与 scMAE 的最佳结合形式**：保留 scMAE corruption；加 gene-specific rank tokens、mask difficulty schedule、ELECTRA/VIME 判别头。
- **建议损失/训练方式**：L = L_scMAE_rec + λ_mask CE(replaced) + λ_tok CE(rank_token) + λ_curriculum adaptive mask。
- **额外注意问题**：避免图像/文本 token 机械迁移；gene token 应按 gene-specific distribution 或 module 生成。
- **重复/备注**：未提供

### 56. Multi-Facet Clustering Variational Autoencoders

- **索引推荐**：rank=55，level=可用，score=40
- **来源/年份**：2021 / NeurIPS
- **PDF**：`055_可用_Multi-Facet_Clustering_Variational_Autoencoders.pdf`
- **代码**：https://github.com/adrianjav/heterogeneous_vaes
- **方法类型**：Distribution / VAE uncertainty
- **采用部分**：采用 NB/ZINB/variational uncertainty 作为重构和不确定性模块，不建议替换整个 scMAE backbone。
- **与 scMAE 的最佳结合形式**：保留 scMAE encoder；decoder 分出 μ/θ/π 或 latent μ/logσ，边界细胞可根据 posterior variance 降低 mix 强度。
- **建议损失/训练方式**：L = L_scMAE + λ_NB NLL_NB/ZINB(raw_counts|μ,θ,π) + β KL(q(z|x)||p(z))。
- **额外注意问题**：必须使用 raw counts + size factor 分支；在 scaled expression 上做 ZINB 是结构性错误。
- **重复/备注**：未提供

### 57. ScInfoVAE: interpretable dimensional reduction of single cell transcription data with variational autoencoders and extended mutual information

- **索引推荐**：rank=56，level=可用，score=39
- **来源/年份**：2023 / BioData Mining
- **PDF**：`056_可用_ScInfoVAE_interpretable_dimensional_reduction_of_single_cell_transcription_data_with_variational_autoencoders_and_extended_mutual_information.pdf`
- **代码**：未提供
- **方法类型**：Distribution / VAE uncertainty
- **采用部分**：采用 NB/ZINB/variational uncertainty 作为重构和不确定性模块，不建议替换整个 scMAE backbone。
- **与 scMAE 的最佳结合形式**：保留 scMAE encoder；decoder 分出 μ/θ/π 或 latent μ/logσ，边界细胞可根据 posterior variance 降低 mix 强度。
- **建议损失/训练方式**：L = L_scMAE + λ_NB NLL_NB/ZINB(raw_counts|μ,θ,π) + β KL(q(z|x)||p(z))。
- **额外注意问题**：必须使用 raw counts + size factor 分支；在 scaled expression 上做 ZINB 是结构性错误。
- **重复/备注**：未提供

### 58. scVAEDer: integrating deep diffusion models and variational autoencoders for single-cell transcriptomics analysis

- **索引推荐**：rank=57，level=可用，score=38
- **来源/年份**：2025 / Genome Biology
- **PDF**：`057_可用_scVAEDer_integrating_deep_diffusion_models_and_variational_autoencoders_for_single-cell_transcriptomics_analysis.pdf`
- **代码**：https://github.com/zhyu-lab/scVAEDer
- **方法类型**：Distribution / VAE uncertainty
- **采用部分**：采用 NB/ZINB/variational uncertainty 作为重构和不确定性模块，不建议替换整个 scMAE backbone。
- **与 scMAE 的最佳结合形式**：保留 scMAE encoder；decoder 分出 μ/θ/π 或 latent μ/logσ，边界细胞可根据 posterior variance 降低 mix 强度。
- **建议损失/训练方式**：L = L_scMAE + λ_NB NLL_NB/ZINB(raw_counts|μ,θ,π) + β KL(q(z|x)||p(z))。
- **额外注意问题**：必须使用 raw counts + size factor 分支；在 scaled expression 上做 ZINB 是结构性错误。
- **重复/备注**：未提供

### 59. Gene selection for single cell RNA-seq data via fuzzy rough iterative computation model

- **索引推荐**：rank=58，level=可用，score=37
- **来源/年份**：2025 / Artificial Intelligence Review
- **PDF**：`058_可用_Gene_selection_for_single_cell_RNA-seq_data_via_fuzzy_rough_iterative_computation_model.pdf`
- **代码**：未提供
- **方法类型**：Boundary / clustering head
- **采用部分**：采用软聚类、模糊 membership、rough lower/upper approximation 或社区发现作为诊断和后期约束。
- **与 scMAE 的最佳结合形式**：先用 scMAE/teacher graph 得到核心细胞；只对核心细胞加强聚类，边界细胞允许高熵。
- **建议损失/训练方式**：L = L_scMAE + λ_core KL(P||Q) + λ_fuzzy entropy/balance + λ_sep center separation。
- **额外注意问题**：不要早期强行 DEC；known-k 假设需在论文中明确。
- **重复/备注**：未提供

### 60. Soft Graph Clustering for single-cell RNA Sequencing Data

- **索引推荐**：rank=59，level=可用，score=36
- **来源/年份**：2025 / arXiv/BMC Bioinformatics
- **PDF**：`059_可用_Soft_Graph_Clustering_for_single-cell_RNA_Sequencing_Data.pdf`
- **代码**：未提供
- **方法类型**：Graph / neighbor reliability
- **采用部分**：采用图方法中的建图、边置信、edge mask、attention、high-order neighbor；不要直接堆多层 GNN。
- **与 scMAE 的最佳结合形式**：图模块作为 NeighborMix controller：输出 edge_confidence、edge_persistence、boundary_score，控制是否 mix。
- **建议损失/训练方式**：L = L_scMAE + λ_edge BCE(masked edge reconstruction) + λ_graph neighbor consistency + λ_dropedge robustness。
- **额外注意问题**：核心风险是过平滑和错边放大；必须 shallow、residual、DropEdge、EMA graph refresh。
- **重复/备注**：未提供

### 61. q-Diffusion leverages the full dimensionality of gene coexpression in single-cell transcriptomics

- **索引推荐**：rank=60，level=可用，score=35
- **来源/年份**：2024 / Communications Biology
- **PDF**：`060_可用_q-Diffusion_leverages_the_full_dimensionality_of_gene_coexpression_in_single-cell_transcriptomics.pdf`
- **代码**：https://github.com/zsteve/qDiffusion
- **方法类型**：Generative / denoising auxiliary
- **采用部分**：采用生成模型的噪声调度、denoising target、rare-cell augmentation 思想；不建议训练重型生成主干。
- **与 scMAE 的最佳结合形式**：在 scMAE latent 上做轻量 denoising/consistency，或生成难 corruption；生成样本只用于辅助，不进入正式真实样本评估。
- **建议损失/训练方式**：L = L_scMAE + λ_denoise ||z_clean-z_denoised|| + λ_score/flow latent regularization。
- **额外注意问题**：生成模型容易过重且难验证真实性；必须用 marker/pathway/DEG 保真作为附加评估。
- **重复/备注**：未提供

### 62. scDiffusion: conditional generation of high-quality single-cell data using diffusion model

- **索引推荐**：rank=61，level=可用，score=34
- **来源/年份**：2024 / Bioinformatics/arXiv
- **PDF**：`061_可用_scDiffusion_conditional_generation_of_high-quality_single-cell_data_using_diffusion_model.pdf`
- **代码**：https://github.com/EperLuo/scDiffusion
- **方法类型**：Generative / denoising auxiliary
- **采用部分**：采用生成模型的噪声调度、denoising target、rare-cell augmentation 思想；不建议训练重型生成主干。
- **与 scMAE 的最佳结合形式**：在 scMAE latent 上做轻量 denoising/consistency，或生成难 corruption；生成样本只用于辅助，不进入正式真实样本评估。
- **建议损失/训练方式**：L = L_scMAE + λ_denoise ||z_clean-z_denoised|| + λ_score/flow latent regularization。
- **额外注意问题**：生成模型容易过重且难验证真实性；必须用 marker/pathway/DEG 保真作为附加评估。
- **重复/备注**：未提供

### 63. Scalable Single-Cell Gene Expression Generation with Latent Diffusion Models

- **索引推荐**：rank=62，level=可用，score=33
- **来源/年份**：2025/2026 / arXiv
- **PDF**：`062_可用_Scalable_Single-Cell_Gene_Expression_Generation_with_Latent_Diffusion_Models.pdf`
- **代码**：https://github.com/OmicsML/scLDM
- **方法类型**：Generative / denoising auxiliary
- **采用部分**：采用生成模型的噪声调度、denoising target、rare-cell augmentation 思想；不建议训练重型生成主干。
- **与 scMAE 的最佳结合形式**：在 scMAE latent 上做轻量 denoising/consistency，或生成难 corruption；生成样本只用于辅助，不进入正式真实样本评估。
- **建议损失/训练方式**：L = L_scMAE + λ_denoise ||z_clean-z_denoised|| + λ_score/flow latent regularization。
- **额外注意问题**：生成模型容易过重且难验证真实性；必须用 marker/pathway/DEG 保真作为附加评估。
- **重复/备注**：未提供

### 64. Attention-Guided Probabilistic Diffusion Model for Generating Cell-Type-Specific Gene Regulatory Networks from Gene Expression Data

- **索引推荐**：rank=63，level=可用，score=32
- **来源/年份**：2025 / Genes
- **PDF**：`063_可用_Attention-Guided_Probabilistic_Diffusion_Model_for_Generating_Cell-Type-Specific_Gene_Regulatory_Networks_from_Gene_Expression_Data.pdf`
- **代码**：未提供
- **方法类型**：Generative / denoising auxiliary
- **采用部分**：采用生成模型的噪声调度、denoising target、rare-cell augmentation 思想；不建议训练重型生成主干。
- **与 scMAE 的最佳结合形式**：在 scMAE latent 上做轻量 denoising/consistency，或生成难 corruption；生成样本只用于辅助，不进入正式真实样本评估。
- **建议损失/训练方式**：L = L_scMAE + λ_denoise ||z_clean-z_denoised|| + λ_score/flow latent regularization。
- **额外注意问题**：生成模型容易过重且难验证真实性；必须用 marker/pathway/DEG 保真作为附加评估。
- **重复/备注**：未提供

### 65. GRouNdGAN: GRN-guided simulation of single-cell RNA-seq data using causal generative adversarial networks

- **索引推荐**：rank=64，level=可用，score=31
- **来源/年份**：2024 / Nature Communications
- **PDF**：`064_可用_GRouNdGAN_GRN-guided_simulation_of_single-cell_RNA-seq_data_using_causal_generative_adversarial_networks.pdf`
- **代码**：https://github.com/Emad-COMBINE-lab/GRouNdGAN
- **方法类型**：Generative / denoising auxiliary
- **采用部分**：采用生成模型的噪声调度、denoising target、rare-cell augmentation 思想；不建议训练重型生成主干。
- **与 scMAE 的最佳结合形式**：在 scMAE latent 上做轻量 denoising/consistency，或生成难 corruption；生成样本只用于辅助，不进入正式真实样本评估。
- **建议损失/训练方式**：L = L_scMAE + λ_denoise ||z_clean-z_denoised|| + λ_score/flow latent regularization。
- **额外注意问题**：生成模型容易过重且难验证真实性；必须用 marker/pathway/DEG 保真作为附加评估。
- **重复/备注**：未提供

### 66. Realistic in silico generation and augmentation of single-cell RNA-seq data using generative adversarial networks

- **索引推荐**：rank=65，level=可用，score=30
- **来源/年份**：2020 / Nature Communications
- **PDF**：`065_可用_Realistic_in_silico_generation_and_augmentation_of_single-cell_RNA-seq_data_using_generative_adversarial_networks.pdf`
- **代码**：https://github.com/imsb-uke/scGAN
- **方法类型**：Generative / denoising auxiliary
- **采用部分**：采用生成模型的噪声调度、denoising target、rare-cell augmentation 思想；不建议训练重型生成主干。
- **与 scMAE 的最佳结合形式**：在 scMAE latent 上做轻量 denoising/consistency，或生成难 corruption；生成样本只用于辅助，不进入正式真实样本评估。
- **建议损失/训练方式**：L = L_scMAE + λ_denoise ||z_clean-z_denoised|| + λ_score/flow latent regularization。
- **额外注意问题**：生成模型容易过重且难验证真实性；必须用 marker/pathway/DEG 保真作为附加评估。
- **重复/备注**：未提供

### 67. scVGATAE: A Variational Graph Attentional Autoencoder Model for Clustering Single-Cell RNA-seq Data

- **索引推荐**：rank=66，level=可用，score=29
- **来源/年份**：2024 / Biology
- **PDF**：`066_可用_scVGATAE_A_Variational_Graph_Attentional_Autoencoder_Model_for_Clustering_Single-Cell_RNA-seq_Data.pdf`
- **代码**：https://github.com/xkmaxidian/scVGATAE
- **方法类型**：Graph / neighbor reliability
- **采用部分**：采用图方法中的建图、边置信、edge mask、attention、high-order neighbor；不要直接堆多层 GNN。
- **与 scMAE 的最佳结合形式**：图模块作为 NeighborMix controller：输出 edge_confidence、edge_persistence、boundary_score，控制是否 mix。
- **建议损失/训练方式**：L = L_scMAE + λ_edge BCE(masked edge reconstruction) + λ_graph neighbor consistency + λ_dropedge robustness。
- **额外注意问题**：核心风险是过平滑和错边放大；必须 shallow、residual、DropEdge、EMA graph refresh。
- **重复/备注**：未提供

### 68. scSiameseClu: A Siamese Clustering Framework for Interpreting Single-cell RNA Sequencing Data

- **索引推荐**：rank=67，level=可用，score=28
- **来源/年份**：2025 / IJCAI
- **PDF**：`067_可用_scSiameseClu_A_Siamese_Clustering_Framework_for_Interpreting_Single-cell_RNA_Sequencing_Data.pdf`
- **代码**：https://github.com/XPgogogo/scSiameseClu
- **方法类型**：Contrastive / pseudo-positive mining
- **采用部分**：采用高置信正样本挖掘和无/弱负样本对比；避免普通 InfoNCE 把同类或连续状态当负样本。
- **与 scMAE 的最佳结合形式**：正样本来自 teacher consensus neighbors、同 pseudo-core、augmentation views；负样本只用远距离高置信不同簇。
- **建议损失/训练方式**：L = L_scMAE + λ_ins InfoNCE_filtered + λ_cluster swapped/cluster consistency。
- **额外注意问题**：必须 false-negative filtering；稀有细胞不能被大簇拉走。
- **重复/备注**：未提供

### 69. Deep single-cell RNA-seq data clustering with graph prototypical contrastive learning

- **索引推荐**：rank=68，level=可用，score=27
- **来源/年份**：2023 / Bioinformatics / ICML CompBio
- **PDF**：`068_可用_Deep_single-cell_RNA-seq_data_clustering_with_graph_prototypical_contrastive_learning.pdf`
- **代码**：https://github.com/Junseok0207/scGPCL
- **方法类型**：Graph / neighbor reliability
- **采用部分**：采用图方法中的建图、边置信、edge mask、attention、high-order neighbor；不要直接堆多层 GNN。
- **与 scMAE 的最佳结合形式**：图模块作为 NeighborMix controller：输出 edge_confidence、edge_persistence、boundary_score，控制是否 mix。
- **建议损失/训练方式**：L = L_scMAE + λ_edge BCE(masked edge reconstruction) + λ_graph neighbor consistency + λ_dropedge robustness。
- **额外注意问题**：核心风险是过平滑和错边放大；必须 shallow、residual、DropEdge、EMA graph refresh。
- **重复/备注**：未提供

### 70. IGCLAPS: an interpretable graph contrastive learning method with adaptive positive sampling for scRNA-seq data analysis

- **索引推荐**：rank=69，level=可用，score=26
- **来源/年份**：2025 / Bioinformatics
- **PDF**：`069_可用_IGCLAPS_an_interpretable_graph_contrastive_learning_method_with_adaptive_positive_sampling_for_scRNA-seq_data_analysis.pdf`
- **代码**：https://github.com/ZhengWeihuaYNU/IGCLAPS
- **方法类型**：Graph / neighbor reliability
- **采用部分**：采用图方法中的建图、边置信、edge mask、attention、high-order neighbor；不要直接堆多层 GNN。
- **与 scMAE 的最佳结合形式**：图模块作为 NeighborMix controller：输出 edge_confidence、edge_persistence、boundary_score，控制是否 mix。
- **建议损失/训练方式**：L = L_scMAE + λ_edge BCE(masked edge reconstruction) + λ_graph neighbor consistency + λ_dropedge robustness。
- **额外注意问题**：核心风险是过平滑和错边放大；必须 shallow、residual、DropEdge、EMA graph refresh。
- **重复/备注**：未提供

### 71. Deep clustering of single-cell RNA-seq using adversarial graph contrastive learning

- **索引推荐**：rank=70，level=可用，score=25
- **来源/年份**：2025 / Briefings in Bioinformatics
- **PDF**：`070_可用_Deep_clustering_of_single-cell_RNA-seq_using_adversarial_graph_contrastive_learning.pdf`
- **代码**：https://github.com/levinhcntt/scAGCL
- **方法类型**：Graph / neighbor reliability
- **采用部分**：采用图方法中的建图、边置信、edge mask、attention、high-order neighbor；不要直接堆多层 GNN。
- **与 scMAE 的最佳结合形式**：图模块作为 NeighborMix controller：输出 edge_confidence、edge_persistence、boundary_score，控制是否 mix。
- **建议损失/训练方式**：L = L_scMAE + λ_edge BCE(masked edge reconstruction) + λ_graph neighbor consistency + λ_dropedge robustness。
- **额外注意问题**：核心风险是过平滑和错边放大；必须 shallow、residual、DropEdge、EMA graph refresh。
- **重复/备注**：未提供

### 72. scCAN: single-cell clustering using autoencoder and network fusion

- **索引推荐**：rank=71，level=可用，score=24
- **来源/年份**：2022 / Scientific Reports
- **PDF**：`071_可用_scCAN_single-cell_clustering_using_autoencoder_and_network_fusion.pdf`
- **代码**：https://github.com/cran/scCAN
- **方法类型**：Background / cautious reference
- **采用部分**：仅作为背景或消融参考；不要直接作为主模块。
- **与 scMAE 的最佳结合形式**：先判断它补 scMAE 哪个缺口，再决定是否落到 mask、graph、loss、semantic target 或 boundary controller。
- **建议损失/训练方式**：默认使用 L_scMAE，新增 loss 必须有明确机制诊断。
- **额外注意问题**：没有明确结构适配前不进入 formal benchmark。
- **重复/备注**：未提供

### 73. scASDC: Attention Enhanced Structural Deep Clustering for Single-cell RNA-seq Data

- **索引推荐**：rank=72，level=可用，score=23
- **来源/年份**：2024 / arXiv
- **PDF**：`072_可用_scASDC_Attention_Enhanced_Structural_Deep_Clustering_for_Single-cell_RNA-seq_Data.pdf`
- **代码**：https://github.com/wenwenmin/scASDC
- **方法类型**：Graph / neighbor reliability
- **采用部分**：采用图方法中的建图、边置信、edge mask、attention、high-order neighbor；不要直接堆多层 GNN。
- **与 scMAE 的最佳结合形式**：图模块作为 NeighborMix controller：输出 edge_confidence、edge_persistence、boundary_score，控制是否 mix。
- **建议损失/训练方式**：L = L_scMAE + λ_edge BCE(masked edge reconstruction) + λ_graph neighbor consistency + λ_dropedge robustness。
- **额外注意问题**：核心风险是过平滑和错边放大；必须 shallow、residual、DropEdge、EMA graph refresh。
- **重复/备注**：未提供

### 74. Attention-based deep clustering method for scRNA-seq cell type identification

- **索引推荐**：rank=73，level=可用，score=22
- **来源/年份**：2023 / PLOS Computational Biology
- **PDF**：`073_可用_Attention-based_deep_clustering_method_for_scRNA-seq_cell_type_identification.pdf`
- **代码**：https://github.com/ttgump/AttentionAE-sc
- **方法类型**：Graph / neighbor reliability
- **采用部分**：采用图方法中的建图、边置信、edge mask、attention、high-order neighbor；不要直接堆多层 GNN。
- **与 scMAE 的最佳结合形式**：图模块作为 NeighborMix controller：输出 edge_confidence、edge_persistence、boundary_score，控制是否 mix。
- **建议损失/训练方式**：L = L_scMAE + λ_edge BCE(masked edge reconstruction) + λ_graph neighbor consistency + λ_dropedge robustness。
- **额外注意问题**：核心风险是过平滑和错边放大；必须 shallow、residual、DropEdge、EMA graph refresh。
- **重复/备注**：未提供

### 75. Attention-based deep clustering method for scRNA-seq cell type identification

- **索引推荐**：rank=73，level=可用，score=22
- **来源/年份**：PLoS Computational Biology 2023
- **PDF**：`073_可用_Attention-based_deep_clustering_method_for_scRNA-seq_cell_type_identification.pdf`
- **代码**：https://github.com/ttgump/AttentionAE-sc
- **方法类型**：Graph / neighbor reliability
- **采用部分**：采用图方法中的建图、边置信、edge mask、attention、high-order neighbor；不要直接堆多层 GNN。
- **与 scMAE 的最佳结合形式**：图模块作为 NeighborMix controller：输出 edge_confidence、edge_persistence、boundary_score，控制是否 mix。
- **建议损失/训练方式**：L = L_scMAE + λ_edge BCE(masked edge reconstruction) + λ_graph neighbor consistency + λ_dropedge robustness。
- **额外注意问题**：核心风险是过平滑和错边放大；必须 shallow、residual、DropEdge、EMA graph refresh。
- **重复/备注**：duplicate of rank 73

### 76. Boosting scRNA-seq data clustering by cluster-aware feature weighting

- **索引推荐**：rank=74，level=可用，score=21
- **来源/年份**：2021 / BMC Bioinformatics
- **PDF**：`074_可用_Boosting_scRNA-seq_data_clustering_by_cluster-aware_feature_weighting.pdf`
- **代码**：https://github.com/LiRuiyi-raptor/CaFew_Project
- **方法类型**：Feature selection / gene gate
- **采用部分**：采用 feature weighting / group lasso / sparse gate 保护 marker 与可解释性。
- **与 scMAE 的最佳结合形式**：在 scMAE 输入端加入非负 gene gate；gate 由 dispersion、zero-rate、cluster-aware importance 初始化。
- **建议损失/训练方式**：L = L_scMAE + λ_gate L1/group-lasso + λ_marker marker-preservation。
- **额外注意问题**：不要过度稀疏，否则稀有细胞 marker 会被删掉；需报告选基因 enrichment。
- **重复/备注**：未提供

### 77. scGTN: Deep Siamese Graph Transformer Network for Single-cell RNA Sequencing Clustering

- **索引推荐**：rank=75，level=可用，score=20
- **来源/年份**：2026 / arXiv
- **PDF**：`075_可用_scGTN_Deep_Siamese_Graph_Transformer_Network_for_Single-cell_RNA_Sequencing_Clustering.pdf`
- **代码**：未提供
- **方法类型**：Graph / neighbor reliability
- **采用部分**：采用图方法中的建图、边置信、edge mask、attention、high-order neighbor；不要直接堆多层 GNN。
- **与 scMAE 的最佳结合形式**：图模块作为 NeighborMix controller：输出 edge_confidence、edge_persistence、boundary_score，控制是否 mix。
- **建议损失/训练方式**：L = L_scMAE + λ_edge BCE(masked edge reconstruction) + λ_graph neighbor consistency + λ_dropedge robustness。
- **额外注意问题**：核心风险是过平滑和错边放大；必须 shallow、residual、DropEdge、EMA graph refresh。
- **重复/备注**：未提供

### 78. scGraphformer: unveiling cellular heterogeneity and interactions in scRNA-seq data using a scalable graph transformer network

- **索引推荐**：rank=76，level=可用，score=19
- **来源/年份**：2024 / Communications Biology
- **PDF**：`076_可用_scGraphformer_unveiling_cellular_heterogeneity_and_interactions_in_scRNA-seq_data_using_a_scalable_graph_transformer_network.pdf`
- **代码**：https://github.com/xyfan22/scGraphformer
- **方法类型**：Graph / neighbor reliability
- **采用部分**：采用图方法中的建图、边置信、edge mask、attention、high-order neighbor；不要直接堆多层 GNN。
- **与 scMAE 的最佳结合形式**：图模块作为 NeighborMix controller：输出 edge_confidence、edge_persistence、boundary_score，控制是否 mix。
- **建议损失/训练方式**：L = L_scMAE + λ_edge BCE(masked edge reconstruction) + λ_graph neighbor consistency + λ_dropedge robustness。
- **额外注意问题**：核心风险是过平滑和错边放大；必须 shallow、residual、DropEdge、EMA graph refresh。
- **重复/备注**：未提供

### 79. scNET: learning context-specific gene and cell embeddings by integrating single-cell gene expression data with protein-protein interactions

- **索引推荐**：rank=77，level=可用，score=18
- **来源/年份**：2025 / Nature Methods
- **PDF**：`077_可用_scNET_learning_context-specific_gene_and_cell_embeddings_by_integrating_single-cell_gene_expression_data_with_protein-protein_interactions.pdf`
- **代码**：https://github.com/madilabcode/scNET
- **方法类型**：Graph / neighbor reliability
- **采用部分**：采用图方法中的建图、边置信、edge mask、attention、high-order neighbor；不要直接堆多层 GNN。
- **与 scMAE 的最佳结合形式**：图模块作为 NeighborMix controller：输出 edge_confidence、edge_persistence、boundary_score，控制是否 mix。
- **建议损失/训练方式**：L = L_scMAE + λ_edge BCE(masked edge reconstruction) + λ_graph neighbor consistency + λ_dropedge robustness。
- **额外注意问题**：核心风险是过平滑和错边放大；必须 shallow、residual、DropEdge、EMA graph refresh。
- **重复/备注**：未提供

### 80. scLong: a billion-parameter foundation model for capturing long-range gene context in single-cell transcriptomics

- **索引推荐**：rank=78，level=可用，score=17
- **来源/年份**：2026 / Nature Communications
- **PDF**：`078_可用_scLong_a_billion-parameter_foundation_model_for_capturing_long-range_gene_context_in_single-cell_transcriptomics.pdf`
- **代码**：https://github.com/BaiDing1234/scLong
- **方法类型**：Knowledge / interpretability only
- **采用部分**：采用 ontology/LLM/gene text 作为弱知识、解释和后验分析，不作为无监督聚类主监督。
- **与 scMAE 的最佳结合形式**：可把 cell ontology/marker/pathway 用于高置信簇命名、错误边解释、rare-cell validation；训练默认不用标签语义。
- **建议损失/训练方式**：训练损失不应包含 supervised label CE；可选 L_knowledge 仅用于 gene-module attention bias。
- **额外注意问题**：使用外部模型或指令数据会破坏公平性；普通计算机一区论文应保持无监督主体。
- **重复/备注**：未提供

### 81. scHNTL: single-cell RNA-seq data clustering augmented by high-order neighbors and triplet loss

- **索引推荐**：rank=79，level=可用，score=16
- **来源/年份**：2025 / Bioinformatics
- **PDF**：`079_可用_scHNTL_single-cell_RNA-seq_data_clustering_augmented_by_high-order_neighbors_and_triplet_loss.pdf`
- **代码**：https://github.com/SWJTU-ML/scHNTL-code
- **方法类型**：Graph / neighbor reliability
- **采用部分**：采用图方法中的建图、边置信、edge mask、attention、high-order neighbor；不要直接堆多层 GNN。
- **与 scMAE 的最佳结合形式**：图模块作为 NeighborMix controller：输出 edge_confidence、edge_persistence、boundary_score，控制是否 mix。
- **建议损失/训练方式**：L = L_scMAE + λ_edge BCE(masked edge reconstruction) + λ_graph neighbor consistency + λ_dropedge robustness。
- **额外注意问题**：核心风险是过平滑和错边放大；必须 shallow、residual、DropEdge、EMA graph refresh。
- **重复/备注**：未提供

### 82. M-band wavelet-based multi-view clustering of cells

- **索引推荐**：rank=80，level=可用，score=15
- **来源/年份**：2025 / PLOS Computational Biology
- **PDF**：`080_可用_M-band_wavelet-based_multi-view_clustering_of_cells.pdf`
- **代码**：未提供
- **方法类型**：Signal / multi-view decomposition
- **采用部分**：采用 wavelet 多尺度分解构造表达视图，特别适合区分全局趋势与局部 marker 信号。
- **与 scMAE 的最佳结合形式**：把 low-frequency view 作为稳定背景，高-frequency view 用于 rare/marker detection；二者分别进入 scMAE multi-view heads。
- **建议损失/训练方式**：L = L_expr + λ_low reconstruction + λ_high marker preservation + λ_fusion consistency。
- **额外注意问题**：分解参数需固定或数据驱动，避免调参搜索过多。
- **重复/备注**：未提供

### 83. OKR-CELL: Open-world Language Knowledge-Aided Robust Single-Cell Foundation Model

- **索引推荐**：rank=81，level=谨慎/部分可用，score=14
- **来源/年份**：2026 / arXiv / bioRxiv withdrawn version noted
- **PDF**：`081_谨慎_部分可用_OKR-CELL_Open-world_Language_Knowledge-Aided_Robust_Single-Cell_Foundation_Model.pdf`
- **代码**：未提供
- **方法类型**：Knowledge / interpretability only
- **采用部分**：采用 ontology/LLM/gene text 作为弱知识、解释和后验分析，不作为无监督聚类主监督。
- **与 scMAE 的最佳结合形式**：可把 cell ontology/marker/pathway 用于高置信簇命名、错误边解释、rare-cell validation；训练默认不用标签语义。
- **建议损失/训练方式**：训练损失不应包含 supervised label CE；可选 L_knowledge 仅用于 gene-module attention bias。
- **额外注意问题**：使用外部模型或指令数据会破坏公平性；普通计算机一区论文应保持无监督主体。
- **重复/备注**：withdrawn version exists; use with caution

### 84. DiT: Scalable Diffusion Models with Transformers

- **索引推荐**：rank=82，level=谨慎/部分可用，score=13
- **来源/年份**：2023 / ICCV
- **PDF**：`082_谨慎_部分可用_DiT_Scalable_Diffusion_Models_with_Transformers.pdf`
- **代码**：https://github.com/facebookresearch/DiT
- **方法类型**：Generative / denoising auxiliary
- **采用部分**：采用生成模型的噪声调度、denoising target、rare-cell augmentation 思想；不建议训练重型生成主干。
- **与 scMAE 的最佳结合形式**：在 scMAE latent 上做轻量 denoising/consistency，或生成难 corruption；生成样本只用于辅助，不进入正式真实样本评估。
- **建议损失/训练方式**：L = L_scMAE + λ_denoise ||z_clean-z_denoised|| + λ_score/flow latent regularization。
- **额外注意问题**：生成模型容易过重且难验证真实性；必须用 marker/pathway/DEG 保真作为附加评估。
- **重复/备注**：general image diffusion backbone; not scRNA-specific

### 85. scGNN is a novel graph neural network framework for single-cell RNA-Seq analyses

- **索引推荐**：rank=83，level=未明确，score=12
- **来源/年份**：2021 / Nature Communications
- **PDF**：`083_未明确_scGNN_is_a_novel_graph_neural_network_framework_for_single-cell_RNA-Seq_analyses.pdf`
- **代码**：https://github.com/juexinwang/scGNN
- **方法类型**：Graph / neighbor reliability
- **采用部分**：采用图方法中的建图、边置信、edge mask、attention、high-order neighbor；不要直接堆多层 GNN。
- **与 scMAE 的最佳结合形式**：图模块作为 NeighborMix controller：输出 edge_confidence、edge_persistence、boundary_score，控制是否 mix。
- **建议损失/训练方式**：L = L_scMAE + λ_edge BCE(masked edge reconstruction) + λ_graph neighbor consistency + λ_dropedge robustness。
- **额外注意问题**：核心风险是过平滑和错边放大；必须 shallow、residual、DropEdge、EMA graph refresh。
- **重复/备注**：未提供

### 86. scGNN is a novel graph neural network framework for single-cell RNA-seq analyses

- **索引推荐**：rank=83，level=未明确，score=12
- **来源/年份**：Nature Communications 2021
- **PDF**：`083_未明确_scGNN_is_a_novel_graph_neural_network_framework_for_single-cell_RNA-seq_analyses.pdf`
- **代码**：https://github.com/juexinwang/scGNN
- **方法类型**：Graph / neighbor reliability
- **采用部分**：采用图方法中的建图、边置信、edge mask、attention、high-order neighbor；不要直接堆多层 GNN。
- **与 scMAE 的最佳结合形式**：图模块作为 NeighborMix controller：输出 edge_confidence、edge_persistence、boundary_score，控制是否 mix。
- **建议损失/训练方式**：L = L_scMAE + λ_edge BCE(masked edge reconstruction) + λ_graph neighbor consistency + λ_dropedge robustness。
- **额外注意问题**：核心风险是过平滑和错边放大；必须 shallow、residual、DropEdge、EMA graph refresh。
- **重复/备注**：duplicate of rank 83

### 87. Unsupervised Deep Embedding for Clustering Analysis

- **索引推荐**：rank=84，level=未明确，score=11
- **来源/年份**：ICML 2016
- **PDF**：`084_未明确_Unsupervised_Deep_Embedding_for_Clustering_Analysis.pdf`
- **代码**：https://github.com/piiswrong/dec
- **方法类型**：Boundary / clustering head
- **采用部分**：采用软聚类、模糊 membership、rough lower/upper approximation 或社区发现作为诊断和后期约束。
- **与 scMAE 的最佳结合形式**：先用 scMAE/teacher graph 得到核心细胞；只对核心细胞加强聚类，边界细胞允许高熵。
- **建议损失/训练方式**：L = L_scMAE + λ_core KL(P||Q) + λ_fuzzy entropy/balance + λ_sep center separation。
- **额外注意问题**：不要早期强行 DEC；known-k 假设需在论文中明确。
- **重复/备注**：未提供

### 88. BERT: Pre-training of Deep Bidirectional Transformers for Language Understanding

- **索引推荐**：rank=85，level=未明确，score=10
- **来源/年份**：NAACL 2019
- **PDF**：`085_未明确_BERT_Pre-training_of_Deep_Bidirectional_Transformers_for_Language_Understanding.pdf`
- **代码**：https://github.com/google-research/bert
- **方法类型**：Masked modeling / discriminator target
- **采用部分**：采用 token 化、replaced-expression detection、自引导 mask、条件多阶段 mask，但必须映射到基因表达结构。
- **与 scMAE 的最佳结合形式**：保留 scMAE corruption；加 gene-specific rank tokens、mask difficulty schedule、ELECTRA/VIME 判别头。
- **建议损失/训练方式**：L = L_scMAE_rec + λ_mask CE(replaced) + λ_tok CE(rank_token) + λ_curriculum adaptive mask。
- **额外注意问题**：避免图像/文本 token 机械迁移；gene token 应按 gene-specific distribution 或 module 生成。
- **重复/备注**：未提供

### 89. A Simple Framework for Contrastive Learning of Visual Representations

- **索引推荐**：rank=86，level=未明确，score=9
- **来源/年份**：ICML 2020
- **PDF**：`086_未明确_A_Simple_Framework_for_Contrastive_Learning_of_Visual_Representations.pdf`
- **代码**：https://github.com/google-research/simclr
- **方法类型**：Contrastive / pseudo-positive mining
- **采用部分**：采用高置信正样本挖掘和无/弱负样本对比；避免普通 InfoNCE 把同类或连续状态当负样本。
- **与 scMAE 的最佳结合形式**：正样本来自 teacher consensus neighbors、同 pseudo-core、augmentation views；负样本只用远距离高置信不同簇。
- **建议损失/训练方式**：L = L_scMAE + λ_ins InfoNCE_filtered + λ_cluster swapped/cluster consistency。
- **额外注意问题**：必须 false-negative filtering；稀有细胞不能被大簇拉走。
- **重复/备注**：未提供

### 90. Bootstrap Your Own Latent: A New Approach to Self-Supervised Learning

- **索引推荐**：rank=87，level=未明确，score=8
- **来源/年份**：NeurIPS 2020
- **PDF**：`087_未明确_Bootstrap_Your_Own_Latent_A_New_Approach_to_Self-Supervised_Learning.pdf`
- **代码**：https://github.com/deepmind/deepmind-research/tree/master/byol
- **方法类型**：Teacher-student self-distillation
- **采用部分**：采用 EMA teacher、center/sharpen 或 BYOL predictor，目标是稳定表示和图，而非外部 foundation weight。
- **与 scMAE 的最佳结合形式**：teacher 输入弱扰动/clean view，student 输入强 mask view；teacher embedding 同时用于 consensus graph。
- **建议损失/训练方式**：L = L_scMAE + λ_distill ||p(z_s)-sg(z_t)|| 或 CE(centered teacher probs)；可加 collapse diagnostics。
- **额外注意问题**：不能只加一个 teacher MSE；需要 warmup、EMA decay schedule、variance/covariance 防塌陷监控。
- **重复/备注**：未提供

### 91. Emerging Properties in Self-Supervised Vision Transformers

- **索引推荐**：rank=88，level=未明确，score=7
- **来源/年份**：ICCV 2021
- **PDF**：`088_未明确_Emerging_Properties_in_Self-Supervised_Vision_Transformers.pdf`
- **代码**：https://github.com/facebookresearch/dino
- **方法类型**：Teacher-student self-distillation
- **采用部分**：采用 EMA teacher、center/sharpen 或 BYOL predictor，目标是稳定表示和图，而非外部 foundation weight。
- **与 scMAE 的最佳结合形式**：teacher 输入弱扰动/clean view，student 输入强 mask view；teacher embedding 同时用于 consensus graph。
- **建议损失/训练方式**：L = L_scMAE + λ_distill ||p(z_s)-sg(z_t)|| 或 CE(centered teacher probs)；可加 collapse diagnostics。
- **额外注意问题**：不能只加一个 teacher MSE；需要 warmup、EMA decay schedule、variance/covariance 防塌陷监控。
- **重复/备注**：未提供

### 92. scziDesk: deep soft K-means clustering with self-training

- **索引推荐**：rank=89，level=未明确，score=6
- **来源/年份**：NAR Genomics and Bioinformatics 2020
- **PDF**：`089_未明确_scziDesk_deep_soft_K-means_clustering_with_self-training.pdf`
- **代码**：https://github.com/xuebaliang/scziDesk
- **方法类型**：Boundary / clustering head
- **采用部分**：采用软聚类、模糊 membership、rough lower/upper approximation 或社区发现作为诊断和后期约束。
- **与 scMAE 的最佳结合形式**：先用 scMAE/teacher graph 得到核心细胞；只对核心细胞加强聚类，边界细胞允许高熵。
- **建议损失/训练方式**：L = L_scMAE + λ_core KL(P||Q) + λ_fuzzy entropy/balance + λ_sep center separation。
- **额外注意问题**：不要早期强行 DEC；known-k 假设需在论文中明确。
- **重复/备注**：未提供

### 93. ELECTRA: Pre-training Text Encoders as Discriminators Rather Than Generators

- **索引推荐**：rank=90，level=未明确，score=5
- **来源/年份**：ICLR 2020
- **PDF**：`090_未明确_ELECTRA_Pre-training_Text_Encoders_as_Discriminators_Rather_Than_Generators.pdf`
- **代码**：https://github.com/google-research/electra
- **方法类型**：Masked modeling / discriminator target
- **采用部分**：采用 token 化、replaced-expression detection、自引导 mask、条件多阶段 mask，但必须映射到基因表达结构。
- **与 scMAE 的最佳结合形式**：保留 scMAE corruption；加 gene-specific rank tokens、mask difficulty schedule、ELECTRA/VIME 判别头。
- **建议损失/训练方式**：L = L_scMAE_rec + λ_mask CE(replaced) + λ_tok CE(rank_token) + λ_curriculum adaptive mask。
- **额外注意问题**：避免图像/文本 token 机械迁移；gene token 应按 gene-specific distribution 或 module 生成。
- **重复/备注**：未提供

### 94. GraphMAE: Self-Supervised Masked Graph Autoencoders

- **索引推荐**：rank=91，level=未明确，score=4
- **来源/年份**：KDD 2022
- **PDF**：`091_未明确_GraphMAE_Self-Supervised_Masked_Graph_Autoencoders.pdf`
- **代码**：https://github.com/THUDM/GraphMAE
- **方法类型**：Graph / neighbor reliability
- **采用部分**：采用图方法中的建图、边置信、edge mask、attention、high-order neighbor；不要直接堆多层 GNN。
- **与 scMAE 的最佳结合形式**：图模块作为 NeighborMix controller：输出 edge_confidence、edge_persistence、boundary_score，控制是否 mix。
- **建议损失/训练方式**：L = L_scMAE + λ_edge BCE(masked edge reconstruction) + λ_graph neighbor consistency + λ_dropedge robustness。
- **额外注意问题**：核心风险是过平滑和错边放大；必须 shallow、residual、DropEdge、EMA graph refresh。
- **重复/备注**：未提供

### 95. scDSC: Deep structural clustering jointly through autoencoder and GNN

- **索引推荐**：rank=92，level=未明确，score=3
- **来源/年份**：Briefings in Bioinformatics 2022
- **PDF**：`092_未明确_scDSC_Deep_structural_clustering_jointly_through_autoencoder_and_GNN.pdf`
- **代码**：https://github.com/DHUDBlab/scDSC
- **方法类型**：Graph / neighbor reliability
- **采用部分**：采用图方法中的建图、边置信、edge mask、attention、high-order neighbor；不要直接堆多层 GNN。
- **与 scMAE 的最佳结合形式**：图模块作为 NeighborMix controller：输出 edge_confidence、edge_persistence、boundary_score，控制是否 mix。
- **建议损失/训练方式**：L = L_scMAE + λ_edge BCE(masked edge reconstruction) + λ_graph neighbor consistency + λ_dropedge robustness。
- **额外注意问题**：核心风险是过平滑和错边放大；必须 shallow、residual、DropEdge、EMA graph refresh。
- **重复/备注**：未提供

### 96. GraphMAE2: A Decoding-Enhanced Masked Self-Supervised Graph Learner

- **索引推荐**：rank=93，level=未明确，score=2
- **来源/年份**：WWW 2023
- **PDF**：`093_未明确_GraphMAE2_A_Decoding-Enhanced_Masked_Self-Supervised_Graph_Learner.pdf`
- **代码**：https://github.com/THUDM/GraphMAE2
- **方法类型**：Graph / neighbor reliability
- **采用部分**：采用图方法中的建图、边置信、edge mask、attention、high-order neighbor；不要直接堆多层 GNN。
- **与 scMAE 的最佳结合形式**：图模块作为 NeighborMix controller：输出 edge_confidence、edge_persistence、boundary_score，控制是否 mix。
- **建议损失/训练方式**：L = L_scMAE + λ_edge BCE(masked edge reconstruction) + λ_graph neighbor consistency + λ_dropedge robustness。
- **额外注意问题**：核心风险是过平滑和错边放大；必须 shallow、residual、DropEdge、EMA graph refresh。
- **重复/备注**：未提供

### 97. Graph Contrastive Learning with Augmentations

- **索引推荐**：rank=94，level=未明确，score=1
- **来源/年份**：NeurIPS 2020
- **PDF**：`094_未明确_Graph_Contrastive_Learning_with_Augmentations.pdf`
- **代码**：https://github.com/Shen-Lab/GraphCL
- **方法类型**：Contrastive / pseudo-positive mining
- **采用部分**：采用高置信正样本挖掘和无/弱负样本对比；避免普通 InfoNCE 把同类或连续状态当负样本。
- **与 scMAE 的最佳结合形式**：正样本来自 teacher consensus neighbors、同 pseudo-core、augmentation views；负样本只用远距离高置信不同簇。
- **建议损失/训练方式**：L = L_scMAE + λ_ins InfoNCE_filtered + λ_cluster swapped/cluster consistency。
- **额外注意问题**：必须 false-negative filtering；稀有细胞不能被大簇拉走。
- **重复/备注**：未提供

### 98. scCDCG: Deep Cut-informed Graph Embedding for scRNA-seq clustering

- **索引推荐**：rank=95，level=未明确，score=0
- **来源/年份**：Recent preprint/relevant scRNA method
- **PDF**：`095_未明确_scCDCG_Deep_Cut-informed_Graph_Embedding_for_scRNA-seq_clustering.pdf`
- **代码**：https://github.com/XPgogogo/scCDCG
- **方法类型**：Graph / neighbor reliability
- **采用部分**：采用图方法中的建图、边置信、edge mask、attention、high-order neighbor；不要直接堆多层 GNN。
- **与 scMAE 的最佳结合形式**：图模块作为 NeighborMix controller：输出 edge_confidence、edge_persistence、boundary_score，控制是否 mix。
- **建议损失/训练方式**：L = L_scMAE + λ_edge BCE(masked edge reconstruction) + λ_graph neighbor consistency + λ_dropedge robustness。
- **额外注意问题**：核心风险是过平滑和错边放大；必须 shallow、residual、DropEdge、EMA graph refresh。
- **重复/备注**：未提供

### 99. From Louvain to Leiden: guaranteeing well-connected communities

- **索引推荐**：rank=96，level=未明确，score=-1
- **来源/年份**：Scientific Reports 2019
- **PDF**：`096_未明确_From_Louvain_to_Leiden_guaranteeing_well-connected_communities.pdf`
- **代码**：https://github.com/vtraag/leidenalg
- **方法类型**：Boundary / clustering head
- **采用部分**：采用软聚类、模糊 membership、rough lower/upper approximation 或社区发现作为诊断和后期约束。
- **与 scMAE 的最佳结合形式**：先用 scMAE/teacher graph 得到核心细胞；只对核心细胞加强聚类，边界细胞允许高熵。
- **建议损失/训练方式**：L = L_scMAE + λ_core KL(P||Q) + λ_fuzzy entropy/balance + λ_sep center separation。
- **额外注意问题**：不要早期强行 DEC；known-k 假设需在论文中明确。
- **重复/备注**：未提供

### 100. Comprehensive integration of single-cell data

- **索引推荐**：rank=97，level=未明确，score=-2
- **来源/年份**：Cell 2019 / preprint
- **PDF**：`097_未明确_Comprehensive_integration_of_single-cell_data.pdf`
- **代码**：https://github.com/seandavi/awesome-single-cell
- **方法类型**：Batch/multi-view integration
- **采用部分**：采用 batch alignment、multi-view fusion、shared/private latent 的思想作为多数据集鲁棒性模块。
- **与 scMAE 的最佳结合形式**：scMAE 主干输出 shared latent；batch/platform view 只用 adversarial/OT alignment 或 view consistency。
- **建议损失/训练方式**：L = L_scMAE + λ_align MMD/OT/adversarial + λ_view multi-view consistency。
- **额外注意问题**：不要把当前单数据集聚类问题复杂化；作为跨批次扩展和附加实验。
- **重复/备注**：未提供

### 101. Graph Attention Networks

- **索引推荐**：rank=98，level=未明确，score=-3
- **来源/年份**：ICLR 2018
- **PDF**：`098_未明确_Graph_Attention_Networks.pdf`
- **代码**：https://github.com/PetarV-/GAT
- **方法类型**：Graph / neighbor reliability
- **采用部分**：采用图方法中的建图、边置信、edge mask、attention、high-order neighbor；不要直接堆多层 GNN。
- **与 scMAE 的最佳结合形式**：图模块作为 NeighborMix controller：输出 edge_confidence、edge_persistence、boundary_score，控制是否 mix。
- **建议损失/训练方式**：L = L_scMAE + λ_edge BCE(masked edge reconstruction) + λ_graph neighbor consistency + λ_dropedge robustness。
- **额外注意问题**：核心风险是过平滑和错边放大；必须 shallow、residual、DropEdge、EMA graph refresh。
- **重复/备注**：未提供

### 102. DropEdge: Towards Deep Graph Convolutional Networks on Node Classification

- **索引推荐**：rank=99，level=未明确，score=-4
- **来源/年份**：ICLR 2020
- **PDF**：`099_未明确_DropEdge_Towards_Deep_Graph_Convolutional_Networks_on_Node_Classification.pdf`
- **代码**：https://github.com/DropEdge/DropEdge
- **方法类型**：Graph / neighbor reliability
- **采用部分**：采用图方法中的建图、边置信、edge mask、attention、high-order neighbor；不要直接堆多层 GNN。
- **与 scMAE 的最佳结合形式**：图模块作为 NeighborMix controller：输出 edge_confidence、edge_persistence、boundary_score，控制是否 mix。
- **建议损失/训练方式**：L = L_scMAE + λ_edge BCE(masked edge reconstruction) + λ_graph neighbor consistency + λ_dropedge robustness。
- **额外注意问题**：核心风险是过平滑和错边放大；必须 shallow、residual、DropEdge、EMA graph refresh。
- **重复/备注**：未提供

### 103. Denoising Diffusion Probabilistic Models

- **索引推荐**：rank=100，level=未明确，score=-5
- **来源/年份**：NeurIPS 2020
- **PDF**：`100_未明确_Denoising_Diffusion_Probabilistic_Models.pdf`
- **代码**：https://github.com/hojonathanho/diffusion
- **方法类型**：Generative / denoising auxiliary
- **采用部分**：采用生成模型的噪声调度、denoising target、rare-cell augmentation 思想；不建议训练重型生成主干。
- **与 scMAE 的最佳结合形式**：在 scMAE latent 上做轻量 denoising/consistency，或生成难 corruption；生成样本只用于辅助，不进入正式真实样本评估。
- **建议损失/训练方式**：L = L_scMAE + λ_denoise ||z_clean-z_denoised|| + λ_score/flow latent regularization。
- **额外注意问题**：生成模型容易过重且难验证真实性；必须用 marker/pathway/DEG 保真作为附加评估。
- **重复/备注**：未提供

### 104. Neural Discrete Representation Learning

- **索引推荐**：rank=101，level=未明确，score=-6
- **来源/年份**：NeurIPS 2017
- **PDF**：`101_未明确_Neural_Discrete_Representation_Learning.pdf`
- **代码**：https://github.com/deepmind/sonnet/blob/v1/sonnet/python/modules/nets/vqvae.py
- **方法类型**：Masked modeling / discriminator target
- **采用部分**：采用 token 化、replaced-expression detection、自引导 mask、条件多阶段 mask，但必须映射到基因表达结构。
- **与 scMAE 的最佳结合形式**：保留 scMAE corruption；加 gene-specific rank tokens、mask difficulty schedule、ELECTRA/VIME 判别头。
- **建议损失/训练方式**：L = L_scMAE_rec + λ_mask CE(replaced) + λ_tok CE(rank_token) + λ_curriculum adaptive mask。
- **额外注意问题**：避免图像/文本 token 机械迁移；gene token 应按 gene-specific distribution 或 module 生成。
- **重复/备注**：未提供

### 105. mixup: Beyond Empirical Risk Minimization

- **索引推荐**：rank=102，level=未明确，score=-7
- **来源/年份**：ICLR 2018
- **PDF**：`102_未明确_mixup_Beyond_Empirical_Risk_Minimization.pdf`
- **代码**：https://github.com/facebookresearch/mixup-cifar10
- **方法类型**：Mixing principle
- **采用部分**：采用 mixup 的局部线性假设，但必须约束在可靠细胞邻域内，这正是 NeighborMix 的理论来源之一。
- **与 scMAE 的最佳结合形式**：只在 teacher-consensus mutual neighbors 间做 expression/latent mix；禁止全局随机 mix。
- **建议损失/训练方式**：L = L_scMAE(mixed input) + λ_cons mixed latent consistency。
- **额外注意问题**：全局 mixup 会产生不真实细胞；必须 boundary/rare-cell veto。
- **重复/备注**：未提供

### 106. AutoAugment: Learning Augmentation Strategies from Data

- **索引推荐**：rank=103，level=未明确，score=-8
- **来源/年份**：CVPR 2019
- **PDF**：`103_未明确_AutoAugment_Learning_Augmentation_Strategies_from_Data.pdf`
- **代码**：https://github.com/tensorflow/models/tree/master/research/autoaugment
- **方法类型**：Policy search for mask/mix
- **采用部分**：采用自适应策略搜索思想，搜索 mask type、mask ratio、mix strength、neighbor refresh interval。
- **与 scMAE 的最佳结合形式**：作为 policy controller 包在 scMAE/NeighborMix 外层；搜索空间必须小：random/dropout/module mask × mix_alpha × graph_refresh。
- **建议损失/训练方式**：L_policy 以验证 ARI proxy、loss plateau、graph purity 或 teacher-student agreement 为 reward；训练主损失仍是 scMAE。
- **额外注意问题**：不要让 RL/AutoAugment 成为主贡献；否则工程复杂且不可解释。
- **重复/备注**：未提供

### 107. VIME: Extending the Success of Self- and Semi-supervised Learning to Tabular Domain

- **索引推荐**：rank=104，level=未明确，score=-9
- **来源/年份**：NeurIPS 2020
- **PDF**：`104_未明确_VIME_Extending_the_Success_of_Self-_and_Semi-supervised_Learning_to_Tabular_Domain.pdf`
- **代码**：https://github.com/jsyoon0823/VIME
- **方法类型**：Masked modeling / discriminator target
- **采用部分**：采用 token 化、replaced-expression detection、自引导 mask、条件多阶段 mask，但必须映射到基因表达结构。
- **与 scMAE 的最佳结合形式**：保留 scMAE corruption；加 gene-specific rank tokens、mask difficulty schedule、ELECTRA/VIME 判别头。
- **建议损失/训练方式**：L = L_scMAE_rec + λ_mask CE(replaced) + λ_tok CE(rank_token) + λ_curriculum adaptive mask。
- **额外注意问题**：避免图像/文本 token 机械迁移；gene token 应按 gene-specific distribution 或 module 生成。
- **重复/备注**：未提供

### 108. Deep generative modeling for single-cell transcriptomics

- **索引推荐**：rank=105，level=未明确，score=-10
- **来源/年份**：Nature Methods 2018
- **PDF**：`105_未明确_Deep_generative_modeling_for_single-cell_transcriptomics.pdf`
- **代码**：https://github.com/scverse/scvi-tools
- **方法类型**：Generative / denoising auxiliary
- **采用部分**：采用生成模型的噪声调度、denoising target、rare-cell augmentation 思想；不建议训练重型生成主干。
- **与 scMAE 的最佳结合形式**：在 scMAE latent 上做轻量 denoising/consistency，或生成难 corruption；生成样本只用于辅助，不进入正式真实样本评估。
- **建议损失/训练方式**：L = L_scMAE + λ_denoise ||z_clean-z_denoised|| + λ_score/flow latent regularization。
- **额外注意问题**：生成模型容易过重且难验证真实性；必须用 marker/pathway/DEG 保真作为附加评估。
- **重复/备注**：未提供

### 109. Score-Based Generative Modeling through Stochastic Differential Equations

- **索引推荐**：rank=106，level=未明确，score=-11
- **来源/年份**：ICLR 2021
- **PDF**：`106_未明确_Score-Based_Generative_Modeling_through_Stochastic_Differential_Equations.pdf`
- **代码**：https://github.com/yang-song/score_sde
- **方法类型**：Generative / denoising auxiliary
- **采用部分**：采用生成模型的噪声调度、denoising target、rare-cell augmentation 思想；不建议训练重型生成主干。
- **与 scMAE 的最佳结合形式**：在 scMAE latent 上做轻量 denoising/consistency，或生成难 corruption；生成样本只用于辅助，不进入正式真实样本评估。
- **建议损失/训练方式**：L = L_scMAE + λ_denoise ||z_clean-z_denoised|| + λ_score/flow latent regularization。
- **额外注意问题**：生成模型容易过重且难验证真实性；必须用 marker/pathway/DEG 保真作为附加评估。
- **重复/备注**：未提供

### 110. TabNet: Attentive Interpretable Tabular Learning

- **索引推荐**：rank=107，level=未明确，score=-12
- **来源/年份**：AAAI 2021
- **PDF**：`107_未明确_TabNet_Attentive_Interpretable_Tabular_Learning.pdf`
- **代码**：https://github.com/google-research/google-research/tree/master/tabnet
- **方法类型**：Tabular architecture / row attention
- **采用部分**：采用表格学习中的 feature tokenizer、row attention、feature importance；不替换生物 mask 主线。
- **与 scMAE 的最佳结合形式**：在 scMAE encoder 前加入 gene-wise tokenizer/gate，或用 row-neighbor attention 替代静态 context。
- **建议损失/训练方式**：L = L_scMAE + λ_gate sparsity + λ_ret row-neighbor consistency。
- **额外注意问题**：表格方法通常面向监督分类；无监督聚类下需避免标签先验。
- **重复/备注**：未提供

### 111. SimMIM: A Simple Framework for Masked Image Modeling

- **索引推荐**：rank=108，level=基线/背景，score=-13
- **来源/年份**：2022 / CVPR
- **PDF**：`108_基线_背景_SimMIM_A_Simple_Framework_for_Masked_Image_Modeling.pdf`
- **代码**：https://github.com/microsoft/SimMIM
- **方法类型**：Masked modeling / discriminator target
- **采用部分**：采用 token 化、replaced-expression detection、自引导 mask、条件多阶段 mask，但必须映射到基因表达结构。
- **与 scMAE 的最佳结合形式**：保留 scMAE corruption；加 gene-specific rank tokens、mask difficulty schedule、ELECTRA/VIME 判别头。
- **建议损失/训练方式**：L = L_scMAE_rec + λ_mask CE(replaced) + λ_tok CE(rank_token) + λ_curriculum adaptive mask。
- **额外注意问题**：避免图像/文本 token 机械迁移；gene token 应按 gene-specific distribution 或 module 生成。
- **重复/备注**：background baseline for masked image modeling

### 112. Revisiting Deep Learning Models for Tabular Data / FT-Transformer

- **索引推荐**：rank=109，level=基线/背景，score=-14
- **来源/年份**：2021 / NeurIPS
- **PDF**：`109_基线_背景_Revisiting_Deep_Learning_Models_for_Tabular_Data_FT-Transformer.pdf`
- **代码**：https://github.com/yandex-research/rtdl-revisiting-models
- **方法类型**：Tabular architecture / row attention
- **采用部分**：采用表格学习中的 feature tokenizer、row attention、feature importance；不替换生物 mask 主线。
- **与 scMAE 的最佳结合形式**：在 scMAE encoder 前加入 gene-wise tokenizer/gate，或用 row-neighbor attention 替代静态 context。
- **建议损失/训练方式**：L = L_scMAE + λ_gate sparsity + λ_ret row-neighbor consistency。
- **额外注意问题**：表格方法通常面向监督分类；无监督聚类下需避免标签先验。
- **重复/备注**：background for tabular transformer design

### 113. AdaMAE: Adaptive Masking for Efficient Spatiotemporal Learning with Masked Autoencoders

- **索引推荐**：rank=110，level=基线/背景，score=-15
- **来源/年份**：2023 / CVPR
- **PDF**：`110_基线_背景_AdaMAE_Adaptive_Masking_for_Efficient_Spatiotemporal_Learning_with_Masked_Autoencoders.pdf`
- **代码**：https://github.com/wgcban/adamae
- **方法类型**：Masked modeling / discriminator target
- **采用部分**：采用 token 化、replaced-expression detection、自引导 mask、条件多阶段 mask，但必须映射到基因表达结构。
- **与 scMAE 的最佳结合形式**：保留 scMAE corruption；加 gene-specific rank tokens、mask difficulty schedule、ELECTRA/VIME 判别头。
- **建议损失/训练方式**：L = L_scMAE_rec + λ_mask CE(replaced) + λ_tok CE(rank_token) + λ_curriculum adaptive mask。
- **额外注意问题**：避免图像/文本 token 机械迁移；gene token 应按 gene-specific distribution 或 module 生成。
- **重复/备注**：background for adaptive masking

### 114. scMAE: a masked autoencoder for single-cell RNA-seq clustering

- **索引推荐**：rank=111，level=基线/背景，score=-16
- **来源/年份**：Bioinformatics 2024
- **PDF**：`111_基线_背景_scMAE_a_masked_autoencoder_for_single-cell_RNA-seq_clustering.pdf`
- **代码**：未提供
- **方法类型**：Masked modeling / discriminator target
- **采用部分**：采用 token 化、replaced-expression detection、自引导 mask、条件多阶段 mask，但必须映射到基因表达结构。
- **与 scMAE 的最佳结合形式**：保留 scMAE corruption；加 gene-specific rank tokens、mask difficulty schedule、ELECTRA/VIME 判别头。
- **建议损失/训练方式**：L = L_scMAE_rec + λ_mask CE(replaced) + λ_tok CE(rank_token) + λ_curriculum adaptive mask。
- **额外注意问题**：避免图像/文本 token 机械迁移；gene token 应按 gene-specific distribution 或 module 生成。
- **重复/备注**：original baseline

### 115. SC3: consensus clustering of single-cell RNA-seq data

- **索引推荐**：rank=112，level=基线/背景，score=-17
- **来源/年份**：Nature Methods 2017
- **PDF**：`112_基线_背景_SC3_consensus_clustering_of_single-cell_RNA-seq_data.pdf`
- **代码**：https://github.com/hemberg-lab/SC3
- **方法类型**：Boundary / clustering head
- **采用部分**：采用软聚类、模糊 membership、rough lower/upper approximation 或社区发现作为诊断和后期约束。
- **与 scMAE 的最佳结合形式**：先用 scMAE/teacher graph 得到核心细胞；只对核心细胞加强聚类，边界细胞允许高熵。
- **建议损失/训练方式**：L = L_scMAE + λ_core KL(P||Q) + λ_fuzzy entropy/balance + λ_sep center separation。
- **额外注意问题**：不要早期强行 DEC；known-k 假设需在论文中明确。
- **重复/备注**：classical scRNA clustering baseline

### 116. Fast unfolding of communities in large networks

- **索引推荐**：rank=113，level=基线/背景，score=-18
- **来源/年份**：J. Stat. Mech. 2008
- **PDF**：`113_基线_背景_Fast_unfolding_of_communities_in_large_networks.pdf`
- **代码**：未提供
- **方法类型**：Boundary / clustering head
- **采用部分**：采用软聚类、模糊 membership、rough lower/upper approximation 或社区发现作为诊断和后期约束。
- **与 scMAE 的最佳结合形式**：先用 scMAE/teacher graph 得到核心细胞；只对核心细胞加强聚类，边界细胞允许高熵。
- **建议损失/训练方式**：L = L_scMAE + λ_core KL(P||Q) + λ_fuzzy entropy/balance + λ_sep center separation。
- **额外注意问题**：不要早期强行 DEC；known-k 假设需在论文中明确。
- **重复/备注**：Louvain baseline background

### 117. scMAE: a masked autoencoder for single-cell RNA-seq clustering

- **索引推荐**：rank=未提供，level=未提供，score=未提供
- **来源/年份**：2024 / Bioinformatics
- **PDF**：`未提供`
- **代码**：未提供
- **方法类型**：Masked modeling / discriminator target
- **采用部分**：采用 token 化、replaced-expression detection、自引导 mask、条件多阶段 mask，但必须映射到基因表达结构。
- **与 scMAE 的最佳结合形式**：保留 scMAE corruption；加 gene-specific rank tokens、mask difficulty schedule、ELECTRA/VIME 判别头。
- **建议损失/训练方式**：L = L_scMAE_rec + λ_mask CE(replaced) + λ_tok CE(rank_token) + λ_curriculum adaptive mask。
- **额外注意问题**：避免图像/文本 token 机械迁移；gene token 应按 gene-specific distribution 或 module 生成。
- **重复/备注**：from previous package

### 118. Variational Masked AutoEncoder topic / related work cluster

- **索引推荐**：rank=未提供，level=未提供，score=未提供
- **来源/年份**：2026 / topic/various
- **PDF**：`未提供`
- **代码**：未提供
- **方法类型**：Distribution / VAE uncertainty
- **采用部分**：采用 NB/ZINB/variational uncertainty 作为重构和不确定性模块，不建议替换整个 scMAE backbone。
- **与 scMAE 的最佳结合形式**：保留 scMAE encoder；decoder 分出 μ/θ/π 或 latent μ/logσ，边界细胞可根据 posterior variance 降低 mix 强度。
- **建议损失/训练方式**：L = L_scMAE + λ_NB NLL_NB/ZINB(raw_counts|μ,θ,π) + β KL(q(z|x)||p(z))。
- **额外注意问题**：必须使用 raw counts + size factor 分支；在 scaled expression 上做 ZINB 是结构性错误。
- **重复/备注**：from previous package

### 119. scVAE: variational auto-encoders for single-cell gene expression data

- **索引推荐**：rank=未提供，level=未提供，score=未提供
- **来源/年份**：2020 / Bioinformatics
- **PDF**：`未提供`
- **代码**：https://github.com/scvae/scvae
- **方法类型**：Distribution / VAE uncertainty
- **采用部分**：采用 NB/ZINB/variational uncertainty 作为重构和不确定性模块，不建议替换整个 scMAE backbone。
- **与 scMAE 的最佳结合形式**：保留 scMAE encoder；decoder 分出 μ/θ/π 或 latent μ/logσ，边界细胞可根据 posterior variance 降低 mix 强度。
- **建议损失/训练方式**：L = L_scMAE + λ_NB NLL_NB/ZINB(raw_counts|μ,θ,π) + β KL(q(z|x)||p(z))。
- **额外注意问题**：必须使用 raw counts + size factor 分支；在 scaled expression 上做 ZINB 是结构性错误。
- **重复/备注**：from previous package

### 120. scGAC: a graph attentional architecture for clustering single-cell RNA-seq data

- **索引推荐**：rank=未提供，level=未提供，score=未提供
- **来源/年份**：2022 / Bioinformatics
- **PDF**：`未提供`
- **代码**：https://github.com/Joye9285/scGAC
- **方法类型**：Graph / neighbor reliability
- **采用部分**：采用图方法中的建图、边置信、edge mask、attention、high-order neighbor；不要直接堆多层 GNN。
- **与 scMAE 的最佳结合形式**：图模块作为 NeighborMix controller：输出 edge_confidence、edge_persistence、boundary_score，控制是否 mix。
- **建议损失/训练方式**：L = L_scMAE + λ_edge BCE(masked edge reconstruction) + λ_graph neighbor consistency + λ_dropedge robustness。
- **额外注意问题**：核心风险是过平滑和错边放大；必须 shallow、residual、DropEdge、EMA graph refresh。
- **重复/备注**：from previous package

### 121. scGAAC: A graph attention autoencoder for clustering single-cell RNA-sequencing data

- **索引推荐**：rank=未提供，level=未提供，score=未提供
- **来源/年份**：2024 / Methods / Elsevier
- **PDF**：`未提供`
- **代码**：https://github.com/labiip/scGAAC
- **方法类型**：Graph / neighbor reliability
- **采用部分**：采用图方法中的建图、边置信、edge mask、attention、high-order neighbor；不要直接堆多层 GNN。
- **与 scMAE 的最佳结合形式**：图模块作为 NeighborMix controller：输出 edge_confidence、edge_persistence、boundary_score，控制是否 mix。
- **建议损失/训练方式**：L = L_scMAE + λ_edge BCE(masked edge reconstruction) + λ_graph neighbor consistency + λ_dropedge robustness。
- **额外注意问题**：核心风险是过平滑和错边放大；必须 shallow、residual、DropEdge、EMA graph refresh。
- **重复/备注**：from previous package

### 122. scSCC: A swapped contrastive learning-based clustering method for single-cell gene expression data

- **索引推荐**：rank=未提供，level=未提供，score=未提供
- **来源/年份**：2025 / Quantitative Biology
- **PDF**：`未提供`
- **代码**：未提供
- **方法类型**：Contrastive / pseudo-positive mining
- **采用部分**：采用高置信正样本挖掘和无/弱负样本对比；避免普通 InfoNCE 把同类或连续状态当负样本。
- **与 scMAE 的最佳结合形式**：正样本来自 teacher consensus neighbors、同 pseudo-core、augmentation views；负样本只用远距离高置信不同簇。
- **建议损失/训练方式**：L = L_scMAE + λ_ins InfoNCE_filtered + λ_cluster swapped/cluster consistency。
- **额外注意问题**：必须 false-negative filtering；稀有细胞不能被大簇拉走。
- **重复/备注**：from previous package

### 123. Decoupled GNNs based on multi-view contrastive learning for scRNA-seq data clustering

- **索引推荐**：rank=未提供，level=未提供，score=未提供
- **来源/年份**：2025 / Briefings in Bioinformatics
- **PDF**：`未提供`
- **代码**：未提供
- **方法类型**：Contrastive / pseudo-positive mining
- **采用部分**：采用高置信正样本挖掘和无/弱负样本对比；避免普通 InfoNCE 把同类或连续状态当负样本。
- **与 scMAE 的最佳结合形式**：正样本来自 teacher consensus neighbors、同 pseudo-core、augmentation views；负样本只用远距离高置信不同簇。
- **建议损失/训练方式**：L = L_scMAE + λ_ins InfoNCE_filtered + λ_cluster swapped/cluster consistency。
- **额外注意问题**：必须 false-negative filtering；稀有细胞不能被大簇拉走。
- **重复/备注**：from previous package

### 124. ScCCL: Single-Cell Data Clustering Based on Self-Supervised Contrastive Learning

- **索引推荐**：rank=未提供，level=未提供，score=未提供
- **来源/年份**：2023 / IEEE/ACM TCBB
- **PDF**：`未提供`
- **代码**：https://github.com/LuckyxiaoLin/ScCCL
- **方法类型**：Contrastive / pseudo-positive mining
- **采用部分**：采用高置信正样本挖掘和无/弱负样本对比；避免普通 InfoNCE 把同类或连续状态当负样本。
- **与 scMAE 的最佳结合形式**：正样本来自 teacher consensus neighbors、同 pseudo-core、augmentation views；负样本只用远距离高置信不同簇。
- **建议损失/训练方式**：L = L_scMAE + λ_ins InfoNCE_filtered + λ_cluster swapped/cluster consistency。
- **额外注意问题**：必须 false-negative filtering；稀有细胞不能被大簇拉走。
- **重复/备注**：from previous package

### 125. Multi-view clustering for single-cell RNA-seq data based on graph fusion

- **索引推荐**：rank=未提供，level=未提供，score=未提供
- **来源/年份**：2025 / Briefings in Bioinformatics
- **PDF**：`未提供`
- **代码**：https://github.com/WJ319/scMCGF
- **方法类型**：Graph / neighbor reliability
- **采用部分**：采用图方法中的建图、边置信、edge mask、attention、high-order neighbor；不要直接堆多层 GNN。
- **与 scMAE 的最佳结合形式**：图模块作为 NeighborMix controller：输出 edge_confidence、edge_persistence、boundary_score，控制是否 mix。
- **建议损失/训练方式**：L = L_scMAE + λ_edge BCE(masked edge reconstruction) + λ_graph neighbor consistency + λ_dropedge robustness。
- **额外注意问题**：核心风险是过平滑和错边放大；必须 shallow、residual、DropEdge、EMA graph refresh。
- **重复/备注**：from previous package

### 126. scMVAF: a multi-view adaptive fusion clustering approach for single-cell RNA-sequencing data

- **索引推荐**：rank=未提供，level=未提供，score=未提供
- **来源/年份**：2026 / Briefings in Bioinformatics
- **PDF**：`未提供`
- **代码**：https://github.com/LQXLE/scMVAF
- **方法类型**：Batch/multi-view integration
- **采用部分**：采用 batch alignment、multi-view fusion、shared/private latent 的思想作为多数据集鲁棒性模块。
- **与 scMAE 的最佳结合形式**：scMAE 主干输出 shared latent；batch/platform view 只用 adversarial/OT alignment 或 view consistency。
- **建议损失/训练方式**：L = L_scMAE + λ_align MMD/OT/adversarial + λ_view multi-view consistency。
- **额外注意问题**：不要把当前单数据集聚类问题复杂化；作为跨批次扩展和附加实验。
- **重复/备注**：from previous package

### 127. scVIC: deep generative modeling of heterogeneity for scRNA-seq data

- **索引推荐**：rank=未提供，level=未提供，score=未提供
- **来源/年份**：2024 / Bioinformatics Advances
- **PDF**：`未提供`
- **代码**：https://github.com/HiBearME/scVIC
- **方法类型**：Distribution / VAE uncertainty
- **采用部分**：采用 NB/ZINB/variational uncertainty 作为重构和不确定性模块，不建议替换整个 scMAE backbone。
- **与 scMAE 的最佳结合形式**：保留 scMAE encoder；decoder 分出 μ/θ/π 或 latent μ/logσ，边界细胞可根据 posterior variance 降低 mix 强度。
- **建议损失/训练方式**：L = L_scMAE + λ_NB NLL_NB/ZINB(raw_counts|μ,θ,π) + β KL(q(z|x)||p(z))。
- **额外注意问题**：必须使用 raw counts + size factor 分支；在 scaled expression 上做 ZINB 是结构性错误。
- **重复/备注**：from previous package

### 128. scGANSL: Graph Attention Network with Subspace Learning for scRNA-seq Data Clustering

- **索引推荐**：rank=未提供，level=未提供，score=未提供
- **来源/年份**：2025 / Journal of Chemical Information and Modeling / ACS
- **PDF**：`未提供`
- **代码**：https://github.com/szq0816/scMMN
- **方法类型**：Graph / neighbor reliability
- **采用部分**：采用图方法中的建图、边置信、edge mask、attention、high-order neighbor；不要直接堆多层 GNN。
- **与 scMAE 的最佳结合形式**：图模块作为 NeighborMix controller：输出 edge_confidence、edge_persistence、boundary_score，控制是否 mix。
- **建议损失/训练方式**：L = L_scMAE + λ_edge BCE(masked edge reconstruction) + λ_graph neighbor consistency + λ_dropedge robustness。
- **额外注意问题**：核心风险是过平滑和错边放大；必须 shallow、residual、DropEdge、EMA graph refresh。
- **重复/备注**：from previous package

### 129. scDFC: A deep fusion clustering method for single-cell RNA-seq data

- **索引推荐**：rank=未提供，level=未提供，score=未提供
- **来源/年份**：2023 / Briefings in Bioinformatics
- **PDF**：`未提供`
- **代码**：https://github.com/DayuHuu/scDFC
- **方法类型**：Background / cautious reference
- **采用部分**：仅作为背景或消融参考；不要直接作为主模块。
- **与 scMAE 的最佳结合形式**：先判断它补 scMAE 哪个缺口，再决定是否落到 mask、graph、loss、semantic target 或 boundary controller。
- **建议损失/训练方式**：默认使用 L_scMAE，新增 loss 必须有明确机制诊断。
- **额外注意问题**：没有明确结构适配前不进入 formal benchmark。
- **重复/备注**：from previous package

### 130. Integrating feature selection with unsupervised deep embedding for clustering single-cell RNA-seq data

- **索引推荐**：rank=未提供，level=未提供，score=未提供
- **来源/年份**：2026 / Briefings in Bioinformatics
- **PDF**：`未提供`
- **代码**：未提供
- **方法类型**：Feature selection / gene gate
- **采用部分**：采用 feature weighting / group lasso / sparse gate 保护 marker 与可解释性。
- **与 scMAE 的最佳结合形式**：在 scMAE 输入端加入非负 gene gate；gate 由 dispersion、zero-rate、cluster-aware importance 初始化。
- **建议损失/训练方式**：L = L_scMAE + λ_gate L1/group-lasso + λ_marker marker-preservation。
- **额外注意问题**：不要过度稀疏，否则稀有细胞 marker 会被删掉；需报告选基因 enrichment。
- **重复/备注**：from previous package




## 5. 本文档的使用方式

后续每实现一个方法，都必须在实现前回答：

1. 它补 scMAE 的哪个缺口？
2. 它是否需要 raw/log/scaled/rank 哪一路输入？
3. 它的接入点是 mask、target、graph、mix、decoder 还是 clustering？
4. 它是否与 NeighborMix 冲突？
5. 它是否保护 boundary / rare cells？
6. 它是否能输出机制诊断？
7. 它是否达到 independent-full 和 formal-verified 的证据等级？

若无法回答，不进入正式实现。
