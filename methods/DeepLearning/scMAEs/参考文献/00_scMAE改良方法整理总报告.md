# scMAE改良方法收集整理总报告

本报告由三个文库的PDF、Markdown报告、逐篇TXT建议和manifest CSV整理生成。PDF文件已按推荐程度从高到低复制到 `01_PDF论文_按推荐程度排序`。

## 整理摘要

- 去重后PDF数量：113
- 有建议但缺少PDF的条目：14
- 去重规则：优先按SHA-256完全去重；对明确同题名论文再做题名级合并。
- 说明来源优先级：原Markdown报告；文库2逐篇TXT；缺少逐篇Markdown说明时使用 `manifest_30_papers.csv` 和 `manifest_extra50_papers.csv`。

## 原始综合报告保留

### 文库1：scMAE_recent50_report.md

# 近5年额外50篇论文：scMAE改良发散调研报告

本包是在上一版30篇资料之外新增的50条记录，优先选择近5年计算机顶会/顶刊、arXiv高相关预印本、单细胞基础模型、掩码建模、图自监督、生成模型、时序/信号、模糊粗糙集、表格学习等方向。

- PDF已下载：47 篇，位于 `papers/`。
- PDF未能下载/需谨慎：3 篇，摘要和方法记录在 `notes/not_downloaded_summaries.md`。
- 完整清单：`manifest_extra50_papers.csv`。

## 对scMAE最值得优先尝试的5条主线

### A. Mask策略升级：从随机shuffle mask到生物/信号感知mask
可参考：ScDiVa, BEiT, SimMIM, MaskFeat, iBOT, I-JEPA, MaskGIT。

1. 是否有效：有效。scMAE的核心瓶颈是mask任务是否足够“难且有生物意义”。随机扰动表达值容易学到技术噪声；预测gene-set activity、rank、TF activity或离散token可更接近细胞状态。
2. 如何改良：设计三类mask：随机基因mask、通路块mask、dropout-like absorbing mask；重构目标从MSE扩展为表达值+基因模块特征+离散token。
3. 可能问题：通路块mask过强导致模型只学通路先验，忽略未知基因关系。
4. 规避：混合mask策略，并使用JOAO式策略搜索或验证集自适应选择；保留随机mask作为基线。
5. 是否可用：可用，推荐作为第一主创新。

### B. 分布建模：从MSE到NB/ZINB/鲁棒损失
可参考：scVGAE, scAGC, Improved Consistency Models。

1. 是否有效：高概率有效。scRNA count具备稀疏、过离散、dropout特征，MSE对零值和高表达离群值敏感。
2. 如何改良：decoder输出均值、离散度和dropout概率；或对log-normalized表达使用Pseudo-Huber重构。
3. 可能问题：ZINB并非对所有平台都成立，模型可能过拟合零值。
4. 规避：NB/ZINB/Huber三种loss可切换，按数据平台自动选择；加入gene-wise zero-rate calibration。
5. 是否可用：可用，推荐作为稳定增益模块。

### C. 细胞图结构：从静态KNN到动态/可学习/掩码图
可参考：scAGC, BGRL, MaskGAE, Graph Barlow Twins, GraphGPS, Graphormer, TabR。

1. 是否有效：有效。scRNA聚类本来常依赖KNN图和Leiden/Louvain，scMAE只做表达mask会弱化细胞间关系。
2. 如何改良：构建cell graph encoder，加入edge mask reconstruction、nearest-neighbor retrieval context、Graphormer distance bias、BGRL/Barlow负样本自由一致性。
3. 可能问题：KNN图有噪声，过平滑会让相邻但不同类型细胞混合。
4. 规避：动态边置信度、edge dropout、Gumbel-Softmax稀疏采样、只用高置信mutual nearest neighbors。
5. 是否可用：可用，推荐作为第二主创新。

### D. 不确定性/模糊粗糙集：处理过渡细胞和边界细胞
可参考：Deep Adaptive Fuzzy Clustering, Fuzzy Rough Sets Based on Fuzzy Quantification。

1. 是否有效：对单细胞尤其有效。许多细胞处于连续谱系或状态转换中，硬聚类会把边界细胞错误固定到某类。
2. 如何改良：clustering head输出membership；用lower approximation核心细胞做强监督，upper boundary细胞做弱监督或一致性正则。
3. 可能问题：硬标签指标可能下降，训练更复杂。
4. 规避：最终输出仍可argmax硬标签；额外报告soft silhouette、entropy、trajectory consistency。
5. 是否可用：可用，适合增强论文解释性。

### E. 生成与扰动：latent diffusion/flow matching/consistency作为增强而不是主线
可参考：Latent Diffusion, DiT, Flow Matching, Consistency Models, ScDiVa。

1. 是否有效：对扰动预测、稀有细胞增强、batch alignment可能有效。
2. 如何改良：先用scMAE编码得到latent，再用flow/diffusion做条件生成或OT alignment。
3. 可能问题：生成模型过重，且合成细胞真实性难验证。
4. 规避：不把生成样本直接当真实样本；用于latent regularization、counterfactual analysis或补充实验。
5. 是否可用：可用，但不建议作为第一版核心，除非算力充足。

## 推荐的一区论文方案草案

**BioGraph-Tabular Masked Autoencoder for Single-Cell Clustering**

核心模块：
1. 生物信号感知mask：random gene + pathway block + dropout-like absorbing mask。
2. 多目标重构：表达值重构 + pathway/TF activity feature target + NB/ZINB或Pseudo-Huber。
3. 图一致性：cell graph retrieval context + MaskGAE edge reconstruction + BGRL/Barlow一致性。
4. 模糊边界聚类：核心细胞强约束，边界细胞弱约束。
5. 可解释性：mask importance、gene-set enrichment、marker recovery、rare-cell preservation。

实验建议：
- Baseline：scMAE、DEC/scDeepCluster/scDSC/scGNN/AttentionAE-sc、Leiden/Seurat、scVGAE/scAGC等。
- 指标：ARI/NMI/ACC、silhouette、batch mixing、rare-cell recall、marker enrichment、运行时间。
- 消融：mask策略、loss、图模块、fuzzy head、retrieval context。
- 数据：Tabula Muris/Baron pancreas/Muraro/PBMC/CITE-seq/多批次数据。

## 文件说明

- `papers/`：已下载PDF。
- `manifest_extra50_papers.csv`：50篇/条记录，逐条包含摘要思想、scMAE改良方式、有效性、问题与规避、可用性。
- `notes/not_downloaded_summaries.md`：未下载PDF的论文详细记录。

### 文库3：scMAE_improvement_report.md

# scMAE 改良方向调研报告（面向普通计算机类一区论文）

本报告围绕 scMAE（masked autoencoder for scRNA-seq clustering）的可发表改良点整理。压缩包内包含 30 篇 PDF，以及 `manifest_30_papers.csv`。其中大部分为顶会/顶刊/高影响力期刊论文，另含你已上传的若干直接相关 scRNA-seq 聚类论文，用于建立基线和对照。

## 0. scMAE 基线判断

scMAE 的核心思想可概括为：对每个基因构造 Bernoulli mask，用同一基因在其他细胞中的 shuffled value 替换被 mask 的表达值；网络由 encoder、mask predictor、decoder 组成；mask predictor 用交叉熵预测哪些位置被替换，decoder 用加权 MSE 重构；最后用 K-means 或 Leiden 对 latent embedding 聚类。

这带来三个可改良的切入点：

1. **mask 机制仍较粗糙**：随机基因级 Bernoulli mask 没有显式考虑 dropout、dispersion、基因模块、marker gene、pathway/GRN。
2. **重构损失偏连续信号视角**：weighted MSE 对 scRNA count 的过离散、library size、batch effect、zero inflation 的统计建模不够。
3. **结构信息参与较弱**：scMAE 本体主要在表格式表达矩阵上做 masked modeling，聚类阶段才用 K-means/Leiden，未在预训练阶段强制 latent space 保持细胞图结构、社区结构或伪时间结构。

## 1. 最值得做成论文主线的方案

建议把论文主线控制在 2~3 个强模块，避免把所有技术堆进去。推荐主线：

**Bio-Graph Discriminative scMAE（暂名）**

核心贡献：

1. **生物学/信号感知 mask**：按基因 dropout 率、dispersion、marker score、gene-gene correlation、pathway/module 设计 adaptive mask，而不是均匀随机 mask。
2. **ELECTRA/VIME 式 replaced-expression discrimination**：不仅重构表达值，还判别每个位置是否被替换；对未 mask 位置也学习判别信号，提高样本效率。
3. **图结构一致性正则**：在 encoder latent space 上构造 KNN/Leiden 图，引入 GraphMAE/GraphCL/GAT/DropEdge 思想，防止 cell-cell topology 被重构任务破坏。

可选增强：NB/ZINB likelihood、unbalanced optimal transport/fuzzy clustering head。

一个可实现的损失函数：

```text
L = λ_rec L_count_or_wMSE + λ_mask CE(M, M_hat)
  + λ_graph L_graphMAE + λ_clu KL(P || Q)
  + λ_cons L_teacher_student + λ_reg L_sparse/pathway
```

建议第一版不要同时加入 diffusion/RL/foundation model，除非已有稳定代码和算力。

---

## 2. 候选方法逐项思考

### 方法 1：生物学与信号感知的自适应 Mask

**是否有效**：可用，优先级最高。scRNA-seq 的有效信号不是均匀分布在所有基因上，高变基因、marker genes、pathway genes、dropout-sensitive genes、housekeeping genes 承载的信息不同。随机 mask 会浪费预算，也可能破坏稀有细胞 marker。

**原理**：从信号处理角度，mask 相当于对表达矩阵进行缺失采样。均匀采样假设所有维度同等重要，但 scRNA-seq 是异方差、高稀疏、过离散信号。更合理的是对基因 j 使用 p_j，而非常数 p。

**如何改良 scMAE**：

- 令 mask 概率为：
  ```text
  p_j = sigmoid(a·dispersion_j + b·dropout_j + c·module_score_j - d·marker_risk_j)
  ```
- 对 pathway/gene module 做 block mask：一次 mask 掉一组相关基因，迫使模型学习跨模块依赖。
- 对 rare-cell marker 设置保护项，避免模型在稀有细胞上过度扰动。

**可能问题**：过度依赖外部 marker/pathway 会降低无监督性；mask 策略太复杂会被审稿人认为调参。

**规避**：只用训练数据内统计量构造默认 mask；外部 pathway 作为可选 ablation；提供固定公式和少量超参；报告多数据集鲁棒性。

**结论**：可用，而且最适合作为第一创新点。

---

### 方法 2：ELECTRA/VIME 式 Replaced-expression Detection

**是否有效**：可用，优先级很高。scMAE 已经预测 mask，但可以更彻底地把它设计成“表达值是否被替换”的判别式预训练任务。

**原理**：生成式重构只在 mask 位置有强监督，而 replaced-token/expression detection 可在所有位置提供密集监督。对于高维稀疏表格数据，VIME 也证明 mask estimation + reconstruction 是有效的 tabular SSL 目标。

**如何改良 scMAE**：

- 将 shuffled replacement 替换为更难的 generator replacement：从同 batch/同近邻/同 pathway 条件分布中采样替换值。
- 增加 discriminator head：预测每个 gene-cell entry 是否 original/replaced。
- 对 mask CE 使用 class-balanced focal loss，避免未替换位置太多导致 trivial classifier。

**可能问题**：替换值太容易识别时模型只学技术噪声；替换值太难时训练不稳定。

**规避**：采用 curriculum：从 random shuffle → neighbor-aware shuffle → generator replacement；控制 replacement 与原值分布匹配；加入 batch-aware replacement 防止模型学批次差异。

**结论**：可用，和 scMAE 原始机制贴合，容易形成清晰贡献。

---

### 方法 3：NB/ZINB 计数分布重构替代 MSE

**是否有效**：可用，但需谨慎。scRNA count 常表现为过离散、library-size 影响和 dropout。MSE 对低表达和高表达基因的误差刻画不一定合理。

**原理**：NB/ZINB likelihood 是生物统计建模，能把表达均值、离散度、dropout probability 分开估计，减少重构 loss 对高表达基因的支配。

**如何改良 scMAE**：

- decoder 输出 μ、θ、π，训练负对数似然：`L_NB` 或 `L_ZINB`。
- 对 UMI 数据默认 NB，对非 UMI 或高 dropout 数据启用 ZINB。
- 使用 library size factor 归一化 μ。

**可能问题**：ZINB 是否必要存在争议；分布建模过强可能导致过拟合或过度 imputation。

**规避**：提供 NB/ZINB/MSE 三组 ablation；让模型自动选择 NB/ZINB 或把 ZINB 作为可选；评估 DEG 保真和 marker gene 保真，而不仅是 ARI/NMI。

**结论**：可用，适合作为第二或第三模块，不建议单独作为唯一创新点。

---

### 方法 4：GraphMAE/GraphMAE2 式细胞图 Masked Autoencoder

**是否有效**：可用，优先级高。细胞并非独立样本，相似细胞在流形或 KNN 图上相邻，图 masked modeling 能让 latent embedding 保持拓扑结构。

**原理**：GraphMAE 在图上 mask node attributes 并重构，学习局部邻域与全局拓扑。scRNA 的 cell-cell KNN 图、gene-gene 共表达图、cell-gene bipartite graph 都适合图 SSL。

**如何改良 scMAE**：

- 在 scMAE encoder 后构建 cell graph，加入图编码器或 graph regularizer。
- 同时重构 masked gene expression 和 masked cell-node attributes。
- 引入 edge reconstruction / neighbor consistency loss。

**可能问题**：KNN 图早期质量差；图神经网络容易过平滑；大数据集 O(n²) 邻接矩阵成本高。

**规避**：用 EMA embedding 周期性更新 KNN；使用 approximate nearest neighbor；采用 DropEdge、shallow GAT/GCN、residual connection；只保留稀疏邻接。

**结论**：可用，适合作为结构创新点。

---

### 方法 5：GAT/Graph Transformer 注意力邻居聚合

**是否有效**：可用，但不宜过重。GAT 能学习不同邻居权重，适合细胞图中“近邻有噪声”的问题。

**原理**：KNN 图的每条边不等价，错误邻居或跨批次邻居会伤害表示。attention 可以动态降低低质量边权。

**如何改良 scMAE**：

- 在 latent cell graph 上用 GAT 聚合邻居，得到 topology-aware latent。
- 注意力权重可加入 batch penalty、marker similarity、pathway similarity。
- 输出 attention heatmap 作为可解释性材料。

**可能问题**：attention 不一定等于生物解释；图 transformer 计算量大；小数据集上容易过拟合。

**规避**：使用轻量 GAT，只做 1~2 层；加 DropEdge 和 attention entropy regularization；解释性结论只作为辅助，不作为强生物发现。

**结论**：可用，但建议作为 graph module 内部实现，而不是主标题。

---

### 方法 6：BYOL/DINO/GraphCL/SimCLR 式一致性与对比学习

**是否有效**：可用。scRNA-seq 中细胞状态连续，负样本容易误伤同类型或相邻状态细胞，因此 BYOL/DINO 这类无负样本方法比纯 SimCLR 更安全。

**原理**：同一细胞的两种合理扰动视图应映射到相近 latent。EMA teacher、centering、sharpening 可防止坍塌。

**如何改良 scMAE**：

- 两个 view：mask+shuffle、dropout simulation、library-size perturbation、gene module dropout。
- online encoder 与 EMA teacher encoder 做 latent consistency。
- 若使用 GraphCL，则对 cell graph 做 edge drop/subgraph sampling。

**可能问题**：增强策略不合理会改变细胞类型；对比学习的负样本会伤害稀有细胞或连续谱系。

**规避**：优先 BYOL/DINO 而不是强负样本 SimCLR；如果用 contrastive，使用 batch-aware/cluster-aware false-negative filtering。

**结论**：可用，适合作为稳定性增强模块。

---

### 方法 7：粗糙集与模糊集的不确定性聚类头

**是否有效**：可用，但建议作为“小而精”的辅助创新。细胞类型边界本来模糊，尤其是发育轨迹、过渡态、双细胞、激活状态。

**原理**：模糊集用 membership 表达一个细胞属于多个簇的程度；粗糙集用 lower/upper approximation 表达确定核心和边界样本。这比硬 K-means/Leiden 更符合生物状态连续性。

**如何改良 scMAE**：

- 在 clustering head 上输出 membership `u_ik`。
- 高置信样本进入 lower approximation，用于 DEC/OT self-training；边界样本进入 upper approximation，只参与弱监督或 consistency loss。
- 加入 entropy/compactness loss：核心样本低熵，边界样本允许高熵。

**可能问题**：指标 ARI/NMI 偏硬标签，模糊输出可能不占优势；参数太多。

**规避**：最终仍输出 hard labels 用于指标；同时报告 boundary-cell biological plausibility，如 marker mixedness、pseudotime continuity；只设置一个阈值或用置信度自适应阈值。

**结论**：可用，能形成区别于常规深度聚类的理论亮点。

---

### 方法 8：Optimal Transport / Sinkhorn 平衡或非平衡聚类

**是否有效**：可用，但必须使用非平衡版本。scRNA-seq cell type 分布高度不均衡，强制均衡会伤害真实 rare cell。

**原理**：OT 可把软分配 Q 映射到目标分布 P，避免所有样本坍缩到少数簇。非平衡 OT 允许簇大小偏离先验。

**如何改良 scMAE**：

- 用 Sinkhorn 得到 target assignment P，替代 DEC 的简单平方归一化。
- 簇大小先验来自 Leiden 初始簇或 Dirichlet smoothing。
- 对 rare cluster 设置 lower-bound，而不是 equal-size constraint。

**可能问题**：错误的簇大小先验会强行扭曲真实结构。

**规避**：unbalanced OT；先验随训练 EMA 更新；对低置信样本降低权重。

**结论**：可用，建议和 fuzzy/rough boundary 结合。

---

### 方法 9：Pathway/GRN/LTMG 生物调控信号正则

**是否有效**：可用。仅从表达矩阵学习可能受 dropout 和 batch effect 干扰，基因调控状态和 pathway 先验能提高生物合理性。

**原理**：细胞类型由基因调控网络和功能通路决定。把 gene regulatory signal 离散化或把 pathway module 作为结构先验，可提升信噪比。

**如何改良 scMAE**：

- 对 gene embedding 加 pathway group regularization。
- 用 LTMG/离散表达状态作为辅助预测任务。
- 加入 gene-gene graph，使模型同时学习 cell-cell 和 gene-gene 结构。

**可能问题**：外部数据库不完整；不同物种/组织通路知识偏差大。

**规避**：默认用数据驱动共表达模块；外部 pathway 只作为可选增强；跨物种实验单独报告。

**结论**：可用，若论文偏“计算机+生信”，这是增强可信度的好点。

---

### 方法 10：Diffusion / Score-based 生成模型用于 rare-cell 增强

**是否有效**：作为主创新风险较高，作为辅助可用。生成模型可以做 denoising、uncertainty estimation、rare-cell augmentation，但实现复杂。

**原理**：扩散模型通过逐步加噪/去噪学习数据分布；在 scRNA 上可用于恢复表达、生成稀有细胞邻域样本或构造更困难的 corruption。

**如何改良 scMAE**：

- 不建议一开始训练完整 cell-level diffusion。
- 更稳妥：在 latent space 训练轻量 diffusion/score denoiser，生成 rare-cell neighborhood augmentations。
- 用生成样本只参与 representation consistency，不直接参与 DE analysis。

**可能问题**：合成细胞可能 hallucinate marker；审稿人会质疑生物真实性；算力成本高。

**规避**：只作为辅助增强；用 density filtering、marker consistency、nearest-neighbor biological validation 过滤；报告真实数据指标而非只展示生成样本。

**结论**：可用但不推荐作为第一篇论文主线。

---

### 方法 11：RL/AutoAugment 搜索 Mask/Corruption 策略

**是否有效**：可用但成本高。对于 scMAE，mask ratio、replacement policy、gene block size、dropout augmentation 都是关键超参，自动搜索可能带来提升。

**原理**：AutoAugment 用强化学习搜索数据增强策略。scMAE 可把 policy 定义为 mask/replacement 操作序列，用验证指标或无监督 proxy reward 优化。

**如何改良 scMAE**：

- 搜索空间：mask ratio、gene module mask、neighbor shuffle、batch-aware shuffle、dropout simulation、mixup strength。
- reward：validation reconstruction + silhouette + cluster stability + rare-cell preservation。

**可能问题**：无监督 reward 与最终 ARI/NMI 不一致；搜索耗时；容易过拟合数据集。

**规避**：使用低成本 Bayesian search/Hyperband 替代完整 RL；policy 在多个训练集上共享；只把搜索作为发现默认策略的工具，不作为推理时模块。

**结论**：可用，但更适合作为实验工具和补充贡献。

---

### 方法 12：TabNet 式序列注意力特征选择

**是否有效**：可用。scRNA 表达矩阵是高维稀疏表格，TabNet 的 instance-wise sparse feature selection 可用于发现每个细胞最关键的基因子集。

**原理**：不同细胞类型依赖不同 marker/pathway。全局固定高变基因选择会忽略 cell-state-specific features。

**如何改良 scMAE**：

- 在 encoder 前加入 sparse feature gate，产生 cell-specific gene mask。
- 让自适应 mask 与 TabNet gate 互补：高重要性基因不总是被 mask，而是以 controlled probability 做 hard reconstruction。
- 输出 gene importance 用于可解释性和 marker recovery。

**可能问题**：门控会把低表达稀有 marker 过滤掉；可解释性可能不稳定。

**规避**：加入稀疏但不硬删除的 soft gate；稳定性选择，多随机种子取 consensus gene importance。

**结论**：可用，适合作为可解释性增强。

---

## 3. 推荐实验设计

### 3.1 数据集

复用 scMAE 原文的 15 个公开 scRNA-seq 数据集，并额外选择：

- 小样本数据集：检验过拟合和稳定性。
- 大样本数据集（>10k cells）：检验图模块复杂度。
- 稀有细胞数据集：报告 rare-cell recall / macro-F1。
- 跨平台数据：检验 batch robustness。

### 3.2 Baselines

必须包括：scMAE、scNAME、scGNN、scDSC、AttentionAE-sc、scVI、SC3、Seurat/Leiden、DEC/scziDesk。若引入图模块，还要对比 GraphMAE-style ablation。

### 3.3 指标

- 聚类：ARI、NMI、ACC、macro-F1、rare-cell recall。
- 内部结构：silhouette、Davies-Bouldin、graph modularity、cluster stability。
- 生物学：marker gene recovery、DEG overlap、GO/KEGG enrichment coherence。
- 工程：runtime、GPU memory、scalability。

### 3.4 Ablation

至少做：

1. scMAE 原始 mask vs adaptive mask。
2. MSE vs NB/ZINB。
3. 无 graph vs GAT/GraphMAE vs GraphMAE+DropEdge。
4. 无 discriminator vs replaced-expression discriminator。
5. K-means/Leiden hard head vs DEC/OT/fuzzy head。
6. 不同 mask ratio 与 block mask size。

## 4. 风险控制

- **创新过多但不聚焦**：主文只保留 2~3 个模块，其余放 ablation 或 future work。
- **生物验证不足**：至少做 marker gene 和 pathway enrichment，不要只堆 ARI/NMI。
- **图模块被质疑依赖错误 KNN**：使用动态图、DropEdge、不同 K 的敏感性分析。
- **稀有细胞损伤**：报告 macro-F1 和 rare-cell recall；用非平衡 OT 和 rare-marker保护 mask。
- **算力过高**：给出 sparse graph、mini-batch、ANN KNN 的实现细节。

## 5. 结论

最可发表、实现难度和创新性均衡的方向是：

> 在 scMAE 中加入生物学/信号感知的 adaptive corruption，配合 ELECTRA/VIME 式 replaced-expression detection，并用轻量图结构一致性正则保持细胞拓扑。

这个组合的优点是：

1. 与 scMAE 原始贡献自然衔接，不是硬拼模块。
2. 从数学上有采样/判别式预训练/图正则支撑。
3. 从生物学上能解释为什么 mask 不应均匀随机。
4. 从计算机论文角度有明确模型结构、损失函数和消融空间。
5. 相比完整 diffusion/foundation model，工程风险较低。

## 按推荐程度排序的PDF与逐篇说明

### 001. ScDiVa: Masked Discrete Diffusion for Joint Modeling of Single-Cell Identity and Expression

- **排序文件**：`01_PDF论文_按推荐程度排序/001_很高_ScDiVa_Masked_Discrete_Diffusion_for_Joint_Modeling_of_Single-Cell_Identity_and_Expression.pdf`
- **推荐等级**：很高（score=92）
- **SHA-256**：`B7F8751A154473B1E21613EA6DC433C5A5A3A33245EC97E476F8E8D9F12E0454`
- **来源记录**：文库1#1 01_ScDiVa_Wang_2026_arxiv.pdf

#### 来源说明：文库1 #1

说明来源：文库1/manifest_extra50_papers.csv
- **标题**：ScDiVa: Masked Discrete Diffusion for Joint Modeling of Single-Cell Identity and Expression
- **年份/会议**：2026 / arXiv / ICML-style preprint
- **方向/领域**：single-cell foundation model / masked discrete diffusion
- **核心思想**：用吸收态[MASK]的离散扩散模拟scRNA dropout，双向denoiser同时恢复基因身份与表达值；强调无序集合、稀疏性和深度不变训练。
- **如何改良 scMAE**：把scMAE的随机shuffle mask升级为连续时间masking diffusion：不同t模拟不同测序深度；目标改为gene-id分类+value回归。
- **有效性依据**：高：噪声过程与dropout/缺失检测高度一致。
- **潜在问题**：训练成本高；基因tokenization和表达值离散/连续混合可能复杂。
- **规避方案**：先做小模型：只对高变基因/marker gene做离散扩散，并保留MAE重构作为辅助损失。
- **最终可用性**：可用，适合作为主创新候选。
- **来源链接**：https://arxiv.org/abs/2602.03477

### 002. MaskFeat: Masked Feature Prediction for Self-Supervised Visual Pre-Training

- **排序文件**：`01_PDF论文_按推荐程度排序/002_很高_MaskFeat_Masked_Feature_Prediction_for_Self-Supervised_Visual_Pre-Training.pdf`
- **推荐等级**：很高（score=92）
- **SHA-256**：`BD2D2BDB8BDB31C52BED056D65E36850E2E6D99691402C9FC2ACC9E092B02E43`
- **来源记录**：文库1#19 19_MaskFeat_Wei_2022_CVPR_arxiv.pdf

#### 来源说明：文库1 #19

说明来源：文库1/manifest_extra50_papers.csv
- **标题**：MaskFeat: Masked Feature Prediction for Self-Supervised Visual Pre-Training
- **年份/会议**：2022 / CVPR
- **方向/领域**：feature target prediction / signal descriptors
- **核心思想**：预测HOG等人工特征而不是原始像素，强调局部对比归一化。
- **如何改良 scMAE**：把scMAE目标从原始表达重构扩展为“信号特征重构”：gene-set scores、TF activity、pathway activity、rank statistics。
- **有效性依据**：高：减少对噪声count的过拟合。
- **潜在问题**：人工特征可能丢失细节。
- **规避方案**：多目标：原始表达+通路/排名/局部相关特征。
- **最终可用性**：可用，推荐。
- **来源链接**：https://arxiv.org/abs/2112.09133

### 003. Improved Techniques for Training Consistency Models

- **排序文件**：`01_PDF论文_按推荐程度排序/003_很高_Improved_Techniques_for_Training_Consistency_Models.pdf`
- **推荐等级**：很高（score=92）
- **SHA-256**：`2C77B5AE689BE478B2AFE9F3A8149B4E4755A357C67BE0A0C96FA06AC0D45B31`
- **来源记录**：文库1#32 32_Improved_Consistency_Models_Song_2024_ICLR.pdf

#### 来源说明：文库1 #32

说明来源：文库1/manifest_extra50_papers.csv
- **标题**：Improved Techniques for Training Consistency Models
- **年份/会议**：2024 / ICLR
- **方向/领域**：fast generative modeling / robust loss
- **核心思想**：改进一致性模型训练，使用Pseudo-Huber等鲁棒统计思想。
- **如何改良 scMAE**：把scMAE重构loss改为Pseudo-Huber/robust loss，降低outlier/dropout影响。
- **有效性依据**：高：实现简单且有效。
- **潜在问题**：过强鲁棒会忽略真实稀有高表达。
- **规避方案**：对marker/高置信基因降低鲁棒截断，对噪声基因增强。
- **最终可用性**：可用，推荐。
- **来源链接**：https://arxiv.org/abs/2310.14189

### 004. JOAO: Automated Data Augmentations for Graph Contrastive Learning

- **排序文件**：`01_PDF论文_按推荐程度排序/004_很高_JOAO_Automated_Data_Augmentations_for_Graph_Contrastive_Learning.pdf`
- **推荐等级**：很高（score=92）
- **SHA-256**：`7616489E730B95FACACFFBB3FDDD12E42FFB6DCAE4EEADBE579A7CC7FF5EE765`
- **来源记录**：文库1#44 44_JOAO_You_2021_ICML_arxiv.pdf

#### 来源说明：文库1 #44

说明来源：文库1/manifest_extra50_papers.csv
- **标题**：JOAO: Automated Data Augmentations for Graph Contrastive Learning
- **年份/会议**：2021 / ICML
- **方向/领域**：augmentation policy search
- **核心思想**：自动选择图增强组合，提高图对比学习。
- **如何改良 scMAE**：自动搜索scMAE mask policy：基因随机/通路/block/dropout/graph edge mask。
- **有效性依据**：高：可以成为算法创新。
- **潜在问题**：搜索成本高，策略过拟合数据集。
- **规避方案**：小策略空间+验证多数据集泛化。
- **最终可用性**：可用，推荐。
- **来源链接**：https://arxiv.org/abs/2106.07594

### 005. TabR: Tabular Deep Learning Meets Nearest Neighbors

- **排序文件**：`01_PDF论文_按推荐程度排序/005_很高_TabR_Tabular_Deep_Learning_Meets_Nearest_Neighbors.pdf`
- **推荐等级**：很高（score=92）
- **SHA-256**：`FEEC976974E73B858466227A749375AC6C2EC9320F98006ACA0678883469F0BD`
- **来源记录**：文库1#50 50_TabR_Gorishniy_2024_ICLR_openreview.pdf

#### 来源说明：文库1 #50

说明来源：文库1/manifest_extra50_papers.csv
- **标题**：TabR: Tabular Deep Learning Meets Nearest Neighbors
- **年份/会议**：2024 / ICLR
- **方向/领域**：retrieval-augmented tabular DL
- **核心思想**：将深度表格模型与最近邻检索结合。
- **如何改良 scMAE**：scMAE重构/聚类时检索相似细胞作为context，避免只依赖单细胞内部mask。
- **有效性依据**：高：单细胞邻域信息关键。
- **潜在问题**：检索错误会污染表征。
- **规避方案**：使用置信度筛选、batch-aware检索、mutual nearest neighbors。
- **最终可用性**：可用，推荐。
- **来源链接**：https://openreview.net/forum?id=rhgIgTSSxW

### 006. scVGAE: ZINB-Based Variational Graph Autoencoder for Single-Cell RNA-Seq Imputation

- **排序文件**：`01_PDF论文_按推荐程度排序/006_高_scVGAE_ZINB-Based_Variational_Graph_Autoencoder_for_Single-Cell_RNA-Seq_Imputation.pdf`
- **推荐等级**：高（score=85）
- **SHA-256**：`F990B6C875583708FEE7F3E7EFD6C58D29117EBD3042D998928D79A9EB89FC19`
- **来源记录**：文库1#2 02_scVGAE_Inoue_2024_arxiv.pdf

#### 来源说明：文库1 #2

说明来源：文库1/manifest_extra50_papers.csv
- **标题**：scVGAE: ZINB-Based Variational Graph Autoencoder for Single-Cell RNA-Seq Imputation
- **年份/会议**：2024 / arXiv
- **方向/领域**：ZINB + graph VAE
- **核心思想**：将GCN嵌入VAE并使用ZINB损失处理dropout，报告14个数据集中的多数任务有优势。
- **如何改良 scMAE**：把scMAE decoder从MSE扩展为NB/ZINB likelihood，并加细胞图VAE分支约束。
- **有效性依据**：高：scRNA原始count确实离散、过度离散且大量零值。
- **潜在问题**：ZINB对UMI数据未必总是必要；图错误会传播。
- **规避方案**：NB/ZINB可切换，使用graph confidence或edge dropout。
- **最终可用性**：可用。
- **来源链接**：https://arxiv.org/abs/2403.08959

### 007. DinoBloom: A Foundation Model for Generalizable Cell Embeddings in Hematology

- **排序文件**：`01_PDF论文_按推荐程度排序/007_高_DinoBloom_A_Foundation_Model_for_Generalizable_Cell_Embeddings_in_Hematology.pdf`
- **推荐等级**：高（score=85）
- **SHA-256**：`EDEE616745D691CF7DA6A39DCC88572CA269A3AE1BB40384255D62C80FD0B3A2`
- **来源记录**：文库1#3 03_DinoBloom_Koch_2024_MICCAI_arxiv.pdf

#### 来源说明：文库1 #3

说明来源：文库1/manifest_extra50_papers.csv
- **标题**：DinoBloom: A Foundation Model for Generalizable Cell Embeddings in Hematology
- **年份/会议**：2024 / MICCAI / arXiv
- **方向/领域**：medical image SSL / DINOv2
- **核心思想**：针对血液单细胞图像构建DINOv2式foundation model，在外部域移位下保持泛化。
- **如何改良 scMAE**：将DINO/iBOT teacher-student思想迁移到scMAE：引入EMA teacher对不同mask view的一致性约束。
- **有效性依据**：中高：DINOv2说明大规模自蒸馏能学到稳健表征。
- **潜在问题**：scRNA不是图像；强增强可能破坏生物信号。
- **规避方案**：使用生物保真增强：library size perturb、dropout模拟、gene-set mask。
- **最终可用性**：可用。
- **来源链接**：https://arxiv.org/abs/2404.05022

### 008. scCello: Cell-ontology guided transcriptome foundation model

- **排序文件**：`01_PDF论文_按推荐程度排序/008_高_scCello_Cell-ontology_guided_transcriptome_foundation_model.pdf`
- **推荐等级**：高（score=85）
- **SHA-256**：`AB6611C2FB3E29911D23F097A393F92656C31F910BF4E4817168E05CB83F30EE`
- **来源记录**：文库1#4 04_scCello_Yuan_2024_NeurIPS_arxiv.pdf

#### 来源说明：文库1 #4

说明来源：文库1/manifest_extra50_papers.csv
- **标题**：scCello: Cell-ontology guided transcriptome foundation model
- **年份/会议**：2024 / NeurIPS 2024 / arXiv
- **方向/领域**：ontology-guided single-cell FM
- **核心思想**：在masked gene prediction外加入cell-type coherence loss和ontology alignment loss，利用细胞本体图增强表征。
- **如何改良 scMAE**：为scMAE加入cell ontology/marker hierarchy先验；聚类后伪标签映射到层次图。
- **有效性依据**：高：类别层次先验可缓解相近细胞类型混淆。
- **潜在问题**：依赖标签/本体质量，可能引入偏差。
- **规避方案**：弱监督开关：只在高置信伪标签与通用ontology一致时使用。
- **最终可用性**：可用。
- **来源链接**：https://arxiv.org/abs/2408.12373

### 009. LangCell: Language-Cell Pre-training for Cell Identity Understanding

- **排序文件**：`01_PDF论文_按推荐程度排序/009_高_LangCell_Language-Cell_Pre-training_for_Cell_Identity_Understanding.pdf`
- **推荐等级**：高（score=85）
- **SHA-256**：`C7B9166AFF046DA65F379B324AD20F200D4216A0F5A03C6F34E094647C755B21`
- **来源记录**：文库1#5 05_LangCell_Zhao_2024_arxiv.pdf

#### 来源说明：文库1 #5

说明来源：文库1/manifest_extra50_papers.csv
- **标题**：LangCell: Language-Cell Pre-training for Cell Identity Understanding
- **年份/会议**：2024 / arXiv
- **方向/领域**：cell-language multimodal pretraining
- **核心思想**：把单细胞表达和自然语言cell identity描述对齐，支持zero-shot/few-shot语义理解。
- **如何改良 scMAE**：把scMAE embedding与文本描述/marker解释对齐，形成CLIP式cell-text contrastive auxiliary task。
- **有效性依据**：中高：可提升可解释性和少样本标注能力。
- **潜在问题**：文本质量和LLM hallucination会影响训练。
- **规避方案**：用Cell Ontology、PanglaoDB、CellMarker等结构化文本，RAG后人工过滤。
- **最终可用性**：可用。
- **来源链接**：https://arxiv.org/abs/2405.06708

### 010. Celler: A Genomic Language Model for Long-Tailed Single-Cell Annotation

- **排序文件**：`01_PDF论文_按推荐程度排序/010_高_Celler_A_Genomic_Language_Model_for_Long-Tailed_Single-Cell_Annotation.pdf`
- **推荐等级**：高（score=85）
- **SHA-256**：`875169878228564D7224B893979B0D516A178119319A9935E8D6161A3948A49A`
- **来源记录**：文库1#6 06_Celler_Zhao_2025_arxiv.pdf

#### 来源说明：文库1 #6

说明来源：文库1/manifest_extra50_papers.csv
- **标题**：Celler: A Genomic Language Model for Long-Tailed Single-Cell Annotation
- **年份/会议**：2025 / arXiv
- **方向/领域**：long-tail single-cell GLM
- **核心思想**：提出Gaussian Inflation loss和Hard Data Mining，应对疾病数据长尾与稀有细胞。
- **如何改良 scMAE**：在scMAE聚类或伪标签训练阶段引入rare-cell reweighting/uncertainty-aware sampling。
- **有效性依据**：高：scRNA聚类常被大类主导，稀有细胞是一区论文亮点。
- **潜在问题**：过度加权会把噪声点当稀有群体。
- **规避方案**：结合silhouette、local density、marker一致性过滤“伪稀有”。
- **最终可用性**：可用。
- **来源链接**：https://arxiv.org/abs/2504.00020

### 011. scMamba: Scalable Foundation Model for Single-Cell Multi-Omics Integration

- **排序文件**：`01_PDF论文_按推荐程度排序/011_高_scMamba_Scalable_Foundation_Model_for_Single-Cell_Multi-Omics_Integration.pdf`
- **推荐等级**：高（score=85）
- **SHA-256**：`C6E0CBC486E05471F5AE85F15415E0E35DC0B82B3DA7A60799CBF3D4CB37A004`
- **来源记录**：文库1#9 09_scMamba_Yuan_2025_arxiv.pdf

#### 来源说明：文库1 #9

说明来源：文库1/manifest_extra50_papers.csv
- **标题**：scMamba: Scalable Foundation Model for Single-Cell Multi-Omics Integration
- **年份/会议**：2025 / arXiv
- **方向/领域**：state-space model / multi-omics
- **核心思想**：用patch-based genomic tokenization和Mamba/状态空间模型处理高维稀疏multi-omics，避免高变基因预筛带来的信息损失。
- **如何改良 scMAE**：将scMAE Transformer替换或并联Mamba block，按染色体/基因组位置patch tokenization。
- **有效性依据**：中高：Mamba对长序列高效，适合大基因集。
- **潜在问题**：scRNA中基因顺序不是天然语言顺序；位置先验对RNA-only不一定强。
- **规避方案**：只在multi-omics或含ATAC位置数据中使用位置patch；RNA-only用gene-set排序。
- **最终可用性**：可用，作为扩展版方向。
- **来源链接**：https://arxiv.org/abs/2506.20697

### 012. A Survey on Foundation Language Models for Single-cell Biology

- **排序文件**：`01_PDF论文_按推荐程度排序/012_高_A_Survey_on_Foundation_Language_Models_for_Single-cell_Biology.pdf`
- **推荐等级**：高（score=85）
- **SHA-256**：`E89326959164EB0867EAEC104E8F38EA0E29C33BE12E8773BE4FFFCAC96B7CE8`
- **来源记录**：文库1#10 10_Survey_Single_Cell_Foundation_Language_Models_Zhang_2025_ACL.pdf

#### 来源说明：文库1 #10

说明来源：文库1/manifest_extra50_papers.csv
- **标题**：A Survey on Foundation Language Models for Single-cell Biology
- **年份/会议**：2025 / ACL 2025
- **方向/领域**：survey / single-cell FLM
- **核心思想**：系统梳理single-cell language/foundation models。
- **如何改良 scMAE**：作为写Related Work与定位scMAE改良空间的主综述。
- **有效性依据**：高：综述不是算法，但能避免漏引关键FM。
- **潜在问题**：可能过宽，不能替代实验。
- **规避方案**：只用于框架定位和表格。
- **最终可用性**：可用。
- **来源链接**：https://aclanthology.org/2025.acl-long.26.pdf

### 013. Masked Modeling for Single-cell Clustering of scRNA-seq Data

- **排序文件**：`01_PDF论文_按推荐程度排序/013_高_Masked_Modeling_for_Single-cell_Clustering_of_scRNA-seq_Data.pdf`
- **推荐等级**：高（score=85）
- **SHA-256**：`07A1AC9FF48A3521C627E9093309CC2820BEB443701EDB9B514E3C3699D06D53`
- **来源记录**：文库1#12 12_Masked_Modeling_scRNA_2023_OpenReview.pdf

#### 来源说明：文库1 #12

说明来源：文库1/manifest_extra50_papers.csv
- **标题**：Masked Modeling for Single-cell Clustering of scRNA-seq Data
- **年份/会议**：2023 / OpenReview
- **方向/领域**：masked modeling / scRNA clustering
- **核心思想**：与scMAE同源的masked modeling单细胞聚类研究。
- **如何改良 scMAE**：直接对比mask策略、重构目标、聚类头设计。
- **有效性依据**：高：领域直接相关。
- **潜在问题**：创新空间可能被占用。
- **规避方案**：强调与其不同：生物信号mask、图/分布/扩散或不确定性。
- **最终可用性**：可用。
- **来源链接**：https://openreview.net/pdf?id=B8w5mh55ZU

### 014. CICL: scRNA-seq Data Clustering by Cluster-aware Iterative Contrastive Learning

- **排序文件**：`01_PDF论文_按推荐程度排序/014_高_CICL_scRNA-seq_Data_Clustering_by_Cluster-aware_Iterative_Contrastive_Learning.pdf`
- **推荐等级**：高（score=85）
- **SHA-256**：`7973459B80599957CC24386F04EEFAF9FF5B245DC96EE3CC0D5965F8909799B9`
- **来源记录**：文库1#13 13_CICL_2023_scRNA_cluster_aware_iterative_contrastive_learning_arxiv.pdf

#### 来源说明：文库1 #13

说明来源：文库1/manifest_extra50_papers.csv
- **标题**：CICL: scRNA-seq Data Clustering by Cluster-aware Iterative Contrastive Learning
- **年份/会议**：2023 / arXiv
- **方向/领域**：cluster-aware contrastive learning
- **核心思想**：迭代式聚类感知对比学习提升scRNA聚类。
- **如何改良 scMAE**：为scMAE增加cluster prototype contrastive loss，正样本由高置信聚类产生。
- **有效性依据**：高：可直接提升聚类友好性。
- **潜在问题**：早期伪标签错误会自我强化。
- **规避方案**：warm-up后使用，且只采高置信样本；EMA原型。
- **最终可用性**：可用。
- **来源链接**：https://arxiv.org/pdf/2312.16600

### 015. scAGC: Learning Adaptive Cell Graphs with Contrastive Guidance for Single-Cell Clustering

- **排序文件**：`01_PDF论文_按推荐程度排序/015_高_scAGC_Learning_Adaptive_Cell_Graphs_with_Contrastive_Guidance_for_Single-Cell_Clustering.pdf`
- **推荐等级**：高（score=85）
- **SHA-256**：`E63386A5EC7B85E74FE71C2102A6318FFFA9A08A0218C78BE9B161AABDAF5D44`
- **来源记录**：文库1#14 14_scAGC_Li_2025_arxiv.pdf

#### 来源说明：文库1 #14

说明来源：文库1/manifest_extra50_papers.csv
- **标题**：scAGC: Learning Adaptive Cell Graphs with Contrastive Guidance for Single-Cell Clustering
- **年份/会议**：2025 / arXiv
- **方向/领域**：adaptive graph + contrastive + ZINB
- **核心思想**：Gumbel-Softmax动态修正细胞图，ZINB重构并用对比学习稳定拓扑。
- **如何改良 scMAE**：把scMAE静态KNN图改为可学习边；mask重构与图重构共同训练。
- **有效性依据**：高：静态图是scRNA图方法弱点。
- **潜在问题**：可学习图可能崩塌或过平滑。
- **规避方案**：用初始KNN先验、边稀疏约束、contrastive topology consistency。
- **最终可用性**：可用。
- **来源链接**：https://arxiv.org/abs/2508.09180

### 016. Computational Methods for Single-Cell Multi-Omics Integration and Alignment

- **排序文件**：`01_PDF论文_按推荐程度排序/016_高_Computational_Methods_for_Single-Cell_Multi-Omics_Integration_and_Alignment.pdf`
- **推荐等级**：高（score=85）
- **SHA-256**：`BEEF2860B63623DF874A1BA8484E70FD13A370A7FC2418BD4061E9B5A3B54FE8`
- **来源记录**：文库1#15 15_SingleCell_MultiOmics_Integration_Review_Stanojevic_2022_arxiv.pdf

#### 来源说明：文库1 #15

说明来源：文库1/manifest_extra50_papers.csv
- **标题**：Computational Methods for Single-Cell Multi-Omics Integration and Alignment
- **年份/会议**：2022 / Genomics Proteomics Bioinformatics / arXiv
- **方向/领域**：multi-omics integration survey
- **核心思想**：综述多组学对齐、集成、网络/机器翻译思想。
- **如何改良 scMAE**：为scMAE扩展到CITE-seq/scATAC提供设计依据：跨模态mask、shared latent alignment。
- **有效性依据**：中高：多组学是一区热点。
- **潜在问题**：数据集与评测复杂度提高。
- **规避方案**：先用RNA+ADT CITE-seq做小规模验证。
- **最终可用性**：可用。
- **来源链接**：https://arxiv.org/abs/2201.06725

### 017. Interpretable Deep Learning in Single-Cell Omics

- **排序文件**：`01_PDF论文_按推荐程度排序/017_高_Interpretable_Deep_Learning_in_Single-Cell_Omics.pdf`
- **推荐等级**：高（score=85）
- **SHA-256**：`E5BF63B6E8DEC072F16963AA1841E8EB587499B0E14016693E9A20CA041CA55C`
- **来源记录**：文库1#16 16_Interpretable_DL_SingleCell_Omics_Wagle_2024_arxiv.pdf

#### 来源说明：文库1 #16

说明来源：文库1/manifest_extra50_papers.csv
- **标题**：Interpretable Deep Learning in Single-Cell Omics
- **年份/会议**：2024 / arXiv / review
- **方向/领域**：interpretability / single-cell
- **核心思想**：综述单细胞深度模型可解释性与分子调控解释。
- **如何改良 scMAE**：为scMAE加入gene attribution、masked gene importance、pathway enrichment。
- **有效性依据**：高：普通一区计算机论文需要不仅指标好，还要解释。
- **潜在问题**：解释方法可能不稳定。
- **规避方案**：用多种解释一致性、marker基因和通路验证。
- **最终可用性**：可用。
- **来源链接**：https://arxiv.org/abs/2401.06823

### 018. BEiT: BERT Pre-Training of Image Transformers

- **排序文件**：`01_PDF论文_按推荐程度排序/018_高_BEiT_BERT_Pre-Training_of_Image_Transformers.pdf`
- **推荐等级**：高（score=85）
- **SHA-256**：`6455F47A0AD99294C9E1EB8A361DA558F87F9E4446F5D6C0B02C20646A93741C`
- **来源记录**：文库1#17 17_BEiT_Bao_2022_ICLR.pdf

#### 来源说明：文库1 #17

说明来源：文库1/manifest_extra50_papers.csv
- **标题**：BEiT: BERT Pre-Training of Image Transformers
- **年份/会议**：2022 / ICLR
- **方向/领域**：masked image modeling / discrete tokenizer
- **核心思想**：用视觉tokenizer把图像patch转为离散token并做mask prediction。
- **如何改良 scMAE**：将scRNA表达值分桶/向量量化为token，做离散表达预测而非单纯回归。
- **有效性依据**：中高：适合表达强度离散化。
- **潜在问题**：分桶损失数值精度。
- **规避方案**：分类+回归双头，或ordinal regression。
- **最终可用性**：可用。
- **来源链接**：https://arxiv.org/pdf/2106.08254

### 019. data2vec: A General Framework for Self-supervised Learning in Speech, Vision and Language

- **排序文件**：`01_PDF论文_按推荐程度排序/019_高_data2vec_A_General_Framework_for_Self-supervised_Learning_in_Speech,_Vision_and_Language.pdf`
- **推荐等级**：高（score=85）
- **SHA-256**：`9AB545577568481959F5BA5596385CDF1E6344F5B1A5A821A9779996EBEB0E18`
- **来源记录**：文库1#20 20_data2vec_Baevski_2022_ICML_arxiv.pdf

#### 来源说明：文库1 #20

说明来源：文库1/manifest_extra50_papers.csv
- **标题**：data2vec: A General Framework for Self-supervised Learning in Speech, Vision and Language
- **年份/会议**：2022 / ICML
- **方向/领域**：teacher target / modality-agnostic SSL
- **核心思想**：用teacher网络产生上下文化latent target，统一语音/图像/NLP自监督。
- **如何改良 scMAE**：scMAE可不直接重构噪声表达，而预测teacher latent representation。
- **有效性依据**：高：对噪声数据更稳健。
- **潜在问题**：teacher collapse或目标过平滑。
- **规避方案**：EMA teacher+stop-gradient+保留少量表达重构辅助。
- **最终可用性**：可用。
- **来源链接**：https://arxiv.org/pdf/2202.03555

### 020. MultiMAE: Multi-modal Multi-task Masked Autoencoders

- **排序文件**：`01_PDF论文_按推荐程度排序/020_高_MultiMAE_Multi-modal_Multi-task_Masked_Autoencoders.pdf`
- **推荐等级**：高（score=85）
- **SHA-256**：`1806F7C832716060E0850E83AEAE6404D85D832434E65282E892736698F343D3`
- **来源记录**：文库1#23 23_MultiMAE_Bachmann_2022_ECCV.pdf

#### 来源说明：文库1 #23

说明来源：文库1/manifest_extra50_papers.csv
- **标题**：MultiMAE: Multi-modal Multi-task Masked Autoencoders
- **年份/会议**：2022 / ECCV
- **方向/领域**：multi-modal MAE
- **核心思想**：多模态/多任务masked autoencoding。
- **如何改良 scMAE**：scMAE扩展RNA+ADT+ATAC，随机mask不同模态并交叉重构。
- **有效性依据**：高：若有多模态数据，创新强。
- **潜在问题**：多模态缺失与batch效应复杂。
- **规避方案**：modality dropout、shared+private latent、OT alignment。
- **最终可用性**：可用。
- **来源链接**：https://www.ecva.net/papers/eccv_2022/papers_ECCV/papers/136970341.pdf

### 021. AudioMAE: Masked Autoencoders that Listen

- **排序文件**：`01_PDF论文_按推荐程度排序/021_高_AudioMAE_Masked_Autoencoders_that_Listen.pdf`
- **推荐等级**：高（score=85）
- **SHA-256**：`410E0515933E535FC96F4EA93DA828B31A116635B0F498DD9DDF1B24526A3762`
- **来源记录**：文库1#25 25_AudioMAE_Huang_2022_NeurIPS.pdf

#### 来源说明：文库1 #25

说明来源：文库1/manifest_extra50_papers.csv
- **标题**：AudioMAE: Masked Autoencoders that Listen
- **年份/会议**：2022 / NeurIPS
- **方向/领域**：spectrogram MAE / local window attention
- **核心思想**：音频谱图局部时频相关，decoder用local window attention。
- **如何改良 scMAE**：scMAE可引入gene-set/local pathway window attention，不必全局注意所有基因。
- **有效性依据**：中高：基因模块具有局部功能相关性。
- **潜在问题**：gene-set划分不唯一。
- **规避方案**：MSigDB/GO多视图gene-set attention并做消融。
- **最终可用性**：可用。
- **来源链接**：https://proceedings.neurips.cc/paper_files/paper/2022/file/b89d5e209990b19e33b418e14f323998-Paper-Conference.pdf

### 022. I-JEPA: Self-Supervised Learning from Images with Joint-Embedding Predictive Architecture

- **排序文件**：`01_PDF论文_按推荐程度排序/022_高_I-JEPA_Self-Supervised_Learning_from_Images_with_Joint-Embedding_Predictive_Architecture.pdf`
- **推荐等级**：高（score=85）
- **SHA-256**：`EDDBDC093EB4D48662BCF4FBD1C6735CCD606205CC9070977CE317D68AFF3941`
- **来源记录**：文库1#26 26_IJEPA_Assran_2023_CVPR_arxiv.pdf

#### 来源说明：文库1 #26

说明来源：文库1/manifest_extra50_papers.csv
- **标题**：I-JEPA: Self-Supervised Learning from Images with Joint-Embedding Predictive Architecture
- **年份/会议**：2023 / CVPR
- **方向/领域**：non-generative SSL / block prediction
- **核心思想**：从context block预测target block的表示而非像素。
- **如何改良 scMAE**：scMAE从“恢复表达值”变成“预测被mask基因模块的语义表示”，减小噪声重构。
- **有效性依据**：高：更偏语义聚类。
- **潜在问题**：下游细粒度基因值可能弱化。
- **规避方案**：保留轻量值重构辅助。
- **最终可用性**：可用。
- **来源链接**：https://arxiv.org/abs/2301.08243

### 023. MaskGIT: Masked Generative Image Transformer

- **排序文件**：`01_PDF论文_按推荐程度排序/023_高_MaskGIT_Masked_Generative_Image_Transformer.pdf`
- **推荐等级**：高（score=85）
- **SHA-256**：`DD716A8A8943955B5A1BABD0EE41D57F6D5CD81206B59ED8EDA4C76038C81ED1`
- **来源记录**：文库1#29 29_MaskGIT_Chang_2022_CVPR_arxiv.pdf

#### 来源说明：文库1 #29

说明来源：文库1/manifest_extra50_papers.csv
- **标题**：MaskGIT: Masked Generative Image Transformer
- **年份/会议**：2022 / CVPR
- **方向/领域**：iterative masked generation
- **核心思想**：并行预测masked tokens并迭代refine，避免自回归慢与顺序偏置。
- **如何改良 scMAE**：scMAE聚类前做iterative expression refinement：先高置信恢复，再迭代低置信基因。
- **有效性依据**：高：符合无序表达集合。
- **潜在问题**：迭代生成可能放大错误。
- **规避方案**：置信度阈值+只在预训练使用，不直接替换原始数据。
- **最终可用性**：可用。
- **来源链接**：https://arxiv.org/abs/2202.04200

### 024. Flow Matching for Generative Modeling

- **排序文件**：`01_PDF论文_按推荐程度排序/024_高_Flow_Matching_for_Generative_Modeling.pdf`
- **推荐等级**：高（score=85）
- **SHA-256**：`5D3FEB9423BEE52D1FC12892ECB584BCBC578A1376972DFFFB6FB77807A46387`
- **来源记录**：文库1#33 33_Flow_Matching_Lipman_2023_ICLR_openreview.pdf

#### 来源说明：文库1 #33

说明来源：文库1/manifest_extra50_papers.csv
- **标题**：Flow Matching for Generative Modeling
- **年份/会议**：2023 / ICLR
- **方向/领域**：flow matching / OT paths
- **核心思想**：用固定概率路径的向量场回归训练连续流，OT路径更高效。
- **如何改良 scMAE**：scMAE latent空间做OT-flow regularization：不同batch/condition流向统一cell manifold。
- **有效性依据**：中高：批次校正/扰动预测有潜力。
- **潜在问题**：理论复杂、实现成本高。
- **规避方案**：先用OT alignment loss，不直接训练完整flow。
- **最终可用性**：可用。
- **来源链接**：https://openreview.net/forum?id=PqvMRDCJT9t

### 025. Anomaly Transformer

- **排序文件**：`01_PDF论文_按推荐程度排序/025_高_Anomaly_Transformer.pdf`
- **推荐等级**：高（score=85）
- **SHA-256**：`B3308D4695DD8F5F8D874AD24076B0AC92E938BF4E5C6B582B2AE585A820029F`
- **来源记录**：文库1#37 37_Anomaly_Transformer_Xu_2022_ICLR.pdf

#### 来源说明：文库1 #37

说明来源：文库1/manifest_extra50_papers.csv
- **标题**：Anomaly Transformer
- **年份/会议**：2022 / ICLR
- **方向/领域**：anomaly association / signal outliers
- **核心思想**：用关联差异识别异常点。
- **如何改良 scMAE**：识别scMAE训练中的doublets、低质量细胞、极端outliers并降低权重。
- **有效性依据**：中高：提高鲁棒性。
- **潜在问题**：稀有细胞可能被误判异常。
- **规避方案**：加入marker一致性和局部密度保护稀有群。
- **最终可用性**：可用。
- **来源链接**：https://openreview.net/forum?id=LzQQ89U1qm_

### 026. BGRL: Large-Scale Representation Learning on Graphs via Bootstrapping

- **排序文件**：`01_PDF论文_按推荐程度排序/026_高_BGRL_Large-Scale_Representation_Learning_on_Graphs_via_Bootstrapping.pdf`
- **推荐等级**：高（score=85）
- **SHA-256**：`E5DAD781D9591D20F665B2BCE39F35218B0FA8DC53AD6B5C7153711961C82EB0`
- **来源记录**：文库1#38 38_BGRL_Thakoor_2022_ICLR.pdf

#### 来源说明：文库1 #38

说明来源：文库1/manifest_extra50_papers.csv
- **标题**：BGRL: Large-Scale Representation Learning on Graphs via Bootstrapping
- **年份/会议**：2022 / ICLR
- **方向/领域**：graph SSL without negatives
- **核心思想**：不用负样本，通过两种图增强的bootstrap学习大图表征。
- **如何改良 scMAE**：scMAE细胞图分支使用BGRL，避免负样本误选相近细胞类型。
- **有效性依据**：高：scRNA负样本定义困难。
- **潜在问题**：BYOL类方法有坍塌风险。
- **规避方案**：EMA target、predictor、variance/covariance正则。
- **最终可用性**：可用。
- **来源链接**：https://arxiv.org/abs/2102.06514

### 027. Graph Barlow Twins

- **排序文件**：`01_PDF论文_按推荐程度排序/027_高_Graph_Barlow_Twins.pdf`
- **推荐等级**：高（score=85）
- **SHA-256**：`866E6229C93BB28050272EE5A6F821D613CD7E85F55E1B5CF473943835E7D1E8`
- **来源记录**：文库1#40 40_Graph_Barlow_Twins_Bielak_2022_KBS_arxiv.pdf

#### 来源说明：文库1 #40

说明来源：文库1/manifest_extra50_papers.csv
- **标题**：Graph Barlow Twins
- **年份/会议**：2022 / Knowledge-Based Systems
- **方向/领域**：negative-free graph SSL
- **核心思想**：用cross-correlation冗余抑制取代负样本。
- **如何改良 scMAE**：scMAE加入Barlow Twins loss约束两个mask view的embedding去冗余且一致。
- **有效性依据**：高：适合无标签聚类。
- **潜在问题**：过强decorrelation可能破坏相关基因模块。
- **规避方案**：只作用cell embedding，不直接作用gene-level expression。
- **最终可用性**：可用。
- **来源链接**：https://arxiv.org/abs/2106.02466

### 028. Graphormer

- **排序文件**：`01_PDF论文_按推荐程度排序/028_高_Graphormer.pdf`
- **推荐等级**：高（score=85）
- **SHA-256**：`1F3B7D34FE425ACB9DA5E955E4E7129CD44730AEB5A357E6A3E509E45C681E86`
- **来源记录**：文库1#42 42_Graphormer_Ying_2021_NeurIPS.pdf

#### 来源说明：文库1 #42

说明来源：文库1/manifest_extra50_papers.csv
- **标题**：Graphormer
- **年份/会议**：2021 / NeurIPS
- **方向/领域**：graph transformer structural bias
- **核心思想**：把中心性、最短路径、边特征编码进Transformer attention。
- **如何改良 scMAE**：将细胞图距离/共表达网络距离作为attention bias注入scMAE。
- **有效性依据**：高：结构先验明确。
- **潜在问题**：图距离对噪声敏感。
- **规避方案**：用多图平均和distance clipping。
- **最终可用性**：可用。
- **来源链接**：https://proceedings.neurips.cc/paper/2021/hash/f1c1592588411002af340cbaedd6fc33-Abstract.html

### 029. Deep Adaptive Fuzzy Clustering for Evolutionary Unsupervised Representation Learning

- **排序文件**：`01_PDF论文_按推荐程度排序/029_高_Deep_Adaptive_Fuzzy_Clustering_for_Evolutionary_Unsupervised_Representation_Learning.pdf`
- **推荐等级**：高（score=85）
- **SHA-256**：`36D60AC770E78AAA81661F097C02BA71C05D3D2E40F6AFC0E1B6F5EE73BFB9AE`
- **来源记录**：文库1#45 45_Deep_Adaptive_Fuzzy_Clustering_Tan_2021_arxiv.pdf

#### 来源说明：文库1 #45

说明来源：文库1/manifest_extra50_papers.csv
- **标题**：Deep Adaptive Fuzzy Clustering for Evolutionary Unsupervised Representation Learning
- **年份/会议**：2021 / arXiv / pattern recognition
- **方向/领域**：fuzzy clustering
- **核心思想**：把深度表征与自适应模糊聚类联合训练，membership表示不确定归属。
- **如何改良 scMAE**：scMAE clustering head输出fuzzy memberships，处理连续谱系/过渡细胞。
- **有效性依据**：高：单细胞过渡状态天然模糊。
- **潜在问题**：模糊聚类可能降低硬标签指标。
- **规避方案**：报告ARI/NMI外增加soft silhouette、trajectory consistency；最终硬标签取最大membership。
- **最终可用性**：可用。
- **来源链接**：https://arxiv.org/abs/2103.17086

### 030. Fuzzy Rough Sets Based on Fuzzy Quantification

- **排序文件**：`01_PDF论文_按推荐程度排序/030_高_Fuzzy_Rough_Sets_Based_on_Fuzzy_Quantification.pdf`
- **推荐等级**：高（score=85）
- **SHA-256**：`BAEC986C62573F3D0FC7E622020482A73B762A67C72E9B906BE66FC6549AF3E8`
- **来源记录**：文库1#46 46_Fuzzy_Rough_Quantification_Theerens_2023_FSS_arxiv.pdf

#### 来源说明：文库1 #46

说明来源：文库1/manifest_extra50_papers.csv
- **标题**：Fuzzy Rough Sets Based on Fuzzy Quantification
- **年份/会议**：2023 / Fuzzy Sets and Systems
- **方向/领域**：fuzzy rough sets / noise robustness
- **核心思想**：用模糊量词增强fuzzy rough set对噪声的鲁棒性。
- **如何改良 scMAE**：对scMAE聚类边界细胞构建lower/upper approximation：核心细胞强监督，边界细胞弱监督。
- **有效性依据**：中高：可解释且适合噪声/模糊边界。
- **潜在问题**：实现与深度模型结合不直接。
- **规避方案**：只作为loss weighting/置信度建模层。
- **最终可用性**：可用。
- **来源链接**：https://arxiv.org/abs/2212.04327

### 031. scPilot: Large Language Model Reasoning Toward Automated Single-Cell Analysis and Discovery

- **排序文件**：`01_PDF论文_按推荐程度排序/031_中高_可用_scPilot_Large_Language_Model_Reasoning_Toward_Automated_Single-Cell_Analysis_and_Discovery.pdf`
- **推荐等级**：中高/可用（score=72）
- **SHA-256**：`7295059CBEAF9556A99F739065D7B9A3191DE66B78254D9FD0CE246D152C73F9`
- **来源记录**：文库1#8 08_scPilot_Gao_2026_arxiv.pdf

#### 来源说明：文库1 #8

说明来源：文库1/manifest_extra50_papers.csv
- **标题**：scPilot: Large Language Model Reasoning Toward Automated Single-Cell Analysis and Discovery
- **年份/会议**：2026 / arXiv / OpenReview
- **方向/领域**：LLM tool-use for single-cell reasoning
- **核心思想**：LLM直接检查单细胞数据、调用工具并迭代解释注释/轨迹/调控推断。
- **如何改良 scMAE**：不是训练backbone，而是作为实验与解释模块：自动生成marker验证、失败案例分析、消融建议。
- **有效性依据**：中：对模型性能本身帮助间接，但提高论文完整性。
- **潜在问题**：LLM输出不可作为最终证据。
- **规避方案**：把LLM只用于候选假设，最终以统计检验和marker数据库验证。
- **最终可用性**：可用，适合作为可解释分析附加模块。
- **来源链接**：https://arxiv.org/abs/2602.11609

### 032. ChatCell: Natural Language Interface for Single-Cell Analysis

- **排序文件**：`01_PDF论文_按推荐程度排序/032_中高_可用_ChatCell_Natural_Language_Interface_for_Single-Cell_Analysis.pdf`
- **推荐等级**：中高/可用（score=72）
- **SHA-256**：`38F83AE62BBE22A2283C9C26625BAAB4DD0D4814E2B52978123074464505EC43`
- **来源记录**：文库1#11 11_ChatCell_2024_arxiv.pdf

#### 来源说明：文库1 #11

说明来源：文库1/manifest_extra50_papers.csv
- **标题**：ChatCell: Natural Language Interface for Single-Cell Analysis
- **年份/会议**：2024 / arXiv
- **方向/领域**：LLM + single-cell interface
- **核心思想**：用自然语言驱动单细胞分析任务。
- **如何改良 scMAE**：可用于构建scMAE结果解释与自动报告模块，强化工程系统性。
- **有效性依据**：中：作为辅助工具，不是核心算法。
- **潜在问题**：评审可能认为偏应用。
- **规避方案**：放在补充实验/可解释性，不作为主贡献。
- **最终可用性**：可用。
- **来源链接**：https://www.arxiv.org/pdf/2402.08303v2.pdf

### 033. iBOT: Image BERT Pre-Training with Online Tokenizer

- **排序文件**：`01_PDF论文_按推荐程度排序/033_中高_可用_iBOT_Image_BERT_Pre-Training_with_Online_Tokenizer.pdf`
- **推荐等级**：中高/可用（score=72）
- **SHA-256**：`E76C32B28234434443CBA27ABDA926EC7C69C64AF240D2D7AA78427131C4273D`
- **来源记录**：文库1#21 21_iBOT_Zhou_2022_ICLR_openreview.pdf

#### 来源说明：文库1 #21

说明来源：文库1/manifest_extra50_papers.csv
- **标题**：iBOT: Image BERT Pre-Training with Online Tokenizer
- **年份/会议**：2022 / ICLR
- **方向/领域**：online tokenizer / self-distillation MIM
- **核心思想**：在线tokenizer+masked self-distillation，无需离线dVAE。
- **如何改良 scMAE**：为scMAE建立online gene-expression tokenizer，让token边训练边演化。
- **有效性依据**：中高。
- **潜在问题**：在线token可能不稳定。
- **规避方案**：EMA teacher、token entropy正则、原型均衡。
- **最终可用性**：可用。
- **来源链接**：https://openreview.net/pdf?id=ydopy-e6Dg

### 034. Masked Siamese Networks for Label-Efficient Learning

- **排序文件**：`01_PDF论文_按推荐程度排序/034_中高_可用_Masked_Siamese_Networks_for_Label-Efficient_Learning.pdf`
- **推荐等级**：中高/可用（score=72）
- **SHA-256**：`91B737FE3FFBBD292B51B82ED4466DE3BF0BF7D5FF27E2ABAF4F4309E60186F8`
- **来源记录**：文库1#22 22_MSN_Assran_2022_ECCV.pdf

#### 来源说明：文库1 #22

说明来源：文库1/manifest_extra50_papers.csv
- **标题**：Masked Siamese Networks for Label-Efficient Learning
- **年份/会议**：2022 / ECCV
- **方向/领域**：masked siamese SSL
- **核心思想**：masked view与unmasked/augmented view做siamese一致性。
- **如何改良 scMAE**：scMAE加两个cell views：masked expression view与dropout/noise view一致。
- **有效性依据**：中高。
- **潜在问题**：增强若不保真会抹掉稀有信号。
- **规避方案**：只用生物合理扰动，保留marker genes。
- **最终可用性**：可用。
- **来源链接**：https://www.ecva.net/papers/eccv_2022/papers_ECCV/papers/136910442.pdf

### 035. VideoMAE: Masked Autoencoders are Data-Efficient Learners for Self-Supervised Video Pre-Training

- **排序文件**：`01_PDF论文_按推荐程度排序/035_中高_可用_VideoMAE_Masked_Autoencoders_are_Data-Efficient_Learners_for_Self-Supervised_Video_Pre-Training.pdf`
- **推荐等级**：中高/可用（score=72）
- **SHA-256**：`C074658068663C36684FA873F920909CADF56747F166C689088A403487352ED4`
- **来源记录**：文库1#24 24_VideoMAE_Tong_2022_NeurIPS.pdf

#### 来源说明：文库1 #24

说明来源：文库1/manifest_extra50_papers.csv
- **标题**：VideoMAE: Masked Autoencoders are Data-Efficient Learners for Self-Supervised Video Pre-Training
- **年份/会议**：2022 / NeurIPS
- **方向/领域**：temporal MIM
- **核心思想**：高mask率视频MAE，利用时空冗余。
- **如何改良 scMAE**：把细胞发育轨迹/伪时间当序列，做trajectory-aware mask。
- **有效性依据**：中：适合发育数据。
- **潜在问题**：大多数数据没有真实时间。
- **规避方案**：只在trajectory benchmark中做附加实验。
- **最终可用性**：部分可用。
- **来源链接**：https://proceedings.neurips.cc/paper_files/paper/2022/file/416f9cb3276121c42eebb86352a4354a-Paper-Conference.pdf

### 036. DINOv2: Learning Robust Visual Features without Supervision

- **排序文件**：`01_PDF论文_按推荐程度排序/036_中高_可用_DINOv2_Learning_Robust_Visual_Features_without_Supervision.pdf`
- **推荐等级**：中高/可用（score=72）
- **SHA-256**：`3ADD205547F3ECE5FB1AC7A91720AA9FABF02D9CDE29A096643FF022449F3E91`
- **来源记录**：文库1#27 27_DINOv2_Oquab_2023_arxiv.pdf

#### 来源说明：文库1 #27

说明来源：文库1/manifest_extra50_papers.csv
- **标题**：DINOv2: Learning Robust Visual Features without Supervision
- **年份/会议**：2023 / arXiv / CV foundation model
- **方向/领域**：large-scale SSL / distillation
- **核心思想**：大规模自监督+数据清洗+蒸馏得到通用视觉特征。
- **如何改良 scMAE**：scMAE借鉴数据curation和teacher-student distillation；对batch/organ数据做balanced sampling。
- **有效性依据**：中高。
- **潜在问题**：需要大量数据和算力。
- **规避方案**：先做小规模teacher-student，报告计算量。
- **最终可用性**：可用。
- **来源链接**：https://arxiv.org/abs/2304.07193

### 037. ImageBind: One Embedding Space To Bind Them All

- **排序文件**：`01_PDF论文_按推荐程度排序/037_中高_可用_ImageBind_One_Embedding_Space_To_Bind_Them_All.pdf`
- **推荐等级**：中高/可用（score=72）
- **SHA-256**：`0DCF74B737B46D82C304FB88FD16DFB3B1911975AC9439A23B9B0E82049E0FFB`
- **来源记录**：文库1#28 28_ImageBind_Girdhar_2023_CVPR_arxiv.pdf

#### 来源说明：文库1 #28

说明来源：文库1/manifest_extra50_papers.csv
- **标题**：ImageBind: One Embedding Space To Bind Them All
- **年份/会议**：2023 / CVPR
- **方向/领域**：multimodal alignment
- **核心思想**：通过成对数据把多种模态绑定到同一embedding空间。
- **如何改良 scMAE**：scMAE可对齐表达、文本、图像/空间组学或蛋白模态，做统一cell embedding。
- **有效性依据**：中：跨模态创新大。
- **潜在问题**：单纯RNA任务可能过重。
- **规避方案**：作为扩展实验；主方法保持RNA可用。
- **最终可用性**：可用。
- **来源链接**：https://arxiv.org/abs/2305.05665

### 038. Latent Diffusion Models

- **排序文件**：`01_PDF论文_按推荐程度排序/038_中高_可用_Latent_Diffusion_Models.pdf`
- **推荐等级**：中高/可用（score=72）
- **SHA-256**：`46EDE043A8DC07CA1F0F445620523FE1AD8B2436BD83856A3835612A47E9F79E`
- **来源记录**：文库1#30 30_Latent_Diffusion_Rombach_2022_CVPR_arxiv.pdf

#### 来源说明：文库1 #30

说明来源：文库1/manifest_extra50_papers.csv
- **标题**：Latent Diffusion Models
- **年份/会议**：2022 / CVPR
- **方向/领域**：latent generative model
- **核心思想**：先用AE压缩，再在latent空间扩散，降低成本。
- **如何改良 scMAE**：在scMAE latent space做diffusion/flow增强稀有细胞或建模扰动。
- **有效性依据**：中高。
- **潜在问题**：生成数据可能不真实。
- **规避方案**：只用于正则/对比增强，不直接作为真样本评测。
- **最终可用性**：可用。
- **来源链接**：https://arxiv.org/abs/2112.10752

### 039. TS2Vec: Towards Universal Representation of Time Series

- **排序文件**：`01_PDF论文_按推荐程度排序/039_中高_可用_TS2Vec_Towards_Universal_Representation_of_Time_Series.pdf`
- **推荐等级**：中高/可用（score=72）
- **SHA-256**：`24783992080266234820E077F29A49DE36D44FC71BF4C9553FFE6AFD1207F4FF`
- **来源记录**：文库1#34 34_TS2Vec_Yue_2022_AAAI_arxiv.pdf

#### 来源说明：文库1 #34

说明来源：文库1/manifest_extra50_papers.csv
- **标题**：TS2Vec: Towards Universal Representation of Time Series
- **年份/会议**：2022 / AAAI
- **方向/领域**：time-series contrastive
- **核心思想**：分层时间序列对比学习，学习不同尺度表征。
- **如何改良 scMAE**：用于伪时间/发育轨迹数据：细胞按轨迹排序后做多尺度一致性。
- **有效性依据**：中。
- **潜在问题**：排序误差影响大。
- **规避方案**：只在trajectory数据附加，不影响主任务。
- **最终可用性**：部分可用。
- **来源链接**：https://arxiv.org/abs/2106.10466

### 040. TimesNet: Temporal 2D-Variation Modeling for General Time Series Analysis

- **排序文件**：`01_PDF论文_按推荐程度排序/040_中高_可用_TimesNet_Temporal_2D-Variation_Modeling_for_General_Time_Series_Analysis.pdf`
- **推荐等级**：中高/可用（score=72）
- **SHA-256**：`12DA8DA7BF1B3DB1BDBD55EF0358967C6B49D84F6183F7DA7176C57CECD4D513`
- **来源记录**：文库1#35 35_TimesNet_Wu_2023_ICLR.pdf

#### 来源说明：文库1 #35

说明来源：文库1/manifest_extra50_papers.csv
- **标题**：TimesNet: Temporal 2D-Variation Modeling for General Time Series Analysis
- **年份/会议**：2023 / ICLR
- **方向/领域**：signal periodicity / 2D variation
- **核心思想**：把1D时间序列变成2D周期结构捕捉多尺度变化。
- **如何改良 scMAE**：把基因表达按gene-set/通路分块形成2D“信号图”，用于mask块重构。
- **有效性依据**：中：信号原理启发。
- **潜在问题**：基因排序人为。
- **规避方案**：使用多种排序策略并ensemble。
- **最终可用性**：探索性可用。
- **来源链接**：https://openreview.net/forum?id=ju_Uqw384Oq

### 041. Non-stationary Transformers

- **排序文件**：`01_PDF论文_按推荐程度排序/041_中高_可用_Non-stationary_Transformers.pdf`
- **推荐等级**：中高/可用（score=72）
- **SHA-256**：`21E555C46592EDC51207DCE42AC666B9F22A568BB4CA219D3CEE3EAE22287638`
- **来源记录**：文库1#36 36_Nonstationary_Transformers_Liu_2022_NeurIPS.pdf

#### 来源说明：文库1 #36

说明来源：文库1/manifest_extra50_papers.csv
- **标题**：Non-stationary Transformers
- **年份/会议**：2022 / NeurIPS
- **方向/领域**：distribution shift / de-stationary attention
- **核心思想**：处理非平稳时间序列，用去平稳/再平稳机制。
- **如何改良 scMAE**：scRNA batch/domain shift可类比非平稳；加入domain-specific normalization和attention debias。
- **有效性依据**：中高。
- **潜在问题**：类比不完全。
- **规避方案**：只采用归一化/域偏移校正机制。
- **最终可用性**：可用。
- **来源链接**：https://openreview.net/forum?id=ucNDIDRNjjv

### 042. MaskGAE: Masked Graph Modeling Meets Graph Autoencoders

- **排序文件**：`01_PDF论文_按推荐程度排序/042_中高_可用_MaskGAE_Masked_Graph_Modeling_Meets_Graph_Autoencoders.pdf`
- **推荐等级**：中高/可用（score=72）
- **SHA-256**：`7BBD9AA08B5113AE0E8B1C56F8398A9B3CD7E101BC62FEFE06486D5BE613120E`
- **来源记录**：文库1#39 39_MaskGAE_Li_2023_KDD_arxiv.pdf

#### 来源说明：文库1 #39

说明来源：文库1/manifest_extra50_papers.csv
- **标题**：MaskGAE: Masked Graph Modeling Meets Graph Autoencoders
- **年份/会议**：2023 / KDD
- **方向/领域**：masked graph autoencoder
- **核心思想**：mask edges并重构，理论连接GAE与contrastive learning。
- **如何改良 scMAE**：把scMAE从gene mask扩展到cell graph edge mask，联合重构表达与邻接。
- **有效性依据**：高。
- **潜在问题**：边mask过多导致图断裂。
- **规避方案**：自适应edge mask比例，保持连通性。
- **最终可用性**：可用。
- **来源链接**：https://arxiv.org/abs/2205.10053

### 043. GraphGPS: General Powerful Scalable Graph Transformers

- **排序文件**：`01_PDF论文_按推荐程度排序/043_中高_可用_GraphGPS_General_Powerful_Scalable_Graph_Transformers.pdf`
- **推荐等级**：中高/可用（score=72）
- **SHA-256**：`B55B6C7D74AA578025A31CF1D546FE67EBF22DCE69A2B7D867556DFF06E77A7C`
- **来源记录**：文库1#41 41_GraphGPS_Rampasek_2022_NeurIPS_arxiv.pdf

#### 来源说明：文库1 #41

说明来源：文库1/manifest_extra50_papers.csv
- **标题**：GraphGPS: General Powerful Scalable Graph Transformers
- **年份/会议**：2022 / NeurIPS
- **方向/领域**：graph transformer / positional encoding
- **核心思想**：结合局部message passing和全局attention，使用图位置编码。
- **如何改良 scMAE**：scMAE图分支用local cell graph + global transformer，平衡局部邻域和全局细胞类型。
- **有效性依据**：中高。
- **潜在问题**：模型复杂。
- **规避方案**：轻量化：只在latent cell graph上2层GraphGPS。
- **最终可用性**：可用。
- **来源链接**：https://arxiv.org/abs/2205.12454

### 044. Mole-BERT: Rethinking Pre-training Graph Neural Networks for Molecules

- **排序文件**：`01_PDF论文_按推荐程度排序/044_中高_可用_Mole-BERT_Rethinking_Pre-training_Graph_Neural_Networks_for_Molecules.pdf`
- **推荐等级**：中高/可用（score=72）
- **SHA-256**：`FDFDAE5C7E557803D2EB725CE730FDAF474F3D94E114984E321D4B4C62DBE146`
- **来源记录**：文库1#43 43_MoleBERT_Xia_2023_ICLR.pdf

#### 来源说明：文库1 #43

说明来源：文库1/manifest_extra50_papers.csv
- **标题**：Mole-BERT: Rethinking Pre-training Graph Neural Networks for Molecules
- **年份/会议**：2023 / ICLR
- **方向/领域**：graph tokenization / masked graph modeling
- **核心思想**：分子图中用上下文感知tokenizer和masked atom prediction。
- **如何改良 scMAE**：把gene module/tokenizer用于基因网络；mask基因模块而非随机单基因。
- **有效性依据**：中高。
- **潜在问题**：基因网络不如分子图确定。
- **规避方案**：使用多个PPI/GRN数据库交集或加权图。
- **最终可用性**：可用。
- **来源链接**：https://openreview.net/forum?id=jevY-DtiZTR

### 045. SAINT: Improved Neural Networks for Tabular Data via Row Attention and Contrastive Pre-training

- **排序文件**：`01_PDF论文_按推荐程度排序/045_中高_可用_SAINT_Improved_Neural_Networks_for_Tabular_Data_via_Row_Attention_and_Contrastive_Pre-training.pdf`
- **推荐等级**：中高/可用（score=72）
- **SHA-256**：`A1B13AF9DF1325634301B7A15DE111BD7C759C025A091EE71682E8E7813858A0`
- **来源记录**：文库1#47 47_SAINT_Somepalli_2021_NeurIPS_Workshop.pdf

#### 来源说明：文库1 #47

说明来源：文库1/manifest_extra50_papers.csv
- **标题**：SAINT: Improved Neural Networks for Tabular Data via Row Attention and Contrastive Pre-training
- **年份/会议**：2021 / NeurIPS Workshop
- **方向/领域**：tabular transformer / contrastive
- **核心思想**：表格数据行列attention和contrastive pretraining。
- **如何改良 scMAE**：scRNA是高维表格；增加cell-cell row attention与gene-gene column attention双轴建模。
- **有效性依据**：中高。
- **潜在问题**：全量row attention O(n^2)成本高。
- **规避方案**：用mini-batch neighbor attention。
- **最终可用性**：可用。
- **来源链接**：https://table-representation-learning.github.io/assets/papers/saint_improved_neural_networks.pdf

### 046. TabPFN: A Transformer That Solves Small Tabular Classification Problems in a Second

- **排序文件**：`01_PDF论文_按推荐程度排序/046_中高_可用_TabPFN_A_Transformer_That_Solves_Small_Tabular_Classification_Problems_in_a_Second.pdf`
- **推荐等级**：中高/可用（score=72）
- **SHA-256**：`97796B6E6B134B200555A7EBBACEA2F6DAECB9AB309279BA50D915D5AC5C9DF6`
- **来源记录**：文库1#49 49_TabPFN_Hollmann_2023_ICLR_arxiv.pdf

#### 来源说明：文库1 #49

说明来源：文库1/manifest_extra50_papers.csv
- **标题**：TabPFN: A Transformer That Solves Small Tabular Classification Problems in a Second
- **年份/会议**：2023 / ICLR
- **方向/领域**：tabular foundation model / in-context learning
- **核心思想**：PFN在合成任务上离线学习贝叶斯推断，推理时in-context完成小表分类。
- **如何改良 scMAE**：可尝试“synthetic single-cell priors”：用模拟scRNA任务预训练scMAE，使其快速适配小数据集。
- **有效性依据**：中高，思想新颖。
- **潜在问题**：合成先验与真实scRNA差距。
- **规避方案**：用Splatter/scDesign模拟并和真实预训练混合。
- **最终可用性**：可用。
- **来源链接**：https://arxiv.org/abs/2207.01848

### 047. Self-Guided Masked Autoencoder

- **排序文件**：`01_PDF论文_按推荐程度排序/047_可用_Self-Guided_Masked_Autoencoder.pdf`
- **推荐等级**：可用（score=70）
- **SHA-256**：`A3929DCBBB9B96B04DD395AA9F9E88112E6536C8120E380CB85B1CAC0798B5B6`
- **来源记录**：文库2#2 02_Self-Guided_MAE.pdf

#### 来源说明：文库2 #2

说明来源：文库2/scmae_improvement_reports/02_SelfGuided_MAE.txt
- **标题**：Self-Guided Masked Autoencoder
- **年份/会议**：2024 / NeurIPS
- **备注**：开放PDF。

Self‑Guided Masked Autoencoder (Self‑Guided MAE)

1. **Efficacy & principle** – This NeurIPS 2024 paper observes that standard MAE pre‑training intrinsically clusters patches based on visual patterns.  It proposes a self‑guided masking strategy that uses early learned patch‑level clusters to inform which patches to mask, accelerating training and improving downstream accuracy【944398912377764†L14-L23】.
2. **Integration strategy** – For scMAE, a similar self‑guided strategy could cluster genes or cell features during early epochs and use these clusters to guide masking.  Instead of random masking, genes within the same cluster could be masked together, encouraging the model to learn inter‑cluster dependencies.
3. **Potential issues** – In scRNA‑seq, early clusters may be noisy due to sparsity; using them to guide masking could propagate errors.  The method also requires additional computation to update cluster assignments.
4. **Mitigation** – Employ robust clustering (e.g., graph‑based) to identify stable gene clusters.  Gradually increase reliance on guided masks as the model stabilises.  Use a hybrid of random and guided masking to avoid overfitting to early clusters.
5. **Usability decision** – With careful implementation, self‑guided masking could make scMAE more efficient and capture biologically meaningful relationships.  It is considered usable.

Citations:【944398912377764†L14-L23】

### 048. Masked Autoencoders Are Scalable Vision Learners

- **排序文件**：`01_PDF论文_按推荐程度排序/048_可用_Masked_Autoencoders_Are_Scalable_Vision_Learners.pdf`
- **推荐等级**：可用（score=70）
- **SHA-256**：`1B490443925C72A2B7C770F90DD797E248729AE34A57E1ABFE9ED36751C4CC5B`
- **去重说明**：合并 2 条来源记录；包含 2 个PDF哈希。
- **来源记录**：文库2#3 03_MAE_general.pdf；文库3#3 03_MAE_He_2022_CVPR.pdf

#### 来源说明：文库2 #3

说明来源：文库2/scmae_improvement_reports/03_MaskedAutoencoders_general.txt
- **标题**：Masked Autoencoders Are Scalable Vision Learners
- **年份/会议**：2022 / CVPR
- **备注**：CVF开放PDF。

Masked Autoencoders (MAE) – general framework

1. **Efficacy & principle** – MAEs use an asymmetric encoder–decoder architecture, where a subset of input tokens is visible to the encoder while a lightweight decoder reconstructs the full input.  High masking ratios (e.g., 75%) reduce redundancy and encourage holistic learning, resulting in efficient training and improved accuracy【157010519668395†L14-L24】【157010519668395†L114-L129】.
2. **Integration strategy** – For scRNA‑seq, the encoder can process unmasked genes, while the decoder reconstructs masked gene expression.  Using a high masking ratio encourages the model to capture global gene dependencies and reduces computational cost.
3. **Potential issues** – Standard MAEs assume continuous patches with spatial locality; gene expression data have no spatial structure.  High masking may remove too much information, hindering learning if gene correlations are weak.
4. **Mitigation** – Employ domain‑specific masking patterns (e.g., mask genes within the same pathway) and adjust masking ratio based on sparsity.  Use zero‑inflated negative binomial (ZINB) reconstruction loss to model count data instead of mean‑squared error.
5. **Usability decision** – MAE provides a powerful template, but domain adaptation is necessary.  It is usable when combined with scRNA‑specific loss functions and masking strategies.

Citations:【157010519668395†L14-L24】【157010519668395†L114-L129】

#### 来源说明：文库3 #3

说明来源：文库3/manifest_30_papers.csv
- **标题**：Masked Autoencoders Are Scalable Vision Learners
- **年份/会议**：CVPR 2022
- **方向/领域**：Masked modeling
- **与 scMAE 改良相关性**：High mask ratio/asymmetric encoder-decoder design inspires efficient scMAE pretraining.

### 049. Adaptive-Masking Policy with Deep Reinforcement Learning for Self-Supervised Medical Image Segmentation

- **排序文件**：`01_PDF论文_按推荐程度排序/049_可用_Adaptive-Masking_Policy_with_Deep_Reinforcement_Learning_for_Self-Supervised_Medical_Image_Segmentation.pdf`
- **推荐等级**：可用（score=70）
- **SHA-256**：`7F8B30DD0816E863F8F41C604B54A38A8AB1C460D0EDF260599055BA72696C26`
- **来源记录**：文库2#4 04_RL_Adaptive_Masking.pdf

#### 来源说明：文库2 #4

说明来源：文库2/scmae_improvement_reports/04_RL_AdaptiveMasking.txt
- **标题**：Adaptive-Masking Policy with Deep Reinforcement Learning for Self-Supervised Medical Image Segmentation
- **年份/会议**：2023 / ICME
- **备注**：作者主页开放PDF。

Adaptive masking via reinforcement learning

1. **Efficacy & principle** – An ICME 2023 paper frames the masking process as a Markov decision process and uses a dueling deep Q‑network (DQN) to adaptively select mask positions for self‑supervised medical image segmentation.  The agent receives feedback from the reconstruction network and learns to mask informative regions, improving reconstruction and downstream performance【615655185214976†L22-L33】【615655185214976†L86-L124】.
2. **Integration strategy** – Apply a similar RL framework to scMAE by defining states as the current mask and network performance and actions as selecting genes to mask.  Reward signals can be based on reconstruction loss or clustering quality.  The agent learns which genes contribute most to representation learning.
3. **Potential issues** – RL can be computationally intensive and unstable.  In scRNA‑seq, reward signals may be noisy due to technical variability.  The search space of gene masks is enormous, which complicates learning.
4. **Mitigation** – Reduce the action space by grouping genes (e.g., pathways) and mask groups rather than individual genes.  Use stable RL algorithms such as proximal policy optimization (PPO) and incorporate prior knowledge to guide exploration.
5. **Usability decision** – Adaptive masking has potential to focus on informative genes and improve scMAE.  With careful design and grouping, it is usable.

Citations:【615655185214976†L22-L33】【615655185214976†L86-L124】

### 050. The Dynamic Duo of Collaborative Masking and Target for Advanced Masked Autoencoder Learning

- **排序文件**：`01_PDF论文_按推荐程度排序/050_可用_The_Dynamic_Duo_of_Collaborative_Masking_and_Target_for_Advanced_Masked_Autoencoder_Learning.pdf`
- **推荐等级**：可用（score=70）
- **SHA-256**：`CC5B5B749A6D47014072E5844F086FB6600C9FBE0E26D540D143A2F701341E39`
- **来源记录**：文库2#6 06_CMT-MAE.pdf

#### 来源说明：文库2 #6

说明来源：文库2/scmae_improvement_reports/06_CMT_MAE.txt
- **标题**：The Dynamic Duo of Collaborative Masking and Target for Advanced Masked Autoencoder Learning
- **年份/会议**：2025 / AAAI
- **备注**：AAAI开放PDF。

CMT‑MAE – Collaborative Masked Autoencoder

1. **Efficacy & principle** – The AAAI 2025 paper proposes a teacher–student collaborative masking strategy.  Attention maps from a stronger teacher guide the student’s masking; collaborative targets aggregate the teacher’s attention maps to provide informative supervision, improving MAE pre‑training【813482395155896†L14-L25】.
2. **Integration strategy** – For scMAE, a pre‑trained model (teacher) could generate importance scores for genes.  A student model uses these scores to choose which genes to mask and to distil knowledge from the teacher via collaborative targets.  This may accelerate training and improve gene representation.
3. **Potential issues** – Requires a pre‑trained teacher; if the teacher is biased, it may reinforce wrong patterns.  Collaborative training may be complex and sensitive to hyperparameters.
4. **Mitigation** – Use an ensemble of teachers (e.g., different initializations or architectures) to reduce bias.  Apply regularization to prevent overfitting to teacher signals.
5. **Usability decision** – With a reliable teacher, collaborative masking can guide scMAE to focus on important genes.  It is usable but demands additional resources.

Citations:【813482395155896†L14-L25】

### 051. Reinforcement Learning meets Masked Video Modeling: Trajectory-Guided Adaptive Token Selection

- **排序文件**：`01_PDF论文_按推荐程度排序/051_可用_Reinforcement_Learning_meets_Masked_Video_Modeling_Trajectory-Guided_Adaptive_Token_Selection.pdf`
- **推荐等级**：可用（score=70）
- **SHA-256**：`6FA09AEC61F20AC98AF4D46E60CEA767574198FA5B9C8FB4A92492EC91954E9F`
- **来源记录**：文库2#7 07_TATS.pdf

#### 来源说明：文库2 #7

说明来源：文库2/scmae_improvement_reports/07_TATS.txt
- **标题**：Reinforcement Learning meets Masked Video Modeling: Trajectory-Guided Adaptive Token Selection
- **年份/会议**：2025 / arXiv
- **备注**：arXiv开放PDF。

TATS – Trajectory‑Aware Adaptive Token Sampler

1. **Efficacy & principle** – TATS is an adaptive token selection method for masked video modeling.  It uses proximal policy optimization (PPO) to learn a token sampler that focuses on motion‑centric tokens, enabling higher masking ratios without performance loss【189159888125753†L41-L53】.
2. **Integration strategy** – For scMAE, adapt TATS by treating gene expression patterns across cells as a trajectory; the token sampler learns to select genes whose expression changes inform the latent structure.  PPO can train the sampler to pick dynamic genes or rare cell‑type markers.
3. **Potential issues** – Genes do not possess temporal dynamics like video frames; mapping them to trajectories may be unnatural.  The RL policy might overemphasize fluctuating genes and ignore stable but informative ones.
4. **Mitigation** – Define trajectories based on developmental or pseudotime ordering of cells.  Combine the adaptive sampler with static gene importance measures to ensure stability.
5. **Usability decision** – If pseudotime ordering is available, TATS can focus on dynamic genes.  It may be usable in developmental datasets but less general.

Citations:【189159888125753†L41-L53】

### 052. Anatomically-guided Masked Autoencoder with Domain-Adaptive Prompting for multimodal cerebral aneurysm detection and segmentation

- **排序文件**：`01_PDF论文_按推荐程度排序/052_可用_Anatomically-guided_Masked_Autoencoder_with_Domain-Adaptive_Prompting_for_multimodal_cerebral_aneurysm_detection_and_segmentation.pdf`
- **推荐等级**：可用（score=70）
- **SHA-256**：`8307477DB74027393CB73D6DDD4CAE5D7EC7A999C446E30EF119FE5FF8C1B0AA`
- **来源记录**：文库2#8 08_AMAP.pdf

#### 来源说明：文库2 #8

说明来源：文库2/scmae_improvement_reports/08_AMAP.txt
- **标题**：Anatomically-guided Masked Autoencoder with Domain-Adaptive Prompting for multimodal cerebral aneurysm detection and segmentation
- **年份/会议**：2026 / npj Digital Medicine
- **备注**：Nature开放PDF。

AMAP – Anatomically‑Guided Masked Autoencoder with Domain‑Adaptive Prompting

1. **Efficacy & principle** – AMAP applies an anatomy‑guided MAE pre‑training strategy for cerebral aneurysm detection.  It uses domain‑adaptive prompts and boundary‑aware contrastive learning to improve segmentation accuracy and reduce false positives【896952137263954†L90-L101】.
2. **Integration strategy** – For scMAE, anatomical guidance translates to biological guidance.  Provide domain prompts based on cell type labels or gene ontology to the encoder and train with contrastive objectives that emphasize boundaries between cell types.  This could help the model learn discriminative features.
3. **Potential issues** – Domain prompts require prior knowledge and labelled data.  Overfitting to prompts might reduce the model’s ability to discover novel cell types.
4. **Mitigation** – Use semi‑supervised learning: apply domain prompts for a subset of cells while maintaining unsupervised learning on the rest.  Regularize the model to prevent overemphasis on prompts.
5. **Usability decision** – In datasets with partial labels or well‑defined cell types, AMAP‑like prompting can enhance scMAE.  It is usable with caution.

Citations:【896952137263954†L90-L101】

### 053. DAP-MAE: Domain-Adaptive Point Cloud Masked Autoencoder for Effective Cross-Domain Learning

- **排序文件**：`01_PDF论文_按推荐程度排序/053_可用_DAP-MAE_Domain-Adaptive_Point_Cloud_Masked_Autoencoder_for_Effective_Cross-Domain_Learning.pdf`
- **推荐等级**：可用（score=70）
- **SHA-256**：`34B6BA6B5F4D6047270826FD8AE6F81FA29CE6DD696A05FB51358828E760F511`
- **来源记录**：文库2#9 09_DAP-MAE.pdf

#### 来源说明：文库2 #9

说明来源：文库2/scmae_improvement_reports/09_DAP_MAE.txt
- **标题**：DAP-MAE: Domain-Adaptive Point Cloud Masked Autoencoder for Effective Cross-Domain Learning
- **年份/会议**：2025 / arXiv
- **备注**：arXiv开放PDF。

DAP‑MAE – Domain‑Adaptive Point‑Cloud Masked Autoencoder

1. **Efficacy & principle** – DAP‑MAE introduces a heterogeneous domain adapter and a domain feature generator to integrate cross‑domain point‑cloud datasets in masked autoencoder pretraining.  It improves cross‑domain performance by adaptively combining features across domains【449733547596744†L41-L52】.
2. **Integration strategy** – For scMAE, consider multi‑omics or cross‑species data as different domains.  Design a domain adapter network that aligns gene distributions across domains before masked reconstruction.  Use a domain feature generator to augment the training data and encourage the model to learn domain‑invariant representations.
3. **Potential issues** – Domain adaptation may blur biologically meaningful differences between species or modalities.  Aligning distributions might remove subtle signals necessary for cell‑type identification.
4. **Mitigation** – Apply partial alignment that preserves key biological pathways.  Use adversarial training to learn domain‑invariant features while enforcing classification performance on domain‑specific tasks.
5. **Usability decision** – DAP‑MAE is promising for integrating multiple datasets; with careful domain alignment, it is usable.

Citations:【449733547596744†L41-L52】

### 054. An Empirical Study of Multiple Masking in Masked Autoencoder / Conditional MAE

- **排序文件**：`01_PDF论文_按推荐程度排序/054_可用_An_Empirical_Study_of_Multiple_Masking_in_Masked_Autoencoder_Conditional_MAE.pdf`
- **推荐等级**：可用（score=70）
- **SHA-256**：`F2A1AECBC6BAD2F50A4E1F4E5BA415AB71884AA2906318309AE5B14049CE041C`
- **来源记录**：文库2#11 11_Conditional_MAE.pdf

#### 来源说明：文库2 #11

说明来源：文库2/scmae_improvement_reports/11_ConditionalMAE.txt
- **标题**：An Empirical Study of Multiple Masking in Masked Autoencoder / Conditional MAE
- **年份/会议**：2025 / OpenReview/withdrawn
- **备注**：OpenReview PDF可下载，但状态为withdrawn，不建议作为核心引用。

Conditional MAE – Multiple Masking Stages

1. **Efficacy & principle** – The ICLR 2025 submission introduces a conditional MAE where subsequent masking stages depend on the representations learned from previous stages.  This sequential masking improves optimisation and representation quality【202810644948851†L24-L42】.
2. **Integration strategy** – Apply conditional masking to scMAE by training the model in stages: initial training with random masking, followed by further training where masks are conditioned on reconstructed gene correlations.  This encourages the model to refine its representations iteratively.
3. **Potential issues** – Sequential stages increase training time.  Conditioning masks on model outputs may reinforce biases if early representations are poor.
4. **Mitigation** – Use a validation set to select when to switch stages.  Combine conditional masking with random masking to maintain diversity.
5. **Usability decision** – Sequential masking can refine scMAE representations, but complexity must be managed.  It is usable with monitoring.

Citations:【202810644948851†L24-L42】

### 055. Multi-Facet Clustering Variational Autoencoders

- **排序文件**：`01_PDF论文_按推荐程度排序/055_可用_Multi-Facet_Clustering_Variational_Autoencoders.pdf`
- **推荐等级**：可用（score=70）
- **SHA-256**：`69B9A90276F9F4D6DF4D7B53FB378479B5F5416D329CF80D556E078B9D8FB649`
- **来源记录**：文库2#12 12_MFCVAE.pdf

#### 来源说明：文库2 #12

说明来源：文库2/scmae_improvement_reports/12_MFCVAE.txt
- **标题**：Multi-Facet Clustering Variational Autoencoders
- **年份/会议**：2021 / NeurIPS
- **备注**：NeurIPS开放PDF。

MFCVAE – Multi‑Facet Clustering Variational Autoencoder

1. **Efficacy & principle** – MFCVAE learns multiple clusterings simultaneously by introducing a hierarchy of latent variables with Mixture‑of‑Gaussians priors.  It disentangles different abstract characteristics and provides compositional latent space【930723102025873†L10-L32】.
2. **Integration strategy** – For scMAE, use a multi‑facet latent space to capture different biological aspects (e.g., cell cycle, differentiation stage, cell type).  Each facet could have its own masking predictor, enabling the model to learn multiple clustering views.
3. **Potential issues** – Determining the number of facets and their interpretation can be challenging.  The model may overfit or produce redundant facets if data do not support the complexity.
4. **Mitigation** – Apply Bayesian non‑parametric approaches to infer the number of facets.  Use domain knowledge to guide facet interpretation and combine facets for final clustering.
5. **Usability decision** – Multi‑facet representation could enrich scMAE but requires careful design.  It is usable for complex datasets with multiple biological factors.

Citations:【930723102025873†L10-L32】

### 056. ScInfoVAE: interpretable dimensional reduction of single cell transcription data with variational autoencoders and extended mutual information regularization

- **排序文件**：`01_PDF论文_按推荐程度排序/056_可用_ScInfoVAE_interpretable_dimensional_reduction_of_single_cell_transcription_data_with_variational_autoencoders_and_extended_mutual.pdf`
- **推荐等级**：可用（score=70）
- **SHA-256**：`66FC6BA60255DAFFFDD3B61EB42A2A35FF2F53259C1BAD9CE48F160BE40179C6`
- **来源记录**：文库2#14 14_ScInfoVAE.pdf

#### 来源说明：文库2 #14

说明来源：文库2/scmae_improvement_reports/14_scInfoVAE.txt
- **标题**：ScInfoVAE: interpretable dimensional reduction of single cell transcription data with variational autoencoders and extended mutual information regularization
- **年份/会议**：2023 / BioData Mining
- **备注**：Springer开放PDF。

ScInfoVAE – Mutual‑Information Variational Autoencoder for scRNA

1. **Efficacy & principle** – ScInfoVAE is based on InfoVAE and extends mutual information regularisation.  It combines InfoVAE with a zero‑inflated negative binomial model to reconstruct noisy scRNA‑seq data and learn efficient low‑dimensional representations.  Experiments on 15 real datasets show high clustering performance and interpretability【81971360896341†L91-L110】.
2. **Integration strategy** – Incorporate mutual information regularisation into scMAE to encourage the latent representation to retain more information about the input genes.  Combine with ZINB decoding to handle dropout.
3. **Potential issues** – Mutual information regularisation may require careful balancing; too strong regularisation may overfit noise.  Implementation complexity increases.
4. **Mitigation** – Use a hyperparameter search to balance mutual information and reconstruction losses.  Employ dropout and early stopping to avoid overfitting.
5. **Usability decision** – Extending scMAE with InfoVAE concepts could enhance interpretability and robustness; it is usable with tuning.

Citations:【81971360896341†L91-L110】

### 057. scVAEDer: integrating deep diffusion models and variational autoencoders for single-cell transcriptomics analysis

- **排序文件**：`01_PDF论文_按推荐程度排序/057_可用_scVAEDer_integrating_deep_diffusion_models_and_variational_autoencoders_for_single-cell_transcriptomics_analysis.pdf`
- **推荐等级**：可用（score=70）
- **SHA-256**：`934978A9EC16F070DCA6D75F729E59C8E7DFDBA5C5A0E393E92E7897C5D27CDA`
- **来源记录**：文库2#15 15_scVAEDer.pdf

#### 来源说明：文库2 #15

说明来源：文库2/scmae_improvement_reports/15_scVAEDer.txt
- **标题**：scVAEDer: integrating deep diffusion models and variational autoencoders for single-cell transcriptomics analysis
- **年份/会议**：2025 / Genome Biology
- **备注**：Springer开放PDF。

scVAEDer – Variational Autoencoder and Deep Diffusion Model

1. **Efficacy & principle** – scVAEDer combines variational autoencoders with deep diffusion models to learn robust embeddings and generate novel scRNA‑seq data.  It enables perturbation prediction and identification of master regulators【616742116960631†L69-L79】.
2. **Integration strategy** – Use a diffusion‑based prior in scMAE: after encoding visible genes, sample latent variables via diffusion and decode with a variational decoder.  This may improve generative quality and capture complex gene dependencies.
3. **Potential issues** – Diffusion models are computationally expensive and require careful denoising schedules.  Combining them with MAE may complicate training and increase runtime.
4. **Mitigation** – Use latent diffusion (scLDM) to reduce dimension and speed up sampling.  Train diffusion separately and fine‑tune within scMAE.
5. **Usability decision** – Diffusion‑VAE integration can enhance generative capability but at high computational cost.  It is usable for tasks requiring generation.

Citations:【616742116960631†L69-L79】

### 058. Gene selection for single cell RNA-seq data via fuzzy rough iterative computation model

- **排序文件**：`01_PDF论文_按推荐程度排序/058_可用_Gene_selection_for_single_cell_RNA-seq_data_via_fuzzy_rough_iterative_computation_model.pdf`
- **推荐等级**：可用（score=70）
- **SHA-256**：`2084D463F7133CCBB56A68364582379DBF942F5B54E66F30FF61904D5B044A5B`
- **来源记录**：文库2#16 16_FRIC.pdf

#### 来源说明：文库2 #16

说明来源：文库2/scmae_improvement_reports/16_FRIC.txt
- **标题**：Gene selection for single cell RNA-seq data via fuzzy rough iterative computation model
- **年份/会议**：2025 / Artificial Intelligence Review
- **备注**：Springer开放PDF。

FRIC – Fuzzy Rough Set Feature Selection

1. **Efficacy & principle** – The FRIC model applies fuzzy rough sets to gene selection.  It replaces equality of gene values with a distance‑based fuzzy symmetric relation and iteratively computes dependency functions to select genes, improving scRNA‑seq clustering【975523413296288†L67-L85】.
2. **Integration strategy** – Use FRIC to pre‑select informative genes before training scMAE.  This reduces dimensionality and emphasises genes with discriminative power, potentially speeding up training and improving clustering.
3. **Potential issues** – Fuzzy rough set methods may depend on threshold parameters and can be sensitive to noise.  They may remove genes that are important for novel cell types.
4. **Mitigation** – Tune thresholds using cross‑validation and incorporate biological knowledge to retain marker genes.  Combine FRIC with unsupervised gene ranking to mitigate bias.
5. **Usability decision** – As a preprocessing step, FRIC can enhance scMAE by reducing noise and focusing on informative genes; it is usable.

Citations:【975523413296288†L67-L85】

### 059. Soft Graph Clustering for single-cell RNA Sequencing Data

- **排序文件**：`01_PDF论文_按推荐程度排序/059_可用_Soft_Graph_Clustering_for_single-cell_RNA_Sequencing_Data.pdf`
- **推荐等级**：可用（score=70）
- **SHA-256**：`E4AFDDD6B09332EDF688DB71C5552CF22DC3F321253D286A2F7C10BF1FEB8A7A`
- **来源记录**：文库2#17 17_scSGC.pdf

#### 来源说明：文库2 #17

说明来源：文库2/scmae_improvement_reports/17_scSGC.txt
- **标题**：Soft Graph Clustering for single-cell RNA Sequencing Data
- **年份/会议**：2025 / arXiv/BMC Bioinformatics
- **备注**：arXiv开放PDF。

scSGC – Soft Graph Clustering

1. **Efficacy & principle** – scSGC addresses limitations of hard graph construction by using a zero‑inflated negative binomial autoencoder and dual‑channel cut‑informed soft graph embedding.  It formulates clustering as an optimal transport problem, providing smoother embeddings and better clustering performance【763472327302984†L49-L70】.
2. **Integration strategy** – Replace scMAE’s fixed graph or similarity matrix with a soft graph embedding module.  Use soft edges and ZINB autoencoder to construct a continuous representation of cell relationships.  This can serve as the input for the masked autoencoder.
3. **Potential issues** – Optimal transport computations can be expensive.  Soft graphs may blur sharp boundaries between cell types.
4. **Mitigation** – Employ approximate optimal transport and regularise soft edges to avoid over‑smoothing.  Combine with cluster‑aware losses to maintain distinct cell types.
5. **Usability decision** – Soft graph clustering provides more flexible cell relationships; it is usable as a pre‑processing step for scMAE.

Citations:【763472327302984†L49-L70】

### 060. q-Diffusion leverages the full dimensionality of gene coexpression in single-cell transcriptomics

- **排序文件**：`01_PDF论文_按推荐程度排序/060_可用_q-Diffusion_leverages_the_full_dimensionality_of_gene_coexpression_in_single-cell_transcriptomics.pdf`
- **推荐等级**：可用（score=70）
- **SHA-256**：`AD4F8B30AC6675B4775551ACC200E9F44DF9AAFEA4E5D8A49C667530FC6B5CCF`
- **来源记录**：文库2#18 18_q-Diffusion.pdf

#### 来源说明：文库2 #18

说明来源：文库2/scmae_improvement_reports/18_q_diffusion.txt
- **标题**：q-Diffusion leverages the full dimensionality of gene coexpression in single-cell transcriptomics
- **年份/会议**：2024 / Communications Biology
- **备注**：Nature开放PDF。

q‑Diffusion for scRNA‑seq

1. **Efficacy & principle** – q‑Diffusion introduces a q‑diffused kernel that captures the coexpression structure of an entire gene library, enabling the analysis of high‑dimensional gene interactions.  It improves differential effect detection, unsupervised clustering and spatial segmentation in case studies【675472413082403†L84-L97】.
2. **Integration strategy** – Incorporate the q‑diffused kernel as a similarity measure in scMAE.  Use it to construct neighbourhoods or to weight reconstruction loss based on transcriptional proximity, thus better capturing coexpression patterns.
3. **Potential issues** – Computation of q‑diffusion may be expensive for large datasets.  Parameter tuning (e.g., choice of q) influences performance.
4. **Mitigation** – Precompute the kernel on downsampled data or use random feature approximations.  Conduct parameter search to determine optimal q values.
5. **Usability decision** – q‑Diffusion offers a powerful way to leverage full gene coexpression; with computational optimisation, it is usable.

Citations:【675472413082403†L84-L97】

### 061. scDiffusion: conditional generation of high-quality single-cell data using diffusion model

- **排序文件**：`01_PDF论文_按推荐程度排序/061_可用_scDiffusion_conditional_generation_of_high-quality_single-cell_data_using_diffusion_model.pdf`
- **推荐等级**：可用（score=70）
- **SHA-256**：`ED3CE8D2D76CEED1B2239363F14421922893CF26625698B703226FF5CE73E2B9`
- **来源记录**：文库2#19 19_scDiffusion.pdf

#### 来源说明：文库2 #19

说明来源：文库2/scmae_improvement_reports/19_scDiffusion.txt
- **标题**：scDiffusion: conditional generation of high-quality single-cell data using diffusion model
- **年份/会议**：2024 / Bioinformatics/arXiv
- **备注**：OUP正式文章页；arXiv开放PDF用于下载。

scDiffusion – Conditional Diffusion Model for scRNA‑seq

1. **Efficacy & principle** – scDiffusion is a generative model that combines a diffusion model with multiple classifiers.  It introduces gradient interpolation to generate continuous cell developmental trajectories and can produce rare or out‑of‑distribution cell types【566158587869794†L139-L167】.
2. **Integration strategy** – Use scDiffusion as an augmentation tool for scMAE: generate synthetic cells to balance classes or to train the autoencoder on more diverse examples.  Alternatively, integrate diffusion prior into scMAE to guide reconstruction toward realistic gene patterns.
3. **Potential issues** – Diffusion models are computationally heavy.  Synthetic cells may introduce noise or unrealistic patterns if not properly controlled.
4. **Mitigation** – Use classifier guidance to steer diffusion toward biologically plausible gene profiles.  Validate synthetic cells using known markers.
5. **Usability decision** – scDiffusion is valuable for data augmentation and generative modelling; it is usable for augmenting scMAE training.

Citations:【566158587869794†L139-L167】

### 062. Scalable Single-Cell Gene Expression Generation with Latent Diffusion Models

- **排序文件**：`01_PDF论文_按推荐程度排序/062_可用_Scalable_Single-Cell_Gene_Expression_Generation_with_Latent_Diffusion_Models.pdf`
- **推荐等级**：可用（score=70）
- **SHA-256**：`9AE1F1FBA91C369F942AA7E57662B94BB9401AF6B735A7F4B11BF5DE6E9F4BDB`
- **来源记录**：文库2#20 20_scLDM.pdf

#### 来源说明：文库2 #20

说明来源：文库2/scmae_improvement_reports/20_scLDM.txt
- **标题**：Scalable Single-Cell Gene Expression Generation with Latent Diffusion Models
- **年份/会议**：2025/2026 / arXiv
- **备注**：arXiv开放PDF。

scLDM – Latent Diffusion Model for scRNA‑seq

1. **Efficacy & principle** – scLDM uses a permutation‑equivariant transformer variational autoencoder with cross‑attention blocks and a latent diffusion prior.  It employs classifier‑free guidance to allow multi‑conditional generation and achieves superior performance on observational and perturbational data【139365790920974†L51-L66】.
2. **Integration strategy** – Use a latent diffusion prior within scMAE: after encoding visible genes, apply diffusion in latent space to generate more realistic reconstructions.  Cross‑attention blocks can model gene–cell interactions within scMAE.
3. **Potential issues** – Diffusion models in latent space still require careful tuning and may slow training.  The complexity of cross‑attention may increase computational cost.
4. **Mitigation** – Use smaller latent dimensions and efficient transformer architectures.  Pretrain the diffusion prior separately and fine‑tune within scMAE.
5. **Usability decision** – scLDM offers strong generative capabilities; with efficiency considerations, it is usable.

Citations:【139365790920974†L51-L66】

### 063. Attention-Guided Probabilistic Diffusion Model for Generating Cell-Type-Specific Gene Regulatory Networks from Gene Expression Profiles

- **排序文件**：`01_PDF论文_按推荐程度排序/063_可用_Attention-Guided_Probabilistic_Diffusion_Model_for_Generating_Cell-Type-Specific_Gene_Regulatory_Networks_from_Gene_Expression_Pro.pdf`
- **推荐等级**：可用（score=70）
- **SHA-256**：`C40866DA25CB205FCCFA6255AF61E3A6026C3D01B4B9ACB798488B2FB49B7596`
- **来源记录**：文库2#21 21_Planet.pdf

#### 来源说明：文库2 #21

说明来源：文库2/scmae_improvement_reports/21_Planet.txt
- **标题**：Attention-Guided Probabilistic Diffusion Model for Generating Cell-Type-Specific Gene Regulatory Networks from Gene Expression Profiles
- **年份/会议**：2025 / Genes
- **备注**：MDPI开放PDF。

Planet – Attention‑Guided Probabilistic Diffusion Model

1. **Efficacy & principle** – Planet uses a triple hybrid attention transformer and diffusion modelling to generate cell‑specific gene regulatory networks with fast sampling and improved global consistency【910468027404436†L174-L199】.
2. **Integration strategy** – Use Planet’s attention mechanisms within scMAE to focus on regulatory gene interactions.  Incorporate diffusion sampling to model gene network generation, potentially improving downstream tasks like gene network inference.
3. **Potential issues** – The method is specifically designed for regulatory network generation, not clustering.  Integrating it may overcomplicate scMAE.
4. **Mitigation** – Use Planet as a post‑processing module: after scMAE embedding, apply Planet to generate regulatory networks for each cluster.  Keep scMAE core simple.
5. **Usability decision** – Planet can be leveraged for downstream tasks but may not directly improve clustering.  Its attention mechanism is usable for network inference.

Citations:【910468027404436†L174-L199】

### 064. GRouNdGAN: GRN-guided simulation of single-cell RNA-seq data using causal generative adversarial networks

- **排序文件**：`01_PDF论文_按推荐程度排序/064_可用_GRouNdGAN_GRN-guided_simulation_of_single-cell_RNA-seq_data_using_causal_generative_adversarial_networks.pdf`
- **推荐等级**：可用（score=70）
- **SHA-256**：`5D65724FE66074AC244AF769CB04B310F1428890C756720E23ECE00DBABA253C`
- **来源记录**：文库2#22 22_GRouNdGAN.pdf

#### 来源说明：文库2 #22

说明来源：文库2/scmae_improvement_reports/22_GRouNdGAN.txt
- **标题**：GRouNdGAN: GRN-guided simulation of single-cell RNA-seq data using causal generative adversarial networks
- **年份/会议**：2024 / Nature Communications
- **备注**：Nature开放PDF。

GRouNdGAN – Gene Regulatory Network‑guided Generative Model

1. **Efficacy & principle** – GRouNdGAN imposes a user‑defined gene regulatory network on a generative adversarial network, enabling simulation of both steady and transient single‑cell data.  It captures transcription factor–gene dependencies and preserves trajectories【580321926828693†L73-L87】.
2. **Integration strategy** – Use GRouNdGAN to generate synthetic training data consistent with known regulatory networks.  Train scMAE on this data to learn biologically realistic gene dependencies.  Alternatively, incorporate a GRN prior into scMAE’s decoder.
3. **Potential issues** – The method relies on accurate gene regulatory networks.  If the GRN is incomplete or incorrect, it may misguide generation.  GANs can be unstable.
4. **Mitigation** – Use ensemble GRNs or integrate multiple sources of network information.  Regularize the GAN and monitor training to prevent mode collapse.
5. **Usability decision** – GRouNdGAN is valuable when reliable GRNs are available; it can augment scMAE training.  It is usable with caution.

Citations:【580321926828693†L73-L87】

### 065. Realistic in silico generation and augmentation of single-cell RNA-seq data using generative adversarial networks

- **排序文件**：`01_PDF论文_按推荐程度排序/065_可用_Realistic_in_silico_generation_and_augmentation_of_single-cell_RNA-seq_data_using_generative_adversarial_networks.pdf`
- **推荐等级**：可用（score=70）
- **SHA-256**：`84686624F66CFABA89077C05B183B7D0A8581C2501B118FFE1C6FC8A3ECE5BCA`
- **来源记录**：文库2#23 23_cscGAN.pdf

#### 来源说明：文库2 #23

说明来源：文库2/scmae_improvement_reports/23_cscGAN.txt
- **标题**：Realistic in silico generation and augmentation of single-cell RNA-seq data using generative adversarial networks
- **年份/会议**：2020 / Nature Communications
- **备注**：Nature开放PDF。

cscGAN – Conditional Single‑cell GAN

1. **Efficacy & principle** – cscGAN learns nonlinear gene dependencies and generates realistic single‑cell RNA‑seq data.  It produces cells of defined types and improves marker detection and classifier robustness【781356676628973†L75-L91】.
2. **Integration strategy** – Use cscGAN for data augmentation: generate additional cells for underrepresented types and train scMAE on a balanced dataset.  Alternatively, incorporate the generator into scMAE’s training pipeline to refine latent representations.
3. **Potential issues** – GAN training is notoriously unstable.  Synthetic cells might not capture rare cell states or subtle gene interactions.
4. **Mitigation** – Employ conditional GAN with label smoothing and spectral normalization to stabilize training.  Validate synthetic cells via downstream analyses.
5. **Usability decision** – cscGAN is useful for augmenting training data and enhancing scMAE’s robustness; it is usable.

Citations:【781356676628973†L75-L91】

### 066. scVGATAE: A Variational Graph Attentional Autoencoder Model for Clustering Single-Cell RNA-seq Data

- **排序文件**：`01_PDF论文_按推荐程度排序/066_可用_scVGATAE_A_Variational_Graph_Attentional_Autoencoder_Model_for_Clustering_Single-Cell_RNA-seq_Data.pdf`
- **推荐等级**：可用（score=70）
- **SHA-256**：`B5C2442C5EC43FA96DAF46431F0E9B7ECF007B5AAC2DB61738E660481AD95480`
- **来源记录**：文库2#24 24_scVGATAE.pdf

#### 来源说明：文库2 #24

说明来源：文库2/scmae_improvement_reports/24_scVGATAE.txt
- **标题**：scVGATAE: A Variational Graph Attentional Autoencoder Model for Clustering Single-Cell RNA-seq Data
- **年份/会议**：2024 / Biology
- **备注**：MDPI开放PDF。

scVGATAE – Variational Graph Attention Autoencoder

1. **Efficacy & principle** – scVGATAE constructs a reliable cell graph via network denoising and uses a variational graph attention autoencoder to learn low‑dimensional cell embeddings.  It achieves superior clustering performance on scRNA‑seq datasets【568304791288546†L550-L558】.
2. **Integration strategy** – Replace scMAE’s encoder with a graph attention variational autoencoder to better capture cell–cell relationships before masking.  Use the learned embeddings as inputs to the masked reconstruction decoder.
3. **Potential issues** – Graph construction relies on initial similarity measures, which may be noisy.  Variational graph models require careful tuning.
4. **Mitigation** – Use adaptive graph construction and robust similarity metrics.  Pretrain the variational graph encoder separately before integrating into scMAE.
5. **Usability decision** – The method offers strong embedding capabilities; integrating its graph attention mechanisms into scMAE is promising and usable.

Citations:【568304791288546†L550-L558】

### 067. scSiameseClu: A Siamese Clustering Framework for Interpreting Single-cell RNA Sequencing Data

- **排序文件**：`01_PDF论文_按推荐程度排序/067_可用_scSiameseClu_A_Siamese_Clustering_Framework_for_Interpreting_Single-cell_RNA_Sequencing_Data.pdf`
- **推荐等级**：可用（score=70）
- **SHA-256**：`016C6AF0A47E3306D8331833D3A3B68B85A08EF885743953A5B1F288D985ABCE`
- **来源记录**：文库2#28 28_scSiameseClu.pdf

#### 来源说明：文库2 #28

说明来源：文库2/scmae_improvement_reports/28_scSiameseClu.txt
- **标题**：scSiameseClu: A Siamese Clustering Framework for Interpreting Single-cell RNA Sequencing Data
- **年份/会议**：2025 / IJCAI
- **备注**：IJCAI开放PDF。

scSiameseClu – Siamese Clustering with Dual Augmentation

1. **Efficacy & principle** – scSiameseClu employs a dual augmentation module (biologically informed gene perturbations and graph perturbations), a Siamese fusion module for cross‑correlation refinement, and optimal transport clustering.  It mitigates over‑smoothing and captures complex intercellular information【609583138125184†L24-L35】【609583138125184†L122-L130】.
2. **Integration strategy** – Use two augmented views of scRNA‑seq data within scMAE: one with masked genes and one with perturbed cell graphs.  Train a Siamese network to maximise agreement between views, then apply masked autoencoding.  Use optimal transport to align cluster assignments.
3. **Potential issues** – Optimal transport clustering is computationally expensive.  Dual augmentations require careful design to avoid destroying biological signals.
4. **Mitigation** – Use approximate optimal transport and restrict perturbations to biologically plausible ranges.  Perform sensitivity analysis to select augmentation intensity.
5. **Usability decision** – scSiameseClu provides robust augmentation and clustering techniques; integrating its Siamese fusion into scMAE is usable.

Citations:【609583138125184†L24-L35】【609583138125184†L122-L130】

### 068. Deep single-cell RNA-seq data clustering with graph prototypical contrastive learning

- **排序文件**：`01_PDF论文_按推荐程度排序/068_可用_Deep_single-cell_RNA-seq_data_clustering_with_graph_prototypical_contrastive_learning.pdf`
- **推荐等级**：可用（score=70）
- **SHA-256**：`C2237828A552D9F3EDA13AE1036D91DA07795DF272E84691DFC894A057B32165`
- **来源记录**：文库2#29 29_scGPCL.pdf

#### 来源说明：文库2 #29

说明来源：文库2/scmae_improvement_reports/29_scGPCL.txt
- **标题**：Deep single-cell RNA-seq data clustering with graph prototypical contrastive learning
- **年份/会议**：2023 / Bioinformatics / ICML CompBio
- **备注**：正式OUP文章页；另有ICML CompBio开放PDF。

scGPCL – Graph Prototypical Contrastive Learning

1. **Efficacy & principle** – scGPCL encodes cell representations using graph neural networks on a cell–gene graph and employs prototypical contrastive learning to push apart dissimilar pairs while pulling together similar ones.  It shows strong performance on real and simulated data【803175867092150†L190-L197】.
2. **Integration strategy** – Incorporate graph prototypical contrastive pretraining into scMAE: build a cell–gene graph, learn embeddings via contrastive loss with prototypes, then use these embeddings as initialization for masked autoencoding.
3. **Potential issues** – Requires graph construction; prototypes may not align with true biological clusters if the data are noisy.  Contrastive loss may dominate reconstruction loss.
4. **Mitigation** – Use high‑quality similarity measures for graph construction and update prototypes dynamically.  Balance contrastive and reconstruction losses.
5. **Usability decision** – scGPCL’s contrastive framework can improve representation learning; it is usable as a pretraining step.

Citations:【803175867092150†L190-L197】

### 069. IGCLAPS: an interpretable graph contrastive learning method with adaptive positive sampling for scRNA-seq data analysis

- **排序文件**：`01_PDF论文_按推荐程度排序/069_可用_IGCLAPS_an_interpretable_graph_contrastive_learning_method_with_adaptive_positive_sampling_for_scRNA-seq_data_analysis.pdf`
- **推荐等级**：可用（score=70）
- **SHA-256**：`91E511FD31B9D9D81162ED2177532601E7E27CA33A4CF23CD111A05E3EC4C5D9`
- **来源记录**：文库2#30 30_IGCLAPS.pdf

#### 来源说明：文库2 #30

说明来源：文库2/scmae_improvement_reports/30_IGCLAPS.txt
- **标题**：IGCLAPS: an interpretable graph contrastive learning method with adaptive positive sampling for scRNA-seq data analysis
- **年份/会议**：2025 / Bioinformatics
- **备注**：OUP开放PDF。

IGCLAPS – Interpretable Graph Contrastive Learning with Adaptive Positive Sampling

1. **Efficacy & principle** – IGCLAPS uses a graph transformer to learn low‑dimensional embeddings and a dual‑head contrastive module with adaptive positive sampling.  It dynamically identifies true positive pairs based on expression similarity and soft cluster labels, improving clustering and interpretability【733287641060391†L125-L138】.
2. **Integration strategy** – Pretrain scMAE with IGCLAPS: build a cell graph, compute adaptive positive pairs, and use contrastive loss alongside masked reconstruction.  Use the graph transformer to capture complex relationships among cells.
3. **Potential issues** – Determining adaptive positives may be sensitive to initial cluster labels.  Graph transformers are resource intensive.
4. **Mitigation** – Use iterative clustering to update positive pairs.  Employ efficient transformer variants or limit the number of neighbours.
5. **Usability decision** – IGCLAPS provides interpretable contrastive pretraining; integrating it into scMAE is usable with computational considerations.

Citations:【733287641060391†L125-L138】

### 070. Deep clustering of single-cell RNA-seq using adversarial graph contrastive learning

- **排序文件**：`01_PDF论文_按推荐程度排序/070_可用_Deep_clustering_of_single-cell_RNA-seq_using_adversarial_graph_contrastive_learning.pdf`
- **推荐等级**：可用（score=70）
- **SHA-256**：`9D16E9C570FAF0128BBF47A846B6099B56A712149F201954DFBA1C0780FC393B`
- **来源记录**：文库2#31 31_scAGCL.pdf

#### 来源说明：文库2 #31

说明来源：文库2/scmae_improvement_reports/31_scAGCL.txt
- **标题**：Deep clustering of single-cell RNA-seq using adversarial graph contrastive learning
- **年份/会议**：2025 / Briefings in Bioinformatics
- **备注**：OUP开放PDF。

scAGCL – Adversarial Graph Contrastive Learning

1. **Efficacy & principle** – scAGCL constructs a cell‑cell graph and applies contrastive learning with adversarial attacks on graph structure and node features.  Subgraph sampling improves scalability.  It outperforms seven state‑of‑the‑art methods and identifies marker genes【447889245850842†L170-L184】.
2. **Integration strategy** – Use adversarial contrastive pretraining in scMAE: perturb the cell graph and gene features adversarially during training to make embeddings robust.  Combine with masked reconstruction to learn stable representations.
3. **Potential issues** – Adversarial attacks require careful design; over‑perturbation may degrade useful signal.  Training with adversarial examples increases computational cost.
4. **Mitigation** – Control perturbation magnitude and use virtual adversarial training.  Limit the number of adversarial examples per batch.
5. **Usability decision** – Adversarial contrastive learning can improve robustness; integrating scAGCL elements into scMAE is usable with careful tuning.

Citations:【447889245850842†L170-L184】

### 071. scCAN: single-cell clustering using autoencoder and network fusion

- **排序文件**：`01_PDF论文_按推荐程度排序/071_可用_scCAN_single-cell_clustering_using_autoencoder_and_network_fusion.pdf`
- **推荐等级**：可用（score=70）
- **SHA-256**：`97F4E5FED42D03831F9DCCF95FCAA4A25DE000F20C481CDF9285F20E72C412A7`
- **来源记录**：文库2#37 37_scCAN.pdf

#### 来源说明：文库2 #37

说明来源：文库2/scmae_improvement_reports/37_scCAN.txt
- **标题**：scCAN: single-cell clustering using autoencoder and network fusion
- **年份/会议**：2022 / Scientific Reports
- **备注**：Nature开放PDF。

scCAN – Autoencoder and Network Fusion

1. **Efficacy & principle** – scCAN uses a non‑negative kernel autoencoder and a stacked variational autoencoder to generate multiple low‑dimensional representations of scRNA‑seq data, which are then fused via network fusion for clustering.  It accurately estimates the number of cell types, is robust to dropouts and scales to millions of cells【795904186525999†L80-L89】.
2. **Integration strategy** – Use scCAN’s network fusion approach to combine representations from scMAE and other autoencoders.  The fused representation can capture complementary features and improve clustering.
3. **Potential issues** – Network fusion can be complex and may require many intermediate representations.  Non‑negative constraints may not be applicable to all gene expression patterns.
4. **Mitigation** – Use a subset of representations and regularize fusion weights.  Ensure that non‑negative constraints are applied only when appropriate.
5. **Usability decision** – scCAN’s fusion strategy can complement scMAE by combining diverse embeddings; it is usable.

Citations:【795904186525999†L80-L89】

### 072. scASDC: Attention Enhanced Structural Deep Clustering for Single-cell RNA-seq Data

- **排序文件**：`01_PDF论文_按推荐程度排序/072_可用_scASDC_Attention_Enhanced_Structural_Deep_Clustering_for_Single-cell_RNA-seq_Data.pdf`
- **推荐等级**：可用（score=70）
- **SHA-256**：`4AD5F934402F332BF38024D43BF2EB73B0B697928F897097A4384BDD8BF13C7B`
- **来源记录**：文库2#38 38_scASDC.pdf

#### 来源说明：文库2 #38

说明来源：文库2/scmae_improvement_reports/38_scASDC.txt
- **标题**：scASDC: Attention Enhanced Structural Deep Clustering for Single-cell RNA-seq Data
- **年份/会议**：2024 / arXiv
- **备注**：arXiv开放PDF。

scASDC – Attention‑Enhanced Structural Deep Clustering

1. **Efficacy & principle** – scASDC combines a multi‑layer graph convolutional network (GCN) to capture high‑order structural relationships, a ZINB‑based autoencoder to extract content information, an attention fusion mechanism to combine these at each layer and a self‑supervised module for robustness.  It outperforms state‑of‑the‑art methods【995177089489918†L39-L53】.
2. **Integration strategy** – Integrate scASDC’s attention fusion into scMAE: use a GCN to compute structural embeddings and a ZINB autoencoder for gene expression, then fuse them via attention before masking.  A self‑supervised loss can be added to enforce consistency.
3. **Potential issues** – GCNs may suffer from over‑smoothing and require graph construction.  Attention fusion increases model complexity.
4. **Mitigation** – Limit GCN depth, use residual connections and carefully tune attention parameters.  Precompute graphs to reduce overhead.
5. **Usability decision** – scASDC’s fusion of structure and content aligns well with scMAE’s goals; it is usable with careful tuning.

Citations:【995177089489918†L39-L53】

### 073. Attention-based deep clustering method for scRNA-seq cell type identification

- **排序文件**：`01_PDF论文_按推荐程度排序/073_可用_Attention-based_deep_clustering_method_for_scRNA-seq_cell_type_identification.pdf`
- **推荐等级**：可用（score=70）
- **SHA-256**：`05093BAA1BE6641F1F8FB81DB69B667270B258B41622BA8585A5B4884ECCBDD3`
- **去重说明**：合并 2 条来源记录；包含 1 个PDF哈希。
- **来源记录**：文库2#39 39_AttentionAE-sc.pdf；文库3#13 10_AttentionAE_sc_Li_2023_Attention_based_deep_clustering_method.pdf

#### 来源说明：文库2 #39

说明来源：文库2/scmae_improvement_reports/39_AttentionAE_sc.txt
- **标题**：Attention-based deep clustering method for scRNA-seq cell type identification
- **年份/会议**：2023 / PLOS Computational Biology
- **备注**：PLOS开放PDF/printable。

AttentionAE‑sc – Attention‑based Deep Clustering

1. **Efficacy & principle** – AttentionAE‑sc combines zero‑inflated negative binomial (ZINB) denoising and graph autoencoder (GAE) embeddings via an attention mechanism.  Iterative fusion between denoising and topological embeddings produces clustering‑friendly representations and achieves superior performance on 16 datasets【693866387649793†L182-L194】.
2. **Integration strategy** – Use an attention mechanism to fuse scMAE’s reconstructed gene expression with graph‑based embeddings.  This can improve representation quality and clustering without specifying the number of clusters.
3. **Potential issues** – Requires careful balancing between ZINB and GAE contributions.  Attention may overfit to certain gene or graph features.
4. **Mitigation** – Use multi‑head attention and regularization.  Evaluate across datasets to set attention weights adaptively.
5. **Usability decision** – AttentionAE‑sc’s fusion mechanism can enhance scMAE; it is usable.

Citations:【693866387649793†L182-L194】

#### 来源说明：文库3 #13

说明来源：文库3/manifest_30_papers.csv
- **标题**：Attention-based deep clustering method for scRNA-seq cell type identification
- **年份/会议**：PLoS Computational Biology 2023
- **方向/领域**：Single-cell attention clustering
- **与 scMAE 改良相关性**：Attention fusion of denoising and graph embeddings; no need to pre-specify cluster number.

### 074. Boosting scRNA-seq data clustering by cluster-aware feature weighting

- **排序文件**：`01_PDF论文_按推荐程度排序/074_可用_Boosting_scRNA-seq_data_clustering_by_cluster-aware_feature_weighting.pdf`
- **推荐等级**：可用（score=70）
- **SHA-256**：`F48A4B44AE9E849BE89B8B7EAAD1BFAD211D656290B3D71E8EE7F82E09497AA6`
- **来源记录**：文库2#40 40_CaFew.pdf

#### 来源说明：文库2 #40

说明来源：文库2/scmae_improvement_reports/40_CaFew.txt
- **标题**：Boosting scRNA-seq data clustering by cluster-aware feature weighting
- **年份/会议**：2021 / BMC Bioinformatics
- **备注**：Springer开放PDF。

CaFew – Cluster‑Aware Feature Weighting

1. **Efficacy & principle** – CaFew selects genes based on cluster‑aware feature weighting.  It optimizes a clustering objective to obtain a feature weight matrix and selects genes with large weights or high variance across clusters.  Experiments show that CaFew improves clustering performance and visualization when combined with existing methods【346770119344838†L80-L90】.
2. **Integration strategy** – Apply CaFew as a feature selection step before scMAE training: compute feature weights using an initial clustering (e.g., SC3), select highly weighted genes and train scMAE on the reduced gene set.
3. **Potential issues** – Requires an initial clustering method; errors in the initial clustering may affect feature selection.  The method may bias towards genes with high variability, neglecting subtle markers.
4. **Mitigation** – Use multiple clustering algorithms and aggregate feature weights to reduce dependence on any single clustering result.  Include known marker genes regardless of weight.
5. **Usability decision** – CaFew provides an effective gene selection mechanism; it is usable as a preprocessing step.

Citations:【346770119344838†L80-L90】

### 075. scGTN: Deep Siamese Graph Transformer Network for Single-cell RNA Sequencing Clustering

- **排序文件**：`01_PDF论文_按推荐程度排序/075_可用_scGTN_Deep_Siamese_Graph_Transformer_Network_for_Single-cell_RNA_Sequencing_Clustering.pdf`
- **推荐等级**：可用（score=70）
- **SHA-256**：`F5B04362374D061A0594936BDBFB816497F03067C05D8457571D76E5F0F30409`
- **来源记录**：文库2#42 42_scGTN.pdf

#### 来源说明：文库2 #42

说明来源：文库2/scmae_improvement_reports/42_scGTN.txt
- **标题**：scGTN: Deep Siamese Graph Transformer Network for Single-cell RNA Sequencing Clustering
- **年份/会议**：2026 / arXiv
- **备注**：arXiv开放PDF。

scGTN – Deep Siamese Graph Transformer Network

1. **Efficacy & principle** – scGTN integrates gene expression profiles and intercellular structural dependencies via a Siamese graph transformer network.  It constructs two augmented graph views, incorporates shortest‑path information and node distances, and uses an optimal transport strategy for clustering.  The method outperforms existing techniques【191727095808050†L67-L84】.
2. **Integration strategy** – Use scGTN’s Siamese transformer as a pre‑encoder for scMAE: build two graph views (e.g., based on different similarity measures) and learn embeddings capturing rich structural relations.  Use optimal transport or contrastive loss to align the embeddings before masked autoencoding.
3. **Potential issues** – Transformers are computationally intensive and may not scale to large datasets.  Optimal transport adds additional cost.
4. **Mitigation** – Use efficient transformer architectures (e.g., sparse attention).  Approximate optimal transport with faster algorithms or limit its use to small batches.
5. **Usability decision** – scGTN’s ability to capture complex graph structure makes it valuable; it is usable with computational optimizations.

Citations:【191727095808050†L67-L84】

### 076. scGraphformer: unveiling cellular heterogeneity and interactions in scRNA-seq data using a scalable graph transformer network

- **排序文件**：`01_PDF论文_按推荐程度排序/076_可用_scGraphformer_unveiling_cellular_heterogeneity_and_interactions_in_scRNA-seq_data_using_a_scalable_graph_transformer_network.pdf`
- **推荐等级**：可用（score=70）
- **SHA-256**：`7318C8A9234CF9F5643209D7ED3AE68EE7F98FF61259653AC0E135AA44F03051`
- **来源记录**：文库2#45 45_scGraphformer.pdf

#### 来源说明：文库2 #45

说明来源：文库2/scmae_improvement_reports/45_scGraphformer.txt
- **标题**：scGraphformer: unveiling cellular heterogeneity and interactions in scRNA-seq data using a scalable graph transformer network
- **年份/会议**：2024 / Communications Biology
- **备注**：Nature开放PDF。

scGraphformer – Graph Transformer for scRNA‑seq

1. **Efficacy & principle** – scGraphformer learns an all‑encompassing cell‑cell relational network directly from data via iterative refinement, constructing a dense graph capturing full interactions.  It outperforms existing methods in cell type identification and reveals subtle cell patterns and interactions【435945997048765†L76-L91】.
2. **Integration strategy** – Integrate scGraphformer’s dense relational network into scMAE: refine an initial cell graph iteratively and use the resulting dense adjacency to inform masking and reconstruction.  The transformer can model global dependencies among cells.
3. **Potential issues** – Dense graphs are memory‑intensive and may lead to over‑fitting.  Transformer models require large computational resources.
4. **Mitigation** – Use sparse attention or kernel approximations.  Prune the dense graph based on thresholding or sampling.
5. **Usability decision** – scGraphformer’s global interaction modelling is powerful; it is usable with resource considerations.

Citations:【435945997048765†L76-L91】

### 077. scNET: learning context-specific gene and cell embeddings by integrating single-cell gene expression data with protein-protein interactions

- **排序文件**：`01_PDF论文_按推荐程度排序/077_可用_scNET_learning_context-specific_gene_and_cell_embeddings_by_integrating_single-cell_gene_expression_data_with_protein-protein_inte.pdf`
- **推荐等级**：可用（score=70）
- **SHA-256**：`31838710C9DF9DFA2AA9302F85482AD305DC56678A8A1B303FC2D92431797FAF`
- **来源记录**：文库2#46 46_scNET.pdf

#### 来源说明：文库2 #46

说明来源：文库2/scmae_improvement_reports/46_scNET.txt
- **标题**：scNET: learning context-specific gene and cell embeddings by integrating single-cell gene expression data with protein-protein interactions
- **年份/会议**：2025 / Nature Methods
- **备注**：Nature开放PDF。

scNET – Dual‑View GNN with Protein–Protein Interaction Integration

1. **Efficacy & principle** – scNET integrates scRNA‑seq data with protein–protein interaction networks using a dual‑view GNN architecture.  It models gene–gene and cell–cell relationships, capturing gene annotations and pathways and improving clustering and pathway analysis across conditions【530557996525877†L69-L84】.
2. **Integration strategy** – Incorporate protein–protein interaction information into scMAE: use a dual GNN to capture gene–gene network features alongside cell–cell features.  Combine these features with the masked autoencoder to improve biological interpretability.
3. **Potential issues** – Requires high‑quality PPI networks; integrating multiple networks increases complexity.  Some genes may lack PPI information.
4. **Mitigation** – Use curated and confidence‑weighted PPI networks.  For genes without interactions, rely on gene expression similarity.  Apply regularization to prevent overfitting to the PPI network.
5. **Usability decision** – Integrating PPI data can enhance biological relevance; scNET’s dual‑view approach is usable for scMAE.

Citations:【530557996525877†L69-L84】

### 078. scLong: a billion-parameter foundation model for capturing long-range gene context in single-cell transcriptomics

- **排序文件**：`01_PDF论文_按推荐程度排序/078_可用_scLong_a_billion-parameter_foundation_model_for_capturing_long-range_gene_context_in_single-cell_transcriptomics.pdf`
- **推荐等级**：可用（score=70）
- **SHA-256**：`69391EEA53971B91D571A1730225EB71BDA91CF11EC1C4D2C5E73939AEC0F5FC`
- **来源记录**：文库2#47 47_scLong.pdf

#### 来源说明：文库2 #47

说明来源：文库2/scmae_improvement_reports/47_scLong.txt
- **标题**：scLong: a billion-parameter foundation model for capturing long-range gene context in single-cell transcriptomics
- **年份/会议**：2026 / Nature Communications
- **备注**：Nature开放PDF。

scLong – Single‑Cell Foundation Model

1. **Efficacy & principle** – scLong is a billion‑parameter transformer pretrained on 48 million cells using self‑attention over 28 000 genes.  It captures long‑range gene dependencies and integrates gene ontology knowledge via graph convolution, outperforming state‑of‑the‑art models【389129090882071†L87-L104】.
2. **Integration strategy** – Use scLong as a pretrained encoder to initialize scMAE.  Its ability to model long‑range dependencies can improve representation quality.  Fine‑tune scMAE on task‑specific data while retaining scLong’s knowledge.
3. **Potential issues** – The model is extremely large and requires significant computational resources.  Fine‑tuning may require careful learning rate schedules.
4. **Mitigation** – Use parameter‑efficient fine‑tuning methods (e.g., adapters or LoRA).  Distil scLong into a smaller model before integrating with scMAE.
5. **Usability decision** – scLong provides powerful representations; with resource‑aware fine‑tuning, it is usable for scMAE improvement.

Citations:【389129090882071†L87-L104】

### 079. scHNTL: single-cell RNA-seq data clustering augmented by high-order neighbors and triplet loss

- **排序文件**：`01_PDF论文_按推荐程度排序/079_可用_scHNTL_single-cell_RNA-seq_data_clustering_augmented_by_high-order_neighbors_and_triplet_loss.pdf`
- **推荐等级**：可用（score=70）
- **SHA-256**：`39AB3C64357544C6E9EEFEE2ECC48F3DF4010B7705C3CA5D5109960637015023`
- **来源记录**：文库2#48 48_scHNTL.pdf

#### 来源说明：文库2 #48

说明来源：文库2/scmae_improvement_reports/48_scHNTL.txt
- **标题**：scHNTL: single-cell RNA-seq data clustering augmented by high-order neighbors and triplet loss
- **年份/会议**：2025 / Bioinformatics
- **备注**：OUP开放PDF。

scHNTL – High‑Order Neighbor Triplet Learning

1. **Efficacy & principle** – scHNTL builds a similarity graph and uses a graph attentional autoencoder to learn embeddings.  It employs a triplet loss to separate dissimilar cell pairs and fuses improved embeddings with clustering objectives, outperforming existing methods【665144320614849†L134-L156】.
2. **Integration strategy** – Use triplet loss within scMAE to enforce that cells of the same type have closer embeddings than dissimilar cells.  Combine graph attentional encoder with masked reconstruction to capture both structural and expression information.
3. **Potential issues** – Triplet loss requires careful selection of positive and negative samples.  Overemphasis on triplet loss may harm reconstruction.
4. **Mitigation** – Use hard negative mining and weight triplet loss relative to reconstruction loss.  Balance sampling to ensure diversity.
5. **Usability decision** – Triplet supervision can improve separation of cell types; scHNTL’s approach is usable for scMAE enhancement.

Citations:【665144320614849†L134-L156】

### 080. M-band wavelet-based multi-view clustering of cells

- **排序文件**：`01_PDF论文_按推荐程度排序/080_可用_M-band_wavelet-based_multi-view_clustering_of_cells.pdf`
- **推荐等级**：可用（score=70）
- **SHA-256**：`C733F7A06AD435EFBFBF6B1C95CE635E04845F2D2EA4726E5DC4A12D533FB669`
- **来源记录**：文库2#49 49_WMC.pdf

#### 来源说明：文库2 #49

说明来源：文库2/scmae_improvement_reports/49_WMC.txt
- **标题**：M-band wavelet-based multi-view clustering of cells
- **年份/会议**：2025 / PLOS Computational Biology
- **备注**：PLOS开放PDF/printable。

WMC – Wavelet‑based Multi‑view Clustering

1. **Efficacy & principle** – WMC integrates M‑band discrete wavelet analysis and UMAP by decomposing the log‑transformed scRNA‑seq data into low‑ and high‑resolution components.  It performs multi‑view clustering to visualise missing cell types and discover new cell types, offering fine resolution and better rare cell identification【75457080070267†L317-L328】.
2. **Integration strategy** – Use wavelet decomposition as a preprocessing step for scMAE: decompose gene expression into multiple frequency components, train separate masked autoencoders on each view and fuse the representations for final clustering.
3. **Potential issues** – Wavelet transforms may be sensitive to noise and require careful choice of wavelet basis.  Multi‑view training increases computation.
4. **Mitigation** – Use robust wavelet bases and denoising.  Select a small number of informative scales to reduce complexity.
5. **Usability decision** – Incorporating wavelet‑based multi‑view representations can enhance scMAE’s ability to capture patterns at different resolutions; it is usable.

Citations:【75457080070267†L317-L328】

### 081. OKR-CELL: Open-world Language Knowledge-Aided Robust Single-Cell Foundation Model

- **排序文件**：`01_PDF论文_按推荐程度排序/081_谨慎_部分可用_OKR-CELL_Open-world_Language_Knowledge-Aided_Robust_Single-Cell_Foundation_Model.pdf`
- **推荐等级**：谨慎/部分可用（score=55）
- **SHA-256**：`3FA3C481FBE15E348C776FEAD5D7C1F831EAB5493EF2F23291E75EE9F3F9FE51`
- **来源记录**：文库1#7 07_OKR_CELL_Wang_2026_arxiv.pdf

#### 来源说明：文库1 #7

说明来源：文库1/manifest_extra50_papers.csv
- **标题**：OKR-CELL: Open-world Language Knowledge-Aided Robust Single-Cell Foundation Model
- **年份/会议**：2026 / arXiv / bioRxiv withdrawn version noted
- **方向/领域**：RAG + robust cell-language alignment
- **核心思想**：用RAG扩展细胞文本描述，并用可靠性评估、课程学习和momentum contrastive实现稳健跨模态对齐。
- **如何改良 scMAE**：借鉴CRA：对scMAE的cell-text或cell-ontology对齐样本估计可靠性，低质量文本/伪标签降权。
- **有效性依据**：中：思想有价值，但论文版本和状态需谨慎。
- **潜在问题**：来源有withdrawn记录；结果需复核。
- **规避方案**：只采纳“可靠性加权对齐”机制，不直接依赖其具体指标。
- **最终可用性**：可用但需谨慎引用。
- **来源链接**：https://arxiv.org/abs/2601.05648

### 082. DiT: Scalable Diffusion Models with Transformers

- **排序文件**：`01_PDF论文_按推荐程度排序/082_谨慎_部分可用_DiT_Scalable_Diffusion_Models_with_Transformers.pdf`
- **推荐等级**：谨慎/部分可用（score=55）
- **SHA-256**：`5FBB0EC35E3AE76826240171CA63ED4C27DE6CB43180367453019E9D28A8B076`
- **来源记录**：文库1#31 31_DiT_Peebles_2023_ICCV_arxiv.pdf

#### 来源说明：文库1 #31

说明来源：文库1/manifest_extra50_papers.csv
- **标题**：DiT: Scalable Diffusion Models with Transformers
- **年份/会议**：2023 / ICCV
- **方向/领域**：diffusion transformer
- **核心思想**：用Transformer替代U-Net做latent diffusion，并分析scaling。
- **如何改良 scMAE**：若做scRNA latent diffusion，可用Transformer/DiT条件生成扰动响应。
- **有效性依据**：中。
- **潜在问题**：计算重，偏生成任务。
- **规避方案**：只作为扩展方向或小模型。
- **最终可用性**：可用但非首选。
- **来源链接**：https://arxiv.org/abs/2212.09748

### 083. scGNN is a novel graph neural network framework for single-cell RNA-seq analyses

- **排序文件**：`01_PDF论文_按推荐程度排序/083_未明确_scGNN_is_a_novel_graph_neural_network_framework_for_single-cell_RNA-seq_analyses.pdf`
- **推荐等级**：未明确（score=50）
- **SHA-256**：`2688130DEA51FE52B575F44E246379BED5E09F1EB846A3227B6019CEE6DF5F60`
- **去重说明**：合并 2 条来源记录；包含 1 个PDF哈希。
- **来源记录**：文库2#44 44_scGNN.pdf；文库3#10 08_scGNN_Wang_2021_Novel_graph_neural_network_framework_for_scRNA_seq.pdf

#### 来源说明：文库2 #44

说明来源：文库2/scmae_improvement_reports/44_scGNN.txt
- **标题**：scGNN is a novel graph neural network framework for single-cell RNA-Seq analyses
- **年份/会议**：2021 / Nature Communications
- **备注**：Nature开放PDF。

scGNN – Graph Neural Network Framework

1. **Efficacy & principle** – scGNN formulates cell–cell relationships with graph neural networks and models gene expression with a left‑truncated mixture Gaussian.  It integrates multi‑modal autoencoders and outperforms other tools for imputation and clustering【910622243095874†L84-L98】.
2. **Integration strategy** – Use scGNN to construct a cell graph and generate embeddings for scMAE.  The left‑truncated Gaussian can model gene expression distribution in the decoder.  Multi‑modal autoencoders provide a template for integrating additional data types.
3. **Potential issues** – Graph neural networks require graph construction and may scale poorly.  Mixture Gaussian assumptions may not fit all genes.
4. **Mitigation** – Use approximate nearest neighbours for graph construction and adopt hierarchical GNN architectures.  Fit gene‑specific distributions.
5. **Usability decision** – scGNN provides a strong baseline for graph representation; its concepts are usable in scMAE.

Citations:【910622243095874†L84-L98】

#### 来源说明：文库3 #10

说明来源：文库3/manifest_30_papers.csv
- **标题**：scGNN is a novel graph neural network framework for single-cell RNA-seq analyses
- **年份/会议**：Nature Communications 2021
- **方向/领域**：Single-cell GNN
- **与 scMAE 改良相关性**：Cell-cell graph

### 084. Unsupervised Deep Embedding for Clustering Analysis

- **排序文件**：`01_PDF论文_按推荐程度排序/084_未明确_Unsupervised_Deep_Embedding_for_Clustering_Analysis.pdf`
- **推荐等级**：未明确（score=50）
- **SHA-256**：`D0D6B668038FD8D365C89BEBE6BFDFBE31E611A71A803E95A0ACFDD8B18F6E3F`
- **来源记录**：文库3#1 01_DEC_Xie_2016_Unsupervised_Deep_Embedding_for_Clustering_Analysis.pdf

#### 来源说明：文库3 #1

说明来源：文库3/manifest_30_papers.csv
- **标题**：Unsupervised Deep Embedding for Clustering Analysis
- **年份/会议**：ICML 2016
- **方向/领域**：Deep clustering
- **与 scMAE 改良相关性**：DEC-style KL self-training and target distribution can improve scMAE clustering head.

### 085. BERT: Pre-training of Deep Bidirectional Transformers for Language Understanding

- **排序文件**：`01_PDF论文_按推荐程度排序/085_未明确_BERT_Pre-training_of_Deep_Bidirectional_Transformers_for_Language_Understanding.pdf`
- **推荐等级**：未明确（score=50）
- **SHA-256**：`5692A5514787A8C6727B4FF3B726A3385798BC68E12138D1D4AF83947E2ACF6E`
- **来源记录**：文库3#4 04_BERT_Devlin_2019_NAACL.pdf

#### 来源说明：文库3 #4

说明来源：文库3/manifest_30_papers.csv
- **标题**：BERT: Pre-training of Deep Bidirectional Transformers for Language Understanding
- **年份/会议**：NAACL 2019
- **方向/领域**：NLP masked modeling
- **与 scMAE 改良相关性**：Masked-token pretraining and bidirectional context inspire gene-token modeling.

### 086. A Simple Framework for Contrastive Learning of Visual Representations

- **排序文件**：`01_PDF论文_按推荐程度排序/086_未明确_A_Simple_Framework_for_Contrastive_Learning_of_Visual_Representations.pdf`
- **推荐等级**：未明确（score=50）
- **SHA-256**：`E1B10DB76ADEB0D21A014B9ABCF42F6A6AF928540C223F010FA1D8C293D9D105`
- **来源记录**：文库3#5 05_SimCLR_Chen_2020_ICML.pdf

#### 来源说明：文库3 #5

说明来源：文库3/manifest_30_papers.csv
- **标题**：A Simple Framework for Contrastive Learning of Visual Representations
- **年份/会议**：ICML 2020
- **方向/领域**：Contrastive SSL
- **与 scMAE 改良相关性**：Augmentation design and instance discrimination for robust cell embeddings.

### 087. Bootstrap Your Own Latent: A New Approach to Self-Supervised Learning

- **排序文件**：`01_PDF论文_按推荐程度排序/087_未明确_Bootstrap_Your_Own_Latent_A_New_Approach_to_Self-Supervised_Learning.pdf`
- **推荐等级**：未明确（score=50）
- **SHA-256**：`873E7B74D58E1E17806CF7A3CF6B80EAD914559F36ED19313C749672EED0AA94`
- **来源记录**：文库3#6 06_BYOL_Grill_2020_NeurIPS.pdf

#### 来源说明：文库3 #6

说明来源：文库3/manifest_30_papers.csv
- **标题**：Bootstrap Your Own Latent: A New Approach to Self-Supervised Learning
- **年份/会议**：NeurIPS 2020
- **方向/领域**：Non-contrastive SSL
- **与 scMAE 改良相关性**：Avoids false negatives

### 088. Emerging Properties in Self-Supervised Vision Transformers

- **排序文件**：`01_PDF论文_按推荐程度排序/088_未明确_Emerging_Properties_in_Self-Supervised_Vision_Transformers.pdf`
- **推荐等级**：未明确（score=50）
- **SHA-256**：`AA464CFD59A428890190BEA1065823A491A853478B4FB2A25F5EB44D442CE296`
- **来源记录**：文库3#7 07_DINO_Caron_2021_ICCV.pdf

#### 来源说明：文库3 #7

说明来源：文库3/manifest_30_papers.csv
- **标题**：Emerging Properties in Self-Supervised Vision Transformers
- **年份/会议**：ICCV 2021
- **方向/领域**：Teacher-student SSL
- **与 scMAE 改良相关性**：EMA teacher and centering/sharpening reduce collapse; useful for stable scMAE representations.

### 089. scziDesk: deep soft K-means clustering with self-training

- **排序文件**：`01_PDF论文_按推荐程度排序/089_未明确_scziDesk_deep_soft_K-means_clustering_with_self-training.pdf`
- **推荐等级**：未明确（score=50）
- **SHA-256**：`4EC21E3C4E7215533B36F4972F807D530FB031BC400F2487D38E9C8B3AF2277E`
- **来源记录**：文库3#8 07_scziDesk_Chen_2020_Deep_soft_K_means_clustering_with_self_training.pdf

#### 来源说明：文库3 #8

说明来源：文库3/manifest_30_papers.csv
- **标题**：scziDesk: deep soft K-means clustering with self-training
- **年份/会议**：NAR Genomics and Bioinformatics 2020
- **方向/领域**：Single-cell deep clustering
- **与 scMAE 改良相关性**：Self-training and ZINB variants for scRNA clustering.

### 090. ELECTRA: Pre-training Text Encoders as Discriminators Rather Than Generators

- **排序文件**：`01_PDF论文_按推荐程度排序/090_未明确_ELECTRA_Pre-training_Text_Encoders_as_Discriminators_Rather_Than_Generators.pdf`
- **推荐等级**：未明确（score=50）
- **SHA-256**：`6219A802DCA3BE0027FC88BB2ECD56806E6610EBFE4D083A4C18DD896F45B9AE`
- **来源记录**：文库3#9 08_ELECTRA_Clark_2020_ICLR.pdf

#### 来源说明：文库3 #9

说明来源：文库3/manifest_30_papers.csv
- **标题**：ELECTRA: Pre-training Text Encoders as Discriminators Rather Than Generators
- **年份/会议**：ICLR 2020
- **方向/领域**：Discriminative pretraining
- **与 scMAE 改良相关性**：Replaced-token detection can be adapted to replaced-expression/gene corruption detection.

### 091. GraphMAE: Self-Supervised Masked Graph Autoencoders

- **排序文件**：`01_PDF论文_按推荐程度排序/091_未明确_GraphMAE_Self-Supervised_Masked_Graph_Autoencoders.pdf`
- **推荐等级**：未明确（score=50）
- **SHA-256**：`F297CB1597B086DE05B46699854E0665A7DC74BFC1F6559B868A8AE6AC3A6DA2`
- **来源记录**：文库3#11 09_GraphMAE_Hou_2022_KDD.pdf

#### 来源说明：文库3 #11

说明来源：文库3/manifest_30_papers.csv
- **标题**：GraphMAE: Self-Supervised Masked Graph Autoencoders
- **年份/会议**：KDD 2022
- **方向/领域**：Graph SSL
- **与 scMAE 改良相关性**：Attribute masking on graphs can add topology-aware pretraining to scMAE.

### 092. scDSC: Deep structural clustering jointly through autoencoder and GNN

- **排序文件**：`01_PDF论文_按推荐程度排序/092_未明确_scDSC_Deep_structural_clustering_jointly_through_autoencoder_and_GNN.pdf`
- **推荐等级**：未明确（score=50）
- **SHA-256**：`69C0C7EFF8B5607EF05A8AA3359876E9AAF288EEF7B4B4A0DDBE1D966191306A`
- **来源记录**：文库3#12 09_scDSC_Gan_2022_Deep_structural_clustering_AE_GNN.pdf

#### 来源说明：文库3 #12

说明来源：文库3/manifest_30_papers.csv
- **标题**：scDSC: Deep structural clustering jointly through autoencoder and GNN
- **年份/会议**：Briefings in Bioinformatics 2022
- **方向/领域**：Single-cell GNN clustering
- **与 scMAE 改良相关性**：ZINB/NB autoencoder + GNN + mutual supervision for structural clustering.

### 093. GraphMAE2: A Decoding-Enhanced Masked Self-Supervised Graph Learner

- **排序文件**：`01_PDF论文_按推荐程度排序/093_未明确_GraphMAE2_A_Decoding-Enhanced_Masked_Self-Supervised_Graph_Learner.pdf`
- **推荐等级**：未明确（score=50）
- **SHA-256**：`B71BFFFE1CDD48C5EB5BCDE9E59AA43072737A96792E07AE9A7ABA6B212580BE`
- **来源记录**：文库3#14 10_GraphMAE2_Hou_2023_WWW.pdf

#### 来源说明：文库3 #14

说明来源：文库3/manifest_30_papers.csv
- **标题**：GraphMAE2: A Decoding-Enhanced Masked Self-Supervised Graph Learner
- **年份/会议**：WWW 2023
- **方向/领域**：Graph SSL
- **与 scMAE 改良相关性**：Enhanced graph decoding and masking strategy for robust graph representation.

### 094. Graph Contrastive Learning with Augmentations

- **排序文件**：`01_PDF论文_按推荐程度排序/094_未明确_Graph_Contrastive_Learning_with_Augmentations.pdf`
- **推荐等级**：未明确（score=50）
- **SHA-256**：`DC72862B278B632D7D31852D00561C3F333B06A05061B7C863DC11E944120FFA`
- **来源记录**：文库3#15 11_GraphCL_You_2020_NeurIPS.pdf

#### 来源说明：文库3 #15

说明来源：文库3/manifest_30_papers.csv
- **标题**：Graph Contrastive Learning with Augmentations
- **年份/会议**：NeurIPS 2020
- **方向/领域**：Graph contrastive SSL
- **与 scMAE 改良相关性**：Cell graph augmentations and graph contrastive loss for robustness.

### 095. scCDCG: Deep Cut-informed Graph Embedding for scRNA-seq clustering

- **排序文件**：`01_PDF论文_按推荐程度排序/095_未明确_scCDCG_Deep_Cut-informed_Graph_Embedding_for_scRNA-seq_clustering.pdf`
- **推荐等级**：未明确（score=50）
- **SHA-256**：`3371F6EDB52541770374F388DA26A05D04DBB2519E3FB64B6FE3A418FA658961`
- **来源记录**：文库3#16 11_scCDCG_Xu_2024_Deep_Cut_Informed_Graph_Embedding.pdf

#### 来源说明：文库3 #16

说明来源：文库3/manifest_30_papers.csv
- **标题**：scCDCG: Deep Cut-informed Graph Embedding for scRNA-seq clustering
- **年份/会议**：Recent preprint/relevant scRNA method
- **方向/领域**：Single-cell graph cuts
- **与 scMAE 改良相关性**：Normalized-cut-style structural regularization and optimal transport assignment.

### 096. From Louvain to Leiden: guaranteeing well-connected communities

- **排序文件**：`01_PDF论文_按推荐程度排序/096_未明确_From_Louvain_to_Leiden_guaranteeing_well-connected_communities.pdf`
- **推荐等级**：未明确（score=50）
- **SHA-256**：`9475A1783879A5B7E33790AB34609267FB80036E5BF1E1869E584BA31A246749`
- **来源记录**：文库3#19 14_Leiden_Traag_2019_From_Louvain_to_Leiden.pdf

#### 来源说明：文库3 #19

说明来源：文库3/manifest_30_papers.csv
- **标题**：From Louvain to Leiden: guaranteeing well-connected communities
- **年份/会议**：Scientific Reports 2019
- **方向/领域**：Community detection
- **与 scMAE 改良相关性**：Leiden guarantees connected communities and better partitions.

### 097. Comprehensive integration of single-cell data

- **排序文件**：`01_PDF论文_按推荐程度排序/097_未明确_Comprehensive_integration_of_single-cell_data.pdf`
- **推荐等级**：未明确（score=50）
- **SHA-256**：`EE0606CBB60B4E0958C23BC356ACFD4040A436D57787C1DDD43D76ED55F4E5BE`
- **来源记录**：文库3#20 15_Seurat_Stuart_2019_Comprehensive_integration_of_single_cell_data_PREPRINT.pdf

#### 来源说明：文库3 #20

说明来源：文库3/manifest_30_papers.csv
- **标题**：Comprehensive integration of single-cell data
- **年份/会议**：Cell 2019 / preprint
- **方向/领域**：Single-cell integration
- **与 scMAE 改良相关性**：Anchors and transfer learning for cross-batch/multimodal integration.

### 098. Graph Attention Networks

- **排序文件**：`01_PDF论文_按推荐程度排序/098_未明确_Graph_Attention_Networks.pdf`
- **推荐等级**：未明确（score=50）
- **SHA-256**：`A7D811A44E1C65BF27DA3A01859F23A8EDBE33B8F96114B4B5066AEDB078F284`
- **来源记录**：文库3#21 16_GAT_Velickovic_2018_ICLR.pdf

#### 来源说明：文库3 #21

说明来源：文库3/manifest_30_papers.csv
- **标题**：Graph Attention Networks
- **年份/会议**：ICLR 2018
- **方向/领域**：Graph neural networks
- **与 scMAE 改良相关性**：Learnable neighbor weighting for cell-cell graph message passing.

### 099. DropEdge: Towards Deep Graph Convolutional Networks on Node Classification

- **排序文件**：`01_PDF论文_按推荐程度排序/099_未明确_DropEdge_Towards_Deep_Graph_Convolutional_Networks_on_Node_Classification.pdf`
- **推荐等级**：未明确（score=50）
- **SHA-256**：`0855A55C3D125DD50705F3F36B03083A468ACE85A3A9EC8425F8B2362DEBE3F7`
- **来源记录**：文库3#22 17_DropEdge_Rong_2020_ICLR.pdf

#### 来源说明：文库3 #22

说明来源：文库3/manifest_30_papers.csv
- **标题**：DropEdge: Towards Deep Graph Convolutional Networks on Node Classification
- **年份/会议**：ICLR 2020
- **方向/领域**：Graph regularization
- **与 scMAE 改良相关性**：Edge dropout prevents over-smoothing and improves generalization.

### 100. Denoising Diffusion Probabilistic Models

- **排序文件**：`01_PDF论文_按推荐程度排序/100_未明确_Denoising_Diffusion_Probabilistic_Models.pdf`
- **推荐等级**：未明确（score=50）
- **SHA-256**：`AEE5E07A802E8DFD2A386374C94FD61D1D056CB7E1E0FEC4F28E8120FF5D8505`
- **来源记录**：文库3#23 18_DDPM_Ho_2020_NeurIPS.pdf

#### 来源说明：文库3 #23

说明来源：文库3/manifest_30_papers.csv
- **标题**：Denoising Diffusion Probabilistic Models
- **年份/会议**：NeurIPS 2020
- **方向/领域**：Generative models
- **与 scMAE 改良相关性**：Diffusion-inspired denoising and rare-cell augmentation possibilities.

### 101. Neural Discrete Representation Learning

- **排序文件**：`01_PDF论文_按推荐程度排序/101_未明确_Neural_Discrete_Representation_Learning.pdf`
- **推荐等级**：未明确（score=50）
- **SHA-256**：`D62036803C190A82718B718F59A87ADB2BE8B586E4887A394AD0CEF15C56301E`
- **来源记录**：文库3#24 19_VQVAE_vanDenOord_2017_NeurIPS.pdf

#### 来源说明：文库3 #24

说明来源：文库3/manifest_30_papers.csv
- **标题**：Neural Discrete Representation Learning
- **年份/会议**：NeurIPS 2017
- **方向/领域**：Generative representation learning
- **与 scMAE 改良相关性**：Discrete latent codebook can model cell-state prototypes.

### 102. mixup: Beyond Empirical Risk Minimization

- **排序文件**：`01_PDF论文_按推荐程度排序/102_未明确_mixup_Beyond_Empirical_Risk_Minimization.pdf`
- **推荐等级**：未明确（score=50）
- **SHA-256**：`D477FC9B3C5232669126D84406F17AE8261D19DB8EBDDFCB25AD995CA09B47E3`
- **来源记录**：文库3#25 20_Mixup_Zhang_2018_ICLR.pdf

#### 来源说明：文库3 #25

说明来源：文库3/manifest_30_papers.csv
- **标题**：mixup: Beyond Empirical Risk Minimization
- **年份/会议**：ICLR 2018
- **方向/领域**：Regularization/augmentation
- **与 scMAE 改良相关性**：Manifold interpolation for robust embeddings and pseudo-label smoothing.

### 103. AutoAugment: Learning Augmentation Strategies from Data

- **排序文件**：`01_PDF论文_按推荐程度排序/103_未明确_AutoAugment_Learning_Augmentation_Strategies_from_Data.pdf`
- **推荐等级**：未明确（score=50）
- **SHA-256**：`63E1B703FF7BD955EBCAE0C3086506F2A092EC5B03F46DCFBD2C658107366C56`
- **来源记录**：文库3#26 21_AutoAugment_Cubuk_2019_CVPR.pdf

#### 来源说明：文库3 #26

说明来源：文库3/manifest_30_papers.csv
- **标题**：AutoAugment: Learning Augmentation Strategies from Data
- **年份/会议**：CVPR 2019
- **方向/领域**：RL/augmentation search
- **与 scMAE 改良相关性**：Search mask/corruption policies for scMAE.

### 104. VIME: Extending the Success of Self- and Semi-supervised Learning to Tabular Domain

- **排序文件**：`01_PDF论文_按推荐程度排序/104_未明确_VIME_Extending_the_Success_of_Self-_and_Semi-supervised_Learning_to_Tabular_Domain.pdf`
- **推荐等级**：未明确（score=50）
- **SHA-256**：`CE2E3863F8E1655CD7E4EEF354C44BB3A918ECEA24ED3E9FDDD7E63AEE858B15`
- **来源记录**：文库3#27 22_VIME_Yoon_2020_NeurIPS.pdf

#### 来源说明：文库3 #27

说明来源：文库3/manifest_30_papers.csv
- **标题**：VIME: Extending the Success of Self- and Semi-supervised Learning to Tabular Domain
- **年份/会议**：NeurIPS 2020
- **方向/领域**：Tabular SSL
- **与 scMAE 改良相关性**：Mask estimation + reconstruction for tabular data; close to scMAE design.

### 105. Deep generative modeling for single-cell transcriptomics

- **排序文件**：`01_PDF论文_按推荐程度排序/105_未明确_Deep_generative_modeling_for_single-cell_transcriptomics.pdf`
- **推荐等级**：未明确（score=50）
- **SHA-256**：`EDFAECE53E86FE4DDE79528020D642DE3B86F389DCD9A665F471F991C1ACE3EB`
- **来源记录**：文库3#28 23_scVI_Lopez_2018_NatureMethods.pdf

#### 来源说明：文库3 #28

说明来源：文库3/manifest_30_papers.csv
- **标题**：Deep generative modeling for single-cell transcriptomics
- **年份/会议**：Nature Methods 2018
- **方向/领域**：Single-cell generative modeling
- **与 scMAE 改良相关性**：Probabilistic latent variable model; uncertainty/batch/count modeling.

### 106. Score-Based Generative Modeling through Stochastic Differential Equations

- **排序文件**：`01_PDF论文_按推荐程度排序/106_未明确_Score-Based_Generative_Modeling_through_Stochastic_Differential_Equations.pdf`
- **推荐等级**：未明确（score=50）
- **SHA-256**：`510EFA34D3FF6DBC93B0A09A58272127D8C6ECF656363890E1A823636F83CCF0`
- **来源记录**：文库3#29 24_ScoreSDE_Song_2021_ICLR.pdf

#### 来源说明：文库3 #29

说明来源：文库3/manifest_30_papers.csv
- **标题**：Score-Based Generative Modeling through Stochastic Differential Equations
- **年份/会议**：ICLR 2021
- **方向/领域**：Generative models
- **与 scMAE 改良相关性**：Continuous-time noise/denoising principle for robust generative augmentation.

### 107. TabNet: Attentive Interpretable Tabular Learning

- **排序文件**：`01_PDF论文_按推荐程度排序/107_未明确_TabNet_Attentive_Interpretable_Tabular_Learning.pdf`
- **推荐等级**：未明确（score=50）
- **SHA-256**：`8310BC2F2C0AAD63FB97C8438F051D93C4496B43FC0123DABD6662CC1EA3FB44`
- **来源记录**：文库3#30 25_TabNet_Arik_2021_AAAI.pdf

#### 来源说明：文库3 #30

说明来源：文库3/manifest_30_papers.csv
- **标题**：TabNet: Attentive Interpretable Tabular Learning
- **年份/会议**：AAAI 2021
- **方向/领域**：Tabular attention
- **与 scMAE 改良相关性**：Sequential feature selection and masked self-supervised tabular learning for gene selection.

### 108. SimMIM: A Simple Framework for Masked Image Modeling

- **排序文件**：`01_PDF论文_按推荐程度排序/108_基线_背景_SimMIM_A_Simple_Framework_for_Masked_Image_Modeling.pdf`
- **推荐等级**：基线/背景（score=25）
- **SHA-256**：`822194C85F12B25A1DD50D19C22444BEB1C22C130B539EBC501BBEFAC084720E`
- **来源记录**：文库1#18 18_SimMIM_Xie_2022_CVPR_arxiv.pdf

#### 来源说明：文库1 #18

说明来源：文库1/manifest_extra50_papers.csv
- **标题**：SimMIM: A Simple Framework for Masked Image Modeling
- **年份/会议**：2022 / CVPR
- **方向/领域**：simple MIM / raw reconstruction
- **核心思想**：证明简单mask+raw reconstruction即可有效。
- **如何改良 scMAE**：为scMAE提供极简强基线：不要过度复杂化，先验证mask比例/decoder深度。
- **有效性依据**：高：可作为基线/消融。
- **潜在问题**：创新性不足。
- **规避方案**：作为baseline而非主贡献。
- **最终可用性**：可用。
- **来源链接**：https://arxiv.org/pdf/2111.09886

### 109. Revisiting Deep Learning Models for Tabular Data / FT-Transformer

- **排序文件**：`01_PDF论文_按推荐程度排序/109_基线_背景_Revisiting_Deep_Learning_Models_for_Tabular_Data_FT-Transformer.pdf`
- **推荐等级**：基线/背景（score=25）
- **SHA-256**：`F2FAA32FFD8C15ED8534677F978DC02730FDBE1F94E9C9557ADA6AA19163478D`
- **来源记录**：文库1#48 48_FTTransformer_Gorishniy_2021_NeurIPS_openreview.pdf

#### 来源说明：文库1 #48

说明来源：文库1/manifest_extra50_papers.csv
- **标题**：Revisiting Deep Learning Models for Tabular Data / FT-Transformer
- **年份/会议**：2021 / NeurIPS
- **方向/领域**：tabular transformer
- **核心思想**：系统评估MLP/ResNet/FT-Transformer在表格数据中的表现。
- **如何改良 scMAE**：作为scMAE backbone替代/基线，检验Transformer是否必要。
- **有效性依据**：中：保证实验严谨。
- **潜在问题**：不一定创新。
- **规避方案**：作为baseline与消融。
- **最终可用性**：可用。
- **来源链接**：https://openreview.net/forum?id=i_Q1yrOegLY

### 110. AdaMAE: Adaptive Masking for Efficient Spatiotemporal Learning with Masked Autoencoders

- **排序文件**：`01_PDF论文_按推荐程度排序/110_基线_背景_AdaMAE_Adaptive_Masking_for_Efficient_Spatiotemporal_Learning_with_Masked_Autoencoders.pdf`
- **推荐等级**：基线/背景（score=25）
- **SHA-256**：`C10CF0105691FC4E32C5CE995C5256823353EDA1F447EDEBDAE863BFED7ECFA6`
- **来源记录**：文库2#5 05_AdaMAE.pdf

#### 来源说明：文库2 #5

说明来源：文库2/scmae_improvement_reports/05_AdaMAE.txt
- **标题**：AdaMAE: Adaptive Masking for Efficient Spatiotemporal Learning with Masked Autoencoders
- **年份/会议**：2023 / CVPR
- **备注**：CVF开放PDF。

AdaMAE – Adaptive Masking for Efficient Spatiotemporal Learning

1. **Efficacy & principle** – AdaMAE introduces an auxiliary sampling network to estimate a probability distribution over tokens and selects informative tokens using policy gradient (REINFORCE).  This adaptive sampling network allows high masking ratios (up to 95%), reducing memory and training time while improving accuracy【537530205682633†L15-L33】【537530205682633†L88-L116】.
2. **Integration strategy** – For scMAE, build a token‑sampling network that assigns higher sampling probability to informative genes (e.g., highly variable or biologically important genes) and sample unmasked genes accordingly.  Use policy gradient to train the sampler jointly with the autoencoder.
3. **Potential issues** – The REINFORCE algorithm has high variance and may slow convergence.  Sampling network introduces additional parameters and complexity.  High masking may still remove critical information if sampling is inaccurate.
4. **Mitigation** – Use variance‑reduction techniques (e.g., baseline subtraction) in policy gradient.  Incorporate gene importance priors to warm start the sampling network.  Limit maximum masking ratio based on dataset sparsity.
5. **Usability decision** – Adaptive sampling can make scMAE more efficient and focus on informative genes.  With proper tuning, it is usable.

Citations:【537530205682633†L15-L33】【537530205682633†L88-L116】

### 111. scMAE: a masked autoencoder for single-cell RNA-seq clustering

- **排序文件**：`01_PDF论文_按推荐程度排序/111_基线_背景_scMAE_a_masked_autoencoder_for_single-cell_RNA-seq_clustering.pdf`
- **推荐等级**：基线/背景（score=25）
- **SHA-256**：`4563BBFD3DCA2A6168FA3050D9961C8D5AAE39FF94077A1CB45C8ABB8E8E7D45`
- **来源记录**：文库3#2 02_scMAE_Fang_2024_Bioinformatics.pdf

#### 来源说明：文库3 #2

说明来源：文库3/manifest_30_papers.csv
- **标题**：scMAE: a masked autoencoder for single-cell RNA-seq clustering
- **年份/会议**：Bioinformatics 2024
- **方向/领域**：Single-cell SSL
- **与 scMAE 改良相关性**：Target baseline: shuffled-value masking + mask prediction + weighted reconstruction + k-means/Leiden.

### 112. SC3: consensus clustering of single-cell RNA-seq data

- **排序文件**：`01_PDF论文_按推荐程度排序/112_基线_背景_SC3_consensus_clustering_of_single-cell_RNA-seq_data.pdf`
- **推荐等级**：基线/背景（score=25）
- **SHA-256**：`6342B5142DF5118EC23622DE2E5D28DA3DB2C9EFAF47ED11F897A2796EA18093`
- **来源记录**：文库3#17 12_SC3_Kiselev_2017_Consensus_clustering_of_single_cell_RNA_seq.pdf

#### 来源说明：文库3 #17

说明来源：文库3/manifest_30_papers.csv
- **标题**：SC3: consensus clustering of single-cell RNA-seq data
- **年份/会议**：Nature Methods 2017
- **方向/领域**：Single-cell consensus clustering
- **与 scMAE 改良相关性**：Consensus clustering and robustness baseline.

### 113. Fast unfolding of communities in large networks

- **排序文件**：`01_PDF论文_按推荐程度排序/113_基线_背景_Fast_unfolding_of_communities_in_large_networks.pdf`
- **推荐等级**：基线/背景（score=25）
- **SHA-256**：`7765151A1C6556A6E24623BCCED75F2C6AC04093B8E27BA8DF50BC8F2ECB5577`
- **来源记录**：文库3#18 13_Louvain_Blondel_2008_Fast_unfolding_of_communities_in_large_networks.pdf

#### 来源说明：文库3 #18

说明来源：文库3/manifest_30_papers.csv
- **标题**：Fast unfolding of communities in large networks
- **年份/会议**：J. Stat. Mech. 2008
- **方向/领域**：Community detection
- **与 scMAE 改良相关性**：Louvain baseline used in many scRNA workflows.

## 有建议但缺少PDF的条目

### 文库2 #50. Integrating feature selection with unsupervised deep embedding for clustering single-cell RNA-seq data

- **推荐等级**：可用（score=70）
- **说明来源**：scmae_improvement_reports/50_FSSC.txt

说明来源：文库2/scmae_improvement_reports/50_FSSC.txt
- **标题**：Integrating feature selection with unsupervised deep embedding for clustering single-cell RNA-seq data
- **年份/会议**：2026 / Briefings in Bioinformatics
- **备注**：OUP开放文章页/PMC；直接PDF链接未稳定定位。

FSSC – Feature Selection and Clustering

1. **Efficacy & principle** – FSSC integrates feature selection with unsupervised deep embedding by using a zero‑inflated negative binomial autoencoder, a group Lasso penalty and a clustering loss.  It learns low‑dimensional representations while selecting cluster‑discriminatory genes, outperforming existing methods【904685930220572†L114-L124】.
2. **Integration strategy** – Incorporate FSSC’s feature selection mechanism into scMAE: add a group Lasso penalty on gene weights in the encoder to promote sparsity and discriminative features.  Use ZINB reconstruction to model counts.
3. **Potential issues** – Lasso penalties may discard genes needed for downstream tasks.  Balancing feature selection and reconstruction may be challenging.
4. **Mitigation** – Use adaptive group Lasso where penalties depend on gene importance.  Monitor clustering quality while adjusting penalty strength.
5. **Usability decision** – FSSC’s joint feature selection and clustering can improve scMAE’s interpretability and efficiency; it is usable with careful tuning.

Citations:【904685930220572†L114-L124】

### 文库2 #43. scDFC: A deep fusion clustering method for single-cell RNA-seq data

- **推荐等级**：可用（score=70）
- **说明来源**：scmae_improvement_reports/43_scDFC.txt

说明来源：文库2/scmae_improvement_reports/43_scDFC.txt
- **标题**：scDFC: A deep fusion clustering method for single-cell RNA-seq data
- **年份/会议**：2023 / Briefings in Bioinformatics
- **备注**：OUP文章页；直接PDF链接未稳定定位。

scDFC – Deep Fusion Clustering

1. **Efficacy & principle** – scDFC proposes a deep fusion clustering model with two modules: an attributed feature clustering module and a structure‑attention feature clustering module.  Two autoencoders handle attribute and structural features, and experiments demonstrate that fusing attributes, structure and attention improves clustering on scRNA‑seq data【603915245475056†L165-L181】.
2. **Integration strategy** – Use scDFC’s dual autoencoder structure with scMAE: one autoencoder processes masked gene expressions, while another processes cell‑cell structural features.  Fuse the embeddings via attention and then perform clustering.
3. **Potential issues** – Requires learning two separate autoencoders; fusion parameters may be difficult to set.  Structural information may be noisy.
4. **Mitigation** – Pretrain each autoencoder and fine‑tune jointly.  Use regularization to balance attribute and structure contributions.  Denoise structural information with graph filtering.
5. **Usability decision** – scDFC’s fusion approach can enrich scMAE’s representations and is usable with careful tuning.

Citations:【603915245475056†L165-L181】

### 文库2 #41. scGANSL: Graph Attention Network with Subspace Learning for scRNA-seq Data Clustering

- **推荐等级**：可用（score=70）
- **说明来源**：scmae_improvement_reports/41_scGANSL.txt

说明来源：文库2/scmae_improvement_reports/41_scGANSL.txt
- **标题**：scGANSL: Graph Attention Network with Subspace Learning for scRNA-seq Data Clustering
- **年份/会议**：2025 / Journal of Chemical Information and Modeling / ACS
- **备注**：ACS JCIM；原文通常需机构权限，附摘要背景方法。

scGANSL – Graph Attention Network with Subspace Learning

1. **Efficacy & principle** – scGANSL constructs two views using highly variable genes and principal component analysis, feeds them into a multiview shared graph autoencoder and integrates a zero‑inflated negative binomial (ZINB) model into a self‑supervised graph attention autoencoder.  It introduces local learning and self‑expression strategies to preserve local and global structures and significantly outperforms other methods【296200983764465†L96-L110】.
2. **Integration strategy** – Use scGANSL’s multi‑view graph autoencoder to produce initial embeddings for scMAE.  The self‑expression strategy can guide mask selection by highlighting genes that capture both local and global structure.  The ZINB model ensures appropriate handling of dropout.
3. **Potential issues** – Multi‑view processing increases complexity.  Self‑expression requires solving a reconstruction problem, which can be expensive.  ZINB integration adds parameters.
4. **Mitigation** – Use dimensionality reduction to limit view size.  Employ sparse self‑expression and efficient solvers.  Share ZINB parameters across views.
5. **Usability decision** – scGANSL’s view integration and attention can enhance scMAE; it is usable with careful resource management.

Citations:【296200983764465†L96-L110】

### 文库2 #36. scVIC: deep generative modeling of heterogeneity for scRNA-seq data

- **推荐等级**：可用（score=70）
- **说明来源**：scmae_improvement_reports/36_scVIC.txt

说明来源：文库2/scmae_improvement_reports/36_scVIC.txt
- **标题**：scVIC: deep generative modeling of heterogeneity for scRNA-seq data
- **年份/会议**：2024 / Bioinformatics Advances
- **备注**：OUP文章页；直接PDF链接未稳定定位。

scVIC – Variational Inference for scRNA‑seq Clustering

1. **Efficacy & principle** – scVIC uses variational inference to explicitly model biological heterogeneity and technical variability (including batch effects).  It handles dropout events and batch effects while inferring cell heterogeneity, leading to improved clustering【582714882660786†L162-L173】.
2. **Integration strategy** – Incorporate scVIC’s probabilistic modelling into scMAE: use latent variables to capture batch effects and heterogeneity, and incorporate a batch‑adjusted reconstruction loss.  This could make scMAE robust to batch effects.
3. **Potential issues** – Modelling batch effects adds complexity; mis‑specifying technical noise distribution may hurt performance.
4. **Mitigation** – Use empirical Bayes to estimate batch effect parameters and apply domain adaptation techniques.  Validate with known batch labels.
5. **Usability decision** – Handling batch effects is critical; scVIC’s approach is usable within scMAE for robust clustering.

Citations:【582714882660786†L162-L173】

### 文库2 #35. scMVAF: a multi-view adaptive fusion clustering approach for single-cell RNA-sequencing data

- **推荐等级**：可用（score=70）
- **说明来源**：scmae_improvement_reports/35_scMVAF.txt

说明来源：文库2/scmae_improvement_reports/35_scMVAF.txt
- **标题**：scMVAF: a multi-view adaptive fusion clustering approach for single-cell RNA-sequencing data
- **年份/会议**：2026 / Briefings in Bioinformatics
- **备注**：OUP开放PDF。

scMVAF – Multi‑View Adaptive Fusion Clustering

1. **Efficacy & principle** – scMVAF down‑samples scRNA‑seq data to generate multiple views, denoises each view using a ZINB‑based model, and uses adaptive weighted fusion to integrate embeddings.  It iteratively refines clustering with pseudo‑labels and outperforms eight clustering methods across 16 datasets【937305929821833†L110-L147】.
2. **Integration strategy** – Use scMVAF to produce denoised embeddings for scMAE: create multiple views, denoise with ZINB, fuse them adaptively and input into the masked autoencoder.  Use pseudo‑labels to refine the model iteratively.
3. **Potential issues** – Down‑sampling may discard important cells; pseudo‑labeling can propagate errors.  Fusion weight estimation can be unstable.
4. **Mitigation** – Ensure that down‑sampling preserves rare cell types.  Use robust pseudo‑labeling and gradually update labels.  Regularize fusion weights.
5. **Usability decision** – scMVAF’s adaptive fusion provides strong representations; it is usable as a preprocessing step for scMAE.

Citations:【937305929821833†L110-L147】

### 文库2 #34. Multi-view clustering for single-cell RNA-seq data based on graph fusion

- **推荐等级**：可用（score=70）
- **说明来源**：scmae_improvement_reports/34_scMCGF.txt

说明来源：文库2/scmae_improvement_reports/34_scMCGF.txt
- **标题**：Multi-view clustering for single-cell RNA-seq data based on graph fusion
- **年份/会议**：2025 / Briefings in Bioinformatics
- **备注**：OUP文章页；直接PDF链接未稳定定位。

scMCGF – Multi‑View Graph Fusion Clustering

1. **Efficacy & principle** – scMCGF constructs multiple data views (gene expression, pathway scores, PCA, diffusion maps), learns a similarity graph for each and fuses them into a unified graph using adaptive weighted fusion.  The fused graph is used for clustering and improves accuracy and stability across datasets【55779319394859†L93-L117】.
2. **Integration strategy** – Use scMCGF’s multi‑view fusion to build an informative cell graph for scMAE.  Incorporate pathway scores and diffusion maps as auxiliary inputs, fuse graphs adaptively and then perform masked autoencoding on the fused graph representation.
3. **Potential issues** – Generating multiple views and fusing graphs may be computationally expensive.  Weight selection can be sensitive.
4. **Mitigation** – Select a few informative views and use efficient graph fusion algorithms.  Learn weights via backpropagation instead of manual tuning.
5. **Usability decision** – Multi‑view fusion captures complementary information; it is usable as a preprocessing step for scMAE.

Citations:【55779319394859†L93-L117】

### 文库2 #33. ScCCL: Single-Cell Data Clustering Based on Self-Supervised Contrastive Learning

- **推荐等级**：可用（score=70）
- **说明来源**：scmae_improvement_reports/33_scCCL.txt

说明来源：文库2/scmae_improvement_reports/33_scCCL.txt
- **标题**：ScCCL: Single-Cell Data Clustering Based on Self-Supervised Contrastive Learning
- **年份/会议**：2023 / IEEE/ACM TCBB
- **备注**：IEEE/ACM TCBB；原文通常需机构权限，附摘要背景方法。

ScCCL – Self‑Supervised Contrastive Clustering

1. **Efficacy & principle** – ScCCL randomly masks gene expressions twice, adds Gaussian noise and employs a momentum encoder to extract features.  Instance‑level and cluster‑level contrastive learning modules yield high‑order embeddings, improving clustering accuracy without relying on graph structure【615012199035452†L262-L268】.
2. **Integration strategy** – Apply ScCCL pretraining before scMAE: generate two noisy views of each cell, train a momentum encoder with contrastive loss, then fine‑tune scMAE on these representations.  Use cluster‑level contrastive loss to encourage grouping.
3. **Potential issues** – Random masking may remove critical genes; Gaussian noise may distort gene expression.  Momentum encoder requires memory of past states.
4. **Mitigation** – Use biologically informed masking and noise levels.  Employ smaller momentum coefficients to adapt quickly.
5. **Usability decision** – ScCCL enhances representation learning; it is usable as a pretraining strategy for scMAE.

Citations:【615012199035452†L262-L268】

### 文库2 #32. Decoupled GNNs based on multi-view contrastive learning for scRNA-seq data clustering

- **推荐等级**：可用（score=70）
- **说明来源**：scmae_improvement_reports/32_scDeGNN.txt

说明来源：文库2/scmae_improvement_reports/32_scDeGNN.txt
- **标题**：Decoupled GNNs based on multi-view contrastive learning for scRNA-seq data clustering
- **年份/会议**：2025 / Briefings in Bioinformatics
- **备注**：OUP文章页；直接PDF链接未稳定定位。

scDeGNN – Decoupled Graph Neural Network for scRNA‑seq

1. **Efficacy & principle** – scDeGNN constructs two adjacency matrices to create distinct graph views, trains decoupled GNNs for each view and refines representations via a multilayer perceptron and contrastive learning layer.  It fuses representations for clustering and addresses high computational complexity【307201161404088†L162-L170】.
2. **Integration strategy** – Use scDeGNN’s multi‑view GNN as a feature extractor for scMAE: obtain two complementary graph embeddings and concatenate them with masked autoencoder features.  Employ contrastive loss to align views.
3. **Potential issues** – Building multiple graphs increases computational and memory requirements.  Fusing views may require complex weighting.
4. **Mitigation** – Use sparse graphs and sampling strategies to reduce cost.  Learn fusion weights adaptively via attention mechanisms.
5. **Usability decision** – Multi‑view graph embedding improves representation diversity; it is usable when resources permit.

Citations:【307201161404088†L162-L170】

### 文库2 #27. scSCC: A swapped contrastive learning-based clustering method for single-cell gene expression data

- **推荐等级**：可用（score=70）
- **说明来源**：scmae_improvement_reports/27_scSCC.txt

说明来源：文库2/scmae_improvement_reports/27_scSCC.txt
- **标题**：scSCC: A swapped contrastive learning-based clustering method for single-cell gene expression data
- **年份/会议**：2025 / Quantitative Biology
- **备注**：文章页开放；直接PDF链接未稳定定位。

scSCC – Swapped Contrastive Clustering

1. **Efficacy & principle** – scSCC combines instance contrastive learning with a swapped prediction module using prototypes.  The swapped prediction encourages cells of the same cluster to converge to shared prototypes, producing clustering‑friendly representations【42213496449950†L125-L146】.
2. **Integration strategy** – Add contrastive pretraining to scMAE: augment cells and apply instance contrastive loss.  Introduce a prototype head that encourages masked embeddings to align with cluster centroids, improving clustering downstream.
3. **Potential issues** – Contrastive learning needs careful design of augmentations; inappropriate augmentations may harm gene semantics.  Prototype assignments may be unstable when clusters are unknown.
4. **Mitigation** – Use biologically meaningful augmentations (e.g., gene masking, dropout simulation).  Update prototypes using moving averages and use pseudo‑labels iteratively.
5. **Usability decision** – Contrastive learning can enhance latent representations.  scSCC’s swapped prediction mechanism is usable when adapted to scRNA augmentations.

Citations:【42213496449950†L125-L146】

### 文库2 #26. scGAAC: A graph attention autoencoder for clustering single-cell RNA-sequencing data

- **推荐等级**：可用（score=70）
- **说明来源**：scmae_improvement_reports/26_scGAAC.txt

说明来源：文库2/scmae_improvement_reports/26_scGAAC.txt
- **标题**：scGAAC: A graph attention autoencoder for clustering single-cell RNA-sequencing data
- **年份/会议**：2024 / Methods / Elsevier
- **备注**：Elsevier Methods；本环境未找到开放PDF，附摘要背景方法。

scGAAC – Graph Attention Autoencoder Clustering

1. **Efficacy & principle** – scGAAC uses a graph attention autoencoder that considers both gene expression and relationships between cells.  An attention fusion module combines features from a graph attention autoencoder and a traditional autoencoder, and self‑supervised learning optimizes clustering.  It outperforms state‑of‑the‑art methods【836303828812709†L125-L147】.
2. **Integration strategy** – For scMAE, embed a graph attention autoencoder as a pre‑encoder to capture cell–cell relationships.  Use the attention fusion mechanism to merge graph‑based features with masked autoencoder features, improving representation quality.
3. **Potential issues** – Requires careful balancing between graph and gene features.  Self‑supervised loss may conflict with reconstruction objective.
4. **Mitigation** – Use multi‑task learning to balance losses.  Normalize and weight graph and gene features to prevent dominance of either.
5. **Usability decision** – scGAAC’s fusion of graph and autoencoder is directly applicable to scMAE and is usable.

Citations:【836303828812709†L125-L147】

### 文库2 #25. scGAC: a graph attentional architecture for clustering single-cell RNA-seq data

- **推荐等级**：可用（score=70）
- **说明来源**：scmae_improvement_reports/25_scGAC.txt

说明来源：文库2/scmae_improvement_reports/25_scGAC.txt
- **标题**：scGAC: a graph attentional architecture for clustering single-cell RNA-seq data
- **年份/会议**：2022 / Bioinformatics
- **备注**：OUP文章页；直接PDF链接未稳定定位。

scGAC – Graph Attentional Clustering

1. **Efficacy & principle** – scGAC constructs a cell graph, refines it via network denoising, and uses a graph attentional autoencoder to learn clustering‑friendly representations.  A self‑optimizing method then produces cell clusters.  Experiments on 16 datasets show excellent performance【701588048509262†L96-L104】.
2. **Integration strategy** – Integrate scGAC’s graph attentional encoder into scMAE.  Use the refined cell graph to inform masking: mask genes in clusters of similar cells and reconstruct them jointly.  The attention mechanism helps weight neighbours.
3. **Potential issues** – Requires graph construction and denoising; if the graph is inaccurate, attention propagation may amplify errors.  Self‑optimizing clustering may not scale.
4. **Mitigation** – Use robust graph construction (e.g., mutual k‑nearest neighbours) and incorporate graph regularization.  Combine self‑optimizing clustering with mini‑batch training.
5. **Usability decision** – scGAC’s attention network can enhance scMAE’s embedding; it is usable with robust graph design.

Citations:【701588048509262†L96-L104】

### 文库2 #10. Variational Masked AutoEncoder topic / related work cluster

- **推荐等级**：可用（score=70）
- **说明来源**：scmae_improvement_reports/10_VMAE.txt

说明来源：文库2/scmae_improvement_reports/10_VMAE.txt
- **标题**：Variational Masked AutoEncoder topic / related work cluster
- **年份/会议**：2026 / topic/various
- **备注**：该条是方法族/主题页，不是单篇原始论文。

VMAE – Variational Masked Autoencoder

1. **Efficacy & principle** – VMAE incorporates variational inference into masked autoencoders, creating a structured latent manifold.  It fits stochastic latent variables for each mask and uses heavy masking to encourage meaningful abstraction【890284605854668†L40-L59】.
2. **Integration strategy** – For scMAE, introduce a variational bottleneck with prior distributions over the latent representation.  The model learns probabilistic embeddings of cells and can quantify uncertainty.  Variational training may regularize the model and improve generalization.
3. **Potential issues** – Variational inference increases complexity and may lead to posterior collapse if the decoder is too flexible.  The KL divergence term needs careful weighting.
4. **Mitigation** – Use annealing schedules or mutual information regularization to prevent posterior collapse.  Employ flexible priors (e.g., mixtures of Gaussians) appropriate for scRNA‑seq data.
5. **Usability decision** – Introducing a variational component can improve scMAE’s robustness and interpretability; it is usable with proper regularization.

Citations:【890284605854668†L40-L59】

### 文库2 #13. scVAE: variational auto-encoders for single-cell gene expression data

- **推荐等级**：未明确（score=50）
- **说明来源**：scmae_improvement_reports/13_scVAE.txt

说明来源：文库2/scmae_improvement_reports/13_scVAE.txt
- **标题**：scVAE: variational auto-encoders for single-cell gene expression data
- **年份/会议**：2020 / Bioinformatics
- **备注**：OUP开放文章页；直接PDF链接可能需站点会话。

scVAE – Variational Autoencoder for scRNA‑seq

1. **Efficacy & principle** – scVAE uses a variational autoencoder that takes raw count data as input and models gene expression using appropriate likelihoods (e.g., negative binomial).  It learns latent representations for each cell and outperforms existing methods in clustering【597422546486935†L163-L179】.
2. **Integration strategy** – Integrate scVAE’s negative‑binomial likelihood and latent modelling into scMAE’s reconstruction loss.  Use variational inference for latent cell embeddings while employing masking to select genes.
3. **Potential issues** – Variational models can suffer from posterior collapse or require careful tuning of priors.  Additionally, combining MAE and VAE objectives may complicate training.
4. **Mitigation** – Use annealing and mutual information regularization to maintain latent information.  Balance reconstruction and KL terms carefully.
5. **Usability decision** – scVAE provides strong modelling of count distributions.  Integrating its ideas into scMAE is feasible and beneficial.

Citations:【597422546486935†L163-L179】

### 文库2 #1. scMAE: a masked autoencoder for single-cell RNA-seq clustering

- **推荐等级**：基线/背景（score=25）
- **说明来源**：scmae_improvement_reports/01_scMAE_baseline.txt

说明来源：文库2/scmae_improvement_reports/01_scMAE_baseline.txt
- **标题**：scMAE: a masked autoencoder for single-cell RNA-seq clustering
- **年份/会议**：2024 / Bioinformatics
- **备注**：OUP开放文章页；本环境未能稳定定位直接PDF URL。

scMAE baseline

1. **Efficacy & principle** – scMAE is a masked autoencoder for single‑cell RNA‑seq clustering.  It perturbs gene expressions and uses a masking predictor to reconstruct the original data, learning latent relationships among genes.  The method captures cell‑level structure by predicting whether gene expression values were masked, improving clustering【828124119550327†L118-L125】.
2. **Integration strategy** – scMAE serves as the starting point.  It introduces a masking predictor capturing gene relationships; subsequent methods can improve on the masking strategy, decoder design or latent regularisation.  One could modify the masking pattern or incorporate domain knowledge (e.g., gene networks) while keeping the reconstruction objective.
3. **Potential issues** – High masking ratios may degrade reconstruction if the model is not trained long enough.  The method lacks explicit biological constraints and may treat all genes equally, ignoring known gene interactions.
4. **Mitigation** – Integrate prior knowledge (PPI networks or gene ontology) into the encoder, or use adaptive masking to select informative genes.  Incorporate decoders that account for zero‑inflated distributions and noise.
5. **Usability decision** – scMAE forms a strong foundation; modifications are needed for better performance.  It remains usable as a baseline to compare improvements.

Citations:【828124119550327†L118-L125】

## 文库2无法直接下载或不确定条目原文

# 无法稳定获得原文 PDF / 非单篇原始论文的摘要、背景、思路和方法

## 01. scMAE baseline

摘要：scMAE 是面向 scRNA-seq 聚类的 masked autoencoder，通过扰动/遮蔽基因表达并重建原始表达来学习细胞表示。背景：scRNA-seq 数据高维、稀疏、dropout 严重，传统聚类和普通 autoencoder 难以充分利用基因相关性。思路：把基因遮蔽作为自监督任务，引入 mask predictor 捕捉基因间依赖。方法：输入预处理后的表达矩阵，随机扰动/遮蔽部分基因，编码器学习低维表示，解码器重建表达，并用预测遮蔽状态的辅助任务增强基因相关性建模。

## 10. VMAE concept

摘要：VMAE 不是单篇固定论文，而是把 masked autoencoder 与 variational inference 结合的方法族。背景：标准 MAE 给出点估计潜变量，难以表达不确定性；scRNA-seq 中 dropout 与技术噪声使不确定性建模更重要。思路：在遮蔽重建任务中引入随机潜变量和 KL 正则，使潜空间更平滑。方法：编码可见基因得到后验分布 q(z|x_v)，采样 z 后重建 masked genes，并加入 KL(q||p)；可进一步用 ZINB/NB likelihood 适配计数数据。

## 13. scVAE

摘要：scVAE 用变分自编码器直接建模单细胞基因表达计数，学习细胞潜表示并可做聚类。背景：scRNA-seq 数据需要同时处理离散计数、dropout、批次差异和高维噪声。思路：用概率生成模型替代单纯 MSE 重建，使潜空间保留细胞群结构。方法：以原始计数为输入，使用 VAE 编码细胞潜变量，解码时用合适的计数 likelihood（如 NB/ZINB）拟合基因表达，并可用混合先验得到聚类结构。

## 25. scGAC

摘要：scGAC 是图注意力聚类方法，构建并去噪细胞图，再通过图注意力自编码器学习聚类友好的表示。背景：普通 autoencoder 只看单个细胞表达，忽略细胞间关系；scRNA-seq 的稀疏性会使简单相似度图不可靠。思路：把细胞作为图节点，用注意力机制为邻居分配不同权重，捕捉潜在细胞关系。方法：先构建 cell graph 并做网络去噪，再用 graph attentional autoencoder 传播邻域信息，最后用 self-optimizing clustering 得到细胞簇。

## 26. scGAAC

摘要：scGAAC 是 Methods 2024 上的图注意力自编码器聚类方法，PubMed 记录显示 DOI 为 10.1016/j.ymeth.2024.06.010。背景：许多 scRNA-seq 聚类方法只提取单细胞表达特征，忽略细胞间关系。思路：利用图注意力自编码器同时建模基因表达和细胞相似图。方法：构建细胞图，使用图注意力编码器学习邻居加权表示，再通过自编码/聚类模块获得细胞簇；其核心价值在于把关系建模与表达特征学习统一起来。

## 27. scSCC

摘要：scSCC 是基于 swapped contrastive learning 的 scRNA-seq 聚类方法，结合 instance contrastive learning 与 swapped prediction。背景：对比学习能缓解稀疏表达的重建依赖，但普通实例级对比不一定显式增强簇结构。思路：通过聚类原型把簇信号注入潜空间，使同簇细胞靠近同一 prototype，不同簇远离。方法：生成增强视图，执行实例级对比；同时进行 swapped prediction，即用一个视图的原型分配监督另一个视图的预测分布，从而学习更利于聚类的表示。

## 32. scDeGNN

摘要：scDeGNN 使用多视图对比学习和 decoupled GNN 改善 scRNA-seq 聚类。背景：GNN 层数增加会使依赖复杂度迅速升高，训练效率下降，并可能出现过平滑。思路：构建两个邻接矩阵作为不同视图，用解耦 GNN 分别学习，再通过对比学习保持一致性与判别性。方法：先生成多视图邻接，分别用 decoupled GNN 得到嵌入，再经过 MLP 与 contrastive layer 精炼表示，最后融合嵌入做聚类。

## 33. ScCCL

摘要：ScCCL 是 IEEE/ACM TCBB 2023 的自监督对比聚类方法，DOI 为 10.1109/TCBB.2023.3241129。背景：scRNA-seq 高维且稀疏，直接重建原始表达容易被噪声主导。思路：构造两个增强视图，通过动量编码器和对比损失学习稳定表示，并在实例级与簇级同时优化。方法：对每个细胞随机遮蔽基因表达并加入高斯噪声形成两个视图，用 momentum encoder 提取特征；随后使用 instance-level contrastive loss 和 cluster-level contrastive loss 共同训练。

## 34. scMCGF

摘要：scMCGF 是基于多视图图融合的 scRNA-seq 聚类方法。背景：单一视角往往无法充分描述细胞关系；基因表达、通路活性、PCA 与 diffusion map 能分别捕捉不同结构。思路：构造多个视图的相似图，并自适应融合为统一图。方法：对表达矩阵预处理后计算通路分数，并提取 PCA/diffusion map 特征；每个视图学习相似图，最终用加权融合和图拉普拉斯约束得到聚类。

## 36. scVIC

摘要：scVIC 是面向 scRNA-seq 异质性建模的深度生成模型，强调同时处理生物异质性、dropout 与批次效应。背景：批次效应和技术噪声会掩盖真实生物异质性。思路：用 variational inference 显式区分生物变异和技术变异。方法：构建概率生成模型，对细胞潜变量、批次效应和 dropout/噪声参数进行联合推断；训练后用去噪且去批次的潜表示进行聚类。

## 41. scGANSL

摘要：scGANSL 是 ACS JCIM 2025 的图注意力网络与子空间学习结合方法，DOI 为 10.1021/acs.jcim.5c00731。背景：多数方法只使用单视图，无法完整解释高维稀疏 scRNA-seq 数据。思路：用 HVG 与 PCA 构建双视图，结合多视图共享图自编码器、ZINB 模型和自表达子空间学习。方法：两个视图分别进入共享图自编码器；聚类标签指导潜表示和系数矩阵学习；ZINB 捕捉 dropout 和过度离散；局部学习与 self-expression 共同约束潜空间局部/全局结构。

## 43. scDFC

摘要：scDFC 是深度融合聚类模型，包含属性特征聚类模块和结构-注意力特征聚类模块。背景：单细胞数据既有表达属性信息，也有细胞间结构信息，单独建模会遗漏一部分信号。思路：分别处理表达属性和结构注意力信息，再融合以获得更完整的表示。方法：构建两个自编码器：一个处理表达属性，一个处理细胞图结构与注意力特征；融合后用于聚类优化。

## 50. FSSC

摘要：FSSC 将基因特征选择与深度嵌入聚类联合优化。背景：先做 HVG/特征选择再聚类可能与最终聚类结构不一致，导致选出的基因不一定最能区分细胞群。思路：在一个统一模型中同时学习低维表示、选择基因和优化聚类目标。方法：ZINB autoencoder 拟合 count/dropout 特征，group Lasso 促使基因层面稀疏选择，专门的 clustering loss 使选中基因具有簇判别性。
