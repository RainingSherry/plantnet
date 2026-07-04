# Introduction 中文草稿 v0

更新时间：2026-06-26

## 第一段：任务背景

单细胞 RNA 测序（single-cell RNA sequencing, scRNA-seq）能够在单细胞分辨率下刻画基因表达状态，已经成为解析细胞异质性、发现新细胞类型、构建细胞图谱和研究疾病微环境的重要技术。细胞聚类是 scRNA-seq 分析流程中的核心步骤之一，其目标是在没有人工标注的情况下，根据表达模式将细胞划分为具有生物学意义的亚群。高质量的聚类结果不仅影响后续的细胞类型注释和 marker gene 发现，也会进一步影响发育轨迹推断、疾病亚型识别和药物响应分析。

## 第二段：数据挑战

然而，scRNA-seq 数据具有高维、高稀疏、高噪声和强批次效应等特点，使得细胞聚类仍然具有挑战性。表达矩阵通常包含成千上万个基因，而每个细胞中只有一部分基因被观测到非零表达，大量零值可能同时来自真实生物学沉默和技术性 dropout。传统聚类方法依赖低维距离度量，难以捕获复杂的非线性表达结构；图聚类方法能够利用细胞间关系，但其性能高度依赖图构建策略，并可能出现过平滑或表示坍缩；通用单细胞 foundation model 虽然具有迁移能力，但其表示学习目标通常不是专门为 clustering 优化，因此在聚类任务上未必优于任务专用模型。

## 第三段：scMAE 的启发与不足

近期，masked autoencoder 为单细胞聚类提供了一条有前景的自监督路线。scMAE 通过 gene-wise shuffle corruption 扰动基因表达，并联合训练表达重建和 mask prediction，使 encoder 学习能够恢复被破坏表达值的细胞表示。已有结果表明，scMAE 在多个真实 scRNA-seq 数据集上取得了较好的聚类性能，并能识别低丰度 rare cell types。这说明 masked reconstruction 不仅可以作为去噪任务，也可以促使模型学习基因间依赖关系和细胞类型相关结构。

尽管如此，现有 masked autoencoder 仍存在两个关键局限。第一，scMAE 的 encoder 主要以 MLP 形式处理每个细胞的整行表达向量，缺少对表达矩阵二维结构的显式建模。scRNA-seq 表达矩阵天然包含两类上下文：同一细胞内不同基因之间的关系，以及同一基因在细胞群体中的分布差异。仅使用 MLP 压缩整行表达，可能无法充分利用跨细胞的 population-level context。第二，随机 mask 策略并不区分被 mask 位置的信息价值。在高稀疏矩阵中，大量 mask 可能落在零值或易恢复位置，使训练信号变弱；相反，一些与细胞类型区分、稀有细胞识别或关键通路相关的位置可能更值得被扰动和重建。

## 第四段：TabPFN 与双轴上下文

TabPFN 在表格数据建模中证明了二维结构的重要性。与将表格简单展开为一维序列不同，TabPFN 为表格中的每个 cell 构造表示，并通过 two-way attention 分别建模行内特征关系和列内样本关系。这一思想对 scRNA-seq 尤其有启发意义：表达矩阵同样可以被视为细胞和基因构成的二维表格，其中行对应细胞，列对应基因。因此，单细胞聚类模型不应只学习

$$
p(x_{ij}\mid X_{i,-j}),
$$

还应尽可能利用同一基因在其他细胞中的表达上下文：

$$
p(x_{ij}\mid X_{i,-j},X_{-i,j}).
$$

这种双轴上下文建模有望在不显式构建细胞图的情况下捕获 population structure，从而缓解图方法对邻接关系构建的依赖。

## 第五段：本文方法

基于上述观察，本文提出 CAAM-scMAE，一种面向 scRNA-seq clustering 的 context-aware adversarial axial masked autoencoder。CAAM-scMAE 以 scMAE 的 masked reconstruction 框架为基础，进一步引入两个关键模块。首先，本文设计轻量级 bi-axial context encoder，通过 gene-axis module attention 建模细胞内基因模块关系，并通过 cell-axis context attention 引入跨细胞上下文。其次，本文提出 constrained adversarial mask selector，在固定 mask budget、无标签泄漏、覆盖度受控和稀疏性约束下，学习选择更有信息量的 mask 位置，而不是完全依赖随机采样。

训练过程中，模型输入扰动后的表达矩阵，并同时优化表达重建损失与 mask prediction 损失；推理阶段仅使用 encoder 输出的 cell embedding 进行聚类。通过这种设计，CAAM-scMAE 试图将 masked reconstruction 从随机扰动任务转化为更贴近聚类需求的上下文恢复任务。

## 第六段：贡献

本文的主要贡献如下：

1. 提出一种面向 scRNA-seq 聚类的 bi-axial masked autoencoder，将 TabPFN 启发的二维上下文建模引入单细胞 masked representation learning。
2. 提出 constrained adversarial mask selector，在合法扰动集合内学习更有训练价值的 mask 位置，避免无约束生成式扰动带来的 shortcut 和训练不稳定。
3. 在 scCluBench-style 多数据集协议下系统评估方法性能，包括 ACC、NMI、ARI、rare cell type discovery、embedding distinguishability、mask diagnostics、runtime 和 memory。
4. 通过消融实验分离 encoder 结构、mask 策略和 reconstruction objective 的贡献，为单细胞聚类中的 masked self-supervised learning 提供可解释证据。

## 需要后续补强的位置

1. 每段都需要替换为正式英文表达。
2. 第二段需要插入 scCluBench 对高稀疏、大规模和 OOM/NaN 的证据。
3. 第三段需要插入 scMAE 的 Bioinformatics 引用。
4. 第四段需要插入 TabPFN 的 Nature 引用，并避免让审稿人误解为直接使用 TabPFN 预训练模型。
5. 第五段需要等实验确定后再决定是否保留 adversarial 字样。如果 adversarial mask 实验不稳定，可改为 informative mask selector 或 curriculum mask selector。
