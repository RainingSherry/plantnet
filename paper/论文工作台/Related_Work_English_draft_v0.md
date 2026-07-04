# Related Work English Draft v0

Updated: 2026-06-26

## 1. Clustering and Benchmarking in scRNA-seq Analysis

Cell clustering is a foundational step in scRNA-seq analysis because it provides the basis for cell-type discovery, marker gene identification, atlas construction, and downstream biological interpretation. The difficulty of this task arises from the high dimensionality, sparsity, noise, and biological heterogeneity of single-cell expression matrices. These challenges have been discussed extensively in surveys and methodological reviews of unsupervised single-cell clustering \citep{kiselev2019challenges}. More recently, scCluBench provided a systematic benchmark covering diverse datasets and methodological families, including traditional clustering, deep learning methods, graph-based methods, and biological foundation models \citep{xu2026scclubench}. This benchmark perspective is important for the present work: a new clustering method should not be evaluated only on a small number of favorable datasets, but should be tested under standardized protocols that include clustering metrics, biological annotation, representation quality, and scalability.

## 2. Classical and Graph-Based Clustering Methods

Classical single-cell clustering pipelines often combine preprocessing, dimensionality reduction, and general-purpose clustering algorithms. SC3 performs consensus clustering and remains an influential reference for unsupervised cell-type discovery \citep{kiselev2017sc3}. Community detection algorithms such as Louvain and Leiden are widely used in graph-based single-cell analysis pipelines \citep{blondel2008louvain,traag2019leiden}, and Seurat provides a widely adopted ecosystem for single-cell integration, preprocessing, and graph clustering \citep{stuart2019seurat}. These methods are robust, interpretable, and computationally attractive, but their performance depends heavily on preprocessing choices, distance metrics, and graph construction. In highly sparse or noisy datasets, low-dimensional neighborhood graphs may fail to preserve subtle cell-state boundaries or rare populations.

Graph neural network methods attempt to improve this limitation by learning representations over cell-cell or gene-cell graphs. scGNN models single-cell data with graph neural networks \citep{wang2021scgnn}, scDSC jointly uses an autoencoder and graph neural network for structural clustering \citep{gan2022scdsc}, AttentionAE-sc introduces attention into deep clustering for scRNA-seq \citep{li2023attentionae}, and scCDCG uses cut-informed graph embedding to improve clustering structure \citep{xu2024sccdcg}. These methods explicitly encode relational information, but they also inherit sensitivity to graph construction and may suffer from over-smoothing or representation collapse. CAAM-scMAE is motivated by a complementary direction: instead of constructing a fixed graph, it injects population-level context through cell-axis attention during masked reconstruction.

## 3. Deep Generative and Autoencoder-Based Representation Learning

Autoencoder-based models have become a major family of representation learning methods for scRNA-seq data. General deep embedded clustering methods such as DEC \citep{xie2016dec} inspired later single-cell adaptations that combine latent representation learning with clustering objectives. DCA introduced a deep count autoencoder for denoising scRNA-seq data under count-based assumptions \citep{eraslan2019dca}, while scVI used deep generative modeling to learn probabilistic latent representations for single-cell transcriptomics \citep{lopez2018scvi}. The broader scvi-tools ecosystem further standardized deep probabilistic analysis for single-cell omics \citep{gayoso2022scvitools}.

Several models adapt deep clustering more directly to single-cell analysis. scDeepCluster combines autoencoding with a model-based clustering formulation \citep{tian2019scdeepcluster}. DESC learns cluster-friendly embeddings while addressing batch effects \citep{li2020desc}. scDCC incorporates constrained clustering into deep embedding learning \citep{tian2021scdcc}, and scziDesk uses deep soft K-means with self-training \citep{chen2020sczidesk}. These methods show that task-specific objectives can improve clustering performance, but many of them rely on latent clustering assumptions, iterative pseudo-labeling, or graph/cluster objectives that may become unstable when cell populations are rare, overlapping, or strongly affected by sparsity.

## 4. Masked and Contrastive Learning for Single-Cell Clustering

Self-supervised learning has recently become an important alternative for single-cell representation learning. Contrastive learning methods construct positive and negative pairs through data augmentation or neighborhood assumptions. scNAME combines neighborhood contrastive clustering with ancillary mask estimation \citep{wan2022scname}, showing that masking-related objectives can improve clustering-oriented representation learning. However, contrastive learning can be sensitive to augmentation design and may incorrectly push apart cells from the same biological population when negative pairs are sampled naively.

Masked modeling has been effective across domains, including language modeling with BERT \citep{devlin2019bert} and vision representation learning with masked autoencoders \citep{he2022mae}. Masked autoencoding offers a related self-supervised route for scRNA-seq data. scMAE perturbs expression values through gene-wise shuffling and jointly optimizes expression reconstruction and mask prediction \citep{fang2024scmae}. Its success suggests that learning to recover corrupted expression patterns can help the encoder capture gene dependencies and cell-type structure. Nevertheless, the standard scMAE design has two limitations that motivate CAAM-scMAE. First, a multilayer perceptron encoder does not explicitly model the two-dimensional organization of the expression matrix. Second, random masking does not distinguish between highly informative entries and weakly informative sparse entries. CAAM-scMAE addresses these limitations by combining bi-axial context encoding with constrained mask selection.

The mask selector in CAAM-scMAE is related to adversarial training only in the restricted sense that it selects harder reconstruction targets. It is not a free-form generative adversarial network. Classical GANs define a generator-discriminator game \citep{goodfellow2014gan}, and later work such as WGAN-GP shows that adversarial objectives require explicit stability constraints \citep{gulrajani2017wgangp}. This motivates our constrained formulation: the selector chooses mask positions under fixed budget, label-free input, coverage regularization, and sparsity-aware diagnostics, while replacement values remain tied to observed expression distributions.

## 5. Foundation Models and the Need for Task-Specific Clustering Objectives

Large-scale foundation models have recently expanded the scope of single-cell representation learning. Geneformer demonstrated transfer learning for network biology \citep{theodoris2023geneformer}, scGPT proposed a generative foundation model for single-cell multi-omics \citep{cui2024scgpt}, scFoundation trained a large-scale foundation model on single-cell transcriptomics \citep{hao2024scfoundation}, and CellFM further scaled pretraining to transcriptomes from 100 million human cells \citep{zeng2025cellfm}. Recent surveys summarize this rapidly expanding area \citep{zhang2025singlecellsurvey,dip2025llm4cell}. These models provide valuable transferable representations and enable broad downstream applications. However, clustering is a task with specific requirements: embeddings must separate biologically meaningful populations, preserve rare cell types, and avoid representation collapse. Benchmark evidence indicates that general-purpose foundation embeddings may perform well in classification or transfer tasks while still requiring task-specific adaptation for clustering \citep{xu2026scclubench,liu2026sceval}.

The architecture of TabPFN also motivates our design, but in a different way. TabPFN shows that tabular data can benefit from modeling both feature-wise and sample-wise interactions through two-dimensional attention \citep{hollmann2025tabpfn}. CAAM-scMAE does not use TabPFN as a pretrained model; instead, it borrows the structural principle that a table entry should be interpreted through both row and column contexts. In scRNA-seq, this corresponds to modeling both within-cell gene relationships and across-cell population context for the same gene or gene module.

## 6. Positioning of CAAM-scMAE

CAAM-scMAE is positioned between masked autoencoding, graph-based clustering, and single-cell foundation modeling. Compared with classical and graph-based clustering, it avoids committing to a fixed cell-cell graph before representation learning. Compared with deep clustering methods, it does not require labels or pseudo-labels during pretraining. Compared with scMAE, it adds explicit bi-axial context and learns a constrained mask policy rather than relying only on random corruption. Compared with general-purpose foundation models, it remains task-specific: its training objective is designed to produce embeddings suitable for unsupervised clustering.

This positioning leads to a concrete hypothesis: clustering-oriented masked reconstruction should be most effective when the model can use both gene-axis and cell-axis context, and when the masked positions are informative enough to provide nontrivial reconstruction signal without introducing shortcut artifacts. The empirical sections of this work should therefore test not only overall clustering accuracy, but also embedding distinguishability, rare cell recovery, mask behavior, runtime, and memory.

## Citation Keys Used

```text
kiselev2019challenges
xu2026scclubench
kiselev2017sc3
blondel2008louvain
traag2019leiden
stuart2019seurat
wang2021scgnn
gan2022scdsc
li2023attentionae
xu2024sccdcg
xie2016dec
eraslan2019dca
lopez2018scvi
gayoso2022scvitools
tian2019scdeepcluster
li2020desc
tian2021scdcc
chen2020sczidesk
wan2022scname
fang2024scmae
theodoris2023geneformer
cui2024scgpt
hao2024scfoundation
hollmann2025tabpfn
devlin2019bert
he2022mae
goodfellow2014gan
gulrajani2017wgangp
zeng2025cellfm
zhang2025singlecellsurvey
dip2025llm4cell
liu2026sceval
```
