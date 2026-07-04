# Introduction English Draft v0

Updated: 2026-06-26

> Note: This draft is written as a working manuscript section. Claims about empirical superiority are intentionally phrased as objectives or hypotheses until our own experiments are complete.

## Draft

Single-cell RNA sequencing (scRNA-seq) has become a central technology for profiling gene expression at single-cell resolution, enabling the discovery of cellular heterogeneity, rare cell populations, developmental programs, and disease-associated transcriptional states. Among the standard analytical steps in scRNA-seq studies, cell clustering is particularly important because it provides the basis for cell-type discovery, marker gene identification, atlas construction, and downstream biological interpretation. A reliable clustering method should therefore learn cell representations that are not only numerically separable, but also biologically meaningful and robust across tissues, protocols, and levels of sparsity.

Despite substantial progress, clustering scRNA-seq data remains challenging. A typical expression matrix contains thousands to tens of thousands of genes, while each cell has only a small fraction of observed nonzero expression values. The resulting data are high-dimensional, sparse, noisy, and often affected by batch effects. Classical clustering pipelines rely on low-dimensional projections and distance metrics that may miss nonlinear transcriptional structure. Graph-based methods can encode cell-cell relationships, but their performance is sensitive to graph construction and may suffer from over-smoothing or representation collapse. Recently, single-cell foundation models such as Geneformer, scGPT, and scFoundation have demonstrated the promise of large-scale pretraining for transferable cell representations \citep{theodoris2023geneformer,cui2024scgpt,hao2024scfoundation}. However, general-purpose representations are not necessarily optimized for clustering, and recent benchmark evidence suggests that task-specific clustering objectives remain important \citep{xu2026scclubench}.

Masked autoencoding provides a promising self-supervised alternative for clustering-oriented representation learning. In particular, scMAE perturbs gene expression values through gene-wise shuffling, trains an encoder to reconstruct the original expression matrix, and uses an auxiliary mask prediction objective to identify corrupted positions \citep{fang2024scmae}. This design encourages the encoder to capture dependencies among genes and has shown strong clustering performance on multiple real scRNA-seq datasets, including the ability to detect rare cell types. These results suggest that masked reconstruction can act as more than a denoising task: it can provide a useful pretext objective for learning cell embeddings that preserve cell-type structure.

However, existing masked autoencoder approaches for scRNA-seq still have two important limitations. First, the encoder is often implemented as a multilayer perceptron that compresses the expression vector of each cell independently. This design treats each cell as a flat feature vector and does not explicitly model the two-dimensional structure of the expression matrix. In scRNA-seq data, each entry is defined jointly by a cell and a gene. Therefore, the expression value of a gene in a cell should be interpreted in relation to both the other genes in the same cell and the distribution of the same gene across other cells. Second, random masking does not distinguish between informative and uninformative positions. In highly sparse matrices, many masked entries may correspond to zero or easily predictable values, producing weak training signals. Conversely, a smaller set of genes and cells may carry greater information for distinguishing cell states, rare populations, or pathway-level variation.

The success of recent tabular foundation models provides a useful architectural cue. TabPFN models tabular data by exploiting its two-dimensional structure, using attention across features within a sample and across samples within a feature \citep{hollmann2025tabpfn}. Although scRNA-seq clustering is an unsupervised biological representation learning problem rather than a supervised tabular prediction task, the structural analogy is useful: cells correspond to rows, genes correspond to columns, and each expression value is embedded in both row-wise and column-wise contexts. A masked reconstruction model for scRNA-seq should therefore move beyond estimating

$$
p(x_{ij}\mid X_{i,-j}),
$$

and should instead approximate a context-aware objective:

$$
p(x_{ij}\mid X_{i,-j}, X_{-i,j}),
$$

where \(X_{i,-j}\) denotes the expression context of other genes in the same cell, and \(X_{-i,j}\) denotes the population-level context of the same gene across other cells. Such a formulation can introduce population information without requiring an explicit cell-cell graph.

In this work, we propose CAAM-scMAE, a context-aware adversarial axial masked autoencoder for scRNA-seq clustering. CAAM-scMAE builds on the masked reconstruction framework of scMAE, but introduces two components designed for sparse expression matrices. First, a lightweight bi-axial context encoder models gene-module interactions within cells and incorporates context from representative cells through cell-axis attention. Second, a constrained adversarial mask selector learns to select more informative masked positions under fixed mask budgets, label-free constraints, coverage regularization, and sparsity-aware diagnostics. During training, the model reconstructs the original expression and predicts the mask from corrupted inputs; during inference, only the encoder output is used for downstream clustering.

The goal of CAAM-scMAE is to convert masked autoencoding from a purely random corruption task into a clustering-oriented context recovery task. We will evaluate this hypothesis under a scCluBench-style protocol across diverse datasets, comparing against classical, deep learning, graph-based, and masked autoencoder baselines. Beyond standard clustering metrics such as ACC, NMI, and ARI, our evaluation will include rare cell type discovery, marker-overlap annotation, embedding distinguishability, mask diagnostics, runtime, and memory usage. This design is intended to clarify not only whether CAAM-scMAE improves clustering performance, but also which components contribute to robust and biologically meaningful single-cell representations.

## Contributions

1. We introduce a bi-axial masked autoencoder for scRNA-seq clustering that explicitly models gene-axis and cell-axis contexts in sparse expression matrices.
2. We formulate constrained adversarial masking as a label-free and budget-controlled pretext task, aiming to select more informative reconstruction targets while avoiding shortcut corruption.
3. We establish an evaluation plan aligned with scCluBench-style benchmarking, covering clustering accuracy, rare cell discovery, embedding quality, biological annotation, mask behavior, and scalability.
4. We separate the effects of encoder architecture, mask strategy, and reconstruction objective through staged ablation experiments.

## Claims To Verify Experimentally

- CAAM-scMAE improves ACC, NMI, and ARI over scMAE-compatible reproduction.
- Bi-axial context attention improves embedding distinguishability without graph-induced collapse.
- Constrained adversarial masking improves informative reconstruction without destabilizing training.
- The selected full model remains computationally feasible under high sparsity and large gene dimensions.
