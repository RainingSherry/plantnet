# Experiments English Draft v0

Updated: 2026-06-26

> This section is written as an experimental design draft. Result-specific claims must be inserted only after the corresponding experiments are completed.

## 1. Experimental Overview

The experiments are designed to test whether CAAM-scMAE learns clustering-oriented single-cell representations through context-aware masked reconstruction. The evaluation has four goals. First, we verify that our implementation can reproduce a strong scMAE-compatible masked autoencoder baseline. Second, we isolate the effect of the proposed bi-axial context encoder. Third, we evaluate whether constrained mask selection provides more informative reconstruction targets than random masking. Fourth, we compare the selected full model against established single-cell clustering baselines under a multi-dataset benchmark protocol.

To avoid conflating engineering complexity with scientific contribution, the evaluation is staged. The full model is not evaluated before the scMAE-compatible baseline and the encoder ablations pass basic sanity checks. Each run saves embeddings, clustering labels, training curves, mask diagnostics, runtime, memory usage, and model size, allowing us to diagnose not only final clustering metrics but also failure modes such as degenerate masking, representation collapse, or excessive computational cost.

Correction-pipeline update after Phase 14:

```text
The constrained AdvMask selector should no longer be treated as a default main component.
In the corrected 3-epoch development triage, AdvMask had real nonzero generator gradients,
but its mean ARI gain over random masking was below the predefined seed-variation threshold.
Therefore subsequent experiments should not run Axial/full as a rescue path for AdvMask.
```

Development smoke update after attention/context check:

```text
An attention/context-only smoke was run after AdvMask was downgraded.
It compared the existing Axial encoder against the standard MLP and a parameter-matched MLP
under random masking, scmae_shuffle corruption, seed 42, and 3 epochs on three development datasets.
The existing Axial encoder underperformed both MLP controls on all three datasets.
This does not invalidate all attention ideas, but it blocks the current Axial implementation
from being used as the main method route without a new mechanism-level redesign.
```

## 2. Datasets

We first use the local development dataset:

```text
数据文件/SRP182008.h5ad
```

for smoke testing, shape validation, and debugging. This file corresponds to the Arabidopsis root dataset SRP182008 \citep{zhang2019arabidopsisroot}. The local matrix contains 13,514 cells and 53,678 genes, with 17,378,932 nonzero entries and an overall sparsity of 97.60%. The file provides `Celltype` labels with 15 categories and `Seurat_clusters` with 24 clusters; these annotations are used only for evaluation. The dataset is therefore well suited for testing sparse masked reconstruction and plant single-cell development, but it is not sufficient by itself for batch-effect or multi-condition claims.

Formal evaluation should follow a scCluBench-style setting \citep{xu2026scclubench}, selecting datasets that cover different tissues, species, cell counts, gene dimensions, sparsity levels, and numbers of cell types. The formal subset should include small datasets with fewer than 5,000 cells, medium datasets with 5,000 to 20,000 cells, large datasets with more than 20,000 cells, high-sparsity datasets with more than 80% zero entries, and at least one dataset containing rare cell populations.

All datasets will be processed using a fixed preprocessing protocol. Unless otherwise specified, labels are excluded from training and used only for post hoc evaluation. This separation is necessary because CAAM-scMAE is intended as an unsupervised clustering representation learner, not a supervised classifier.

## 3. Baselines

We compare CAAM-scMAE against representative methods from several methodological families. Classical and community detection baselines include PCA followed by K-means, PCA followed by Leiden clustering, SC3 \citep{kiselev2017sc3}, Louvain \citep{blondel2008louvain}, Leiden \citep{traag2019leiden}, and Seurat-style graph clustering \citep{stuart2019seurat}. Deep clustering baselines include DEC \citep{xie2016dec}, DESC \citep{li2020desc}, scDeepCluster \citep{tian2019scdeepcluster}, scDCC \citep{tian2021scdcc}, and scziDesk \citep{chen2020sczidesk}. Graph-based and attention-based baselines include scGNN \citep{wang2021scgnn}, scDSC \citep{gan2022scdsc}, AttentionAE-sc \citep{li2023attentionae}, and scCDCG \citep{xu2024sccdcg}. Self-supervised masked or contrastive baselines include scNAME \citep{wan2022scname} and scMAE \citep{fang2024scmae}.

If computational resources are limited, we will first use a reduced baseline set consisting of PCA + K-means, PCA + Leiden, scNAME, scMAE, AttentionAE-sc, scCDCG, and CAAM-scMAE variants. The full benchmark table should include only the selected full CAAM-scMAE model, while internal variants are reported in ablation studies.

## 4. CAAM-scMAE Variants

The original method-paper plan evaluated five internal variants:

```text
V0: scMAE-compatible baseline
    random mask + gene-wise shuffle + MLP encoder

V1: gene-axis encoder
    random mask + gene-wise shuffle + gene-module attention

V2: bi-axial encoder
    random mask + gene-wise shuffle + gene-axis attention + cell-axis context attention

V3: selector-only model
    constrained mask selector + MLP encoder

V4: full CAAM-scMAE
    constrained mask selector + bi-axial context encoder
```

This design separates the contribution of the encoder architecture from the contribution of the mask policy. If V2 improves over V0 but V3 does not, the main contribution should be framed around context-aware axial encoding. If V3 improves over V0 but V2 does not, the paper should be reframed around adaptive or informative masking. If V4 improves over both V2 and V3, the full CAAM-scMAE hypothesis is supported.

Current evidence revises this plan:

```text
V3 did not pass Phase 14 under scmae_shuffle corruption.
V4 should not be run or interpreted as Axial + AdvMask synergy under the current BDD.
AdvMask may be retained only as a supplementary negative-diagnostic variant.
The existing Axial encoder also failed a 3-dataset, seed-42 development smoke
against both a standard MLP and a parameter-matched MLP.
```

If the attention/context hypothesis is revisited, it should be separated from AdvMask and tested against a parameter-matched MLP control under a new, explicit research gate.

## 5. Metrics

### Clustering Accuracy

We report Accuracy (ACC), Normalized Mutual Information (NMI), Adjusted Rand Index (ARI), and macro-F1 after best label mapping. Since ground-truth labels are not used during training, these metrics are used only for evaluation.

### Biological Interpretation

We evaluate whether clustering results correspond to biologically meaningful cell populations using marker-overlap annotation, rare cell recovery, top marker gene overlap, and Sankey diagrams comparing predicted clusters with reference labels. These analyses are included because clustering metrics alone may hide biologically important failure modes.

### Embedding Quality

We measure silhouette score, cell-type ASW, cosine similarity distributions, and representation collapse scores. We also visualize embeddings with UMAP. This is especially important for comparing CAAM-scMAE with graph neural network methods, where over-smoothing can make embeddings overly similar across cell types.

### Mask Diagnostics

For all masked autoencoding variants, we report:

```text
actual_mask_ratio_global
actual_mask_ratio_observed
zero_to_zero_fraction
effective_changed_fraction
gene_mask_frequency
mask_entropy
top_masked_gene_concentration
mask_overlap_across_epochs
```

These diagnostics test whether the masking task is informative. A model may optimize reconstruction loss while masking mostly uninformative zero entries; such behavior should not be interpreted as successful representation learning.

### Scalability

We report runtime, peak GPU memory, parameter count, out-of-memory failures, and NaN-loss failures. A method intended for scRNA-seq clustering must be robust not only in accuracy but also in computational feasibility.

## 6. Ablation Studies

The first ablation isolates the encoder. We compare the MLP encoder, the gene-axis module encoder, and the bi-axial context encoder while keeping random masking and gene-wise shuffle corruption fixed. This tests whether explicit two-dimensional context improves clustering representations without introducing graph-induced collapse.

The second ablation isolates the mask policy. We compare random masking, observed-only random masking, coverage-regularized learned masking, entropy-regularized learned masking, and the full constrained selector. This tests whether learned masks improve downstream clustering rather than merely increasing reconstruction difficulty.

The third ablation controls for model size. If the full model has more parameters than the baseline, we include a parameter-matched MLP variant. This is necessary to show that improvements are due to the proposed inductive biases rather than capacity alone.

## 7. Success and Failure Criteria

We will not claim that CAAM-scMAE improves over scMAE unless the following conditions are met:

```text
1. The V0 scMAE-compatible baseline trains correctly.
2. The full model improves clustering metrics on multiple datasets.
3. The improvement is not explained solely by parameter count.
4. Mask diagnostics show non-degenerate behavior.
5. Runtime and memory are acceptable.
6. Rare cell and marker-overlap analyses do not contradict the clustering metrics.
```

If constrained adversarial masking is unstable or fails to improve clustering, the method should be renamed and reframed as a context-aware axial masked autoencoder. If bi-axial context fails but mask selection works, the method should be reframed as an adaptive masked scMAE. This protects the manuscript from depending on a component that is not empirically supported.

Phase 14 already triggers the first condition:

```text
AdvMask trains but fails to improve clustering beyond seed-level variation.
The current manuscript should downgrade or remove adversarial masking from its main contribution list.
```

The attention/context smoke triggers an additional guardrail:

```text
The current Axial encoder should also be downgraded.
Future attention work must be treated as a new design question, not as evidence that the existing CAAM axial path works.
```

## 8. Required Artifacts

Each training run must save:

```text
config.json
metrics.json
embedding_final.npy
cluster_labels.npy
training_history.json
runtime.json
param_count.json
mask_diagnostics.json
```

Bi-axial variants must additionally save:

```text
context_indices.npy
context_selection.json
attention_stats.json
gene_module_assignment.npz
```

Learned-mask variants must additionally save:

```text
selector_stats.json
gene_mask_frequency.npy
mask_epoch_overlap.json
```

These artifacts allow independent inspection of the training process, representation quality, and mask behavior.
