# Experiments Protocol v0

Updated: 2026-06-26

## 0. Purpose

This protocol defines the experiments required to support the CAAM-scMAE manuscript. The goal is not only to report higher clustering metrics, but to test the central mechanism:

> Does context-aware masked reconstruction with constrained informative masking learn better scRNA-seq clustering embeddings than random masked autoencoding?

The protocol is staged to prevent premature complexity. Each stage must pass basic sanity checks before moving to the next.

## 1. Datasets

### Development Dataset

Current local dataset:

```text
数据文件/SRP182008.h5ad
```

Use this dataset first for smoke tests, shape checks, runtime checks, and debugging.

### Benchmark Datasets

Formal benchmark should follow a scCluBench-style protocol:

```text
multiple tissues
multiple cell counts
multiple sparsity levels
multiple gene dimensions
multiple seeds
```

Minimum formal subset:

```text
small datasets: < 5,000 cells
medium datasets: 5,000-20,000 cells
large datasets: > 20,000 cells
high sparsity: > 80% zeros
many cell types: >= 10 classes
rare cell datasets: at least one class < 1%
```

## 2. Compared Methods

### Required Baselines

```text
PCA + KMeans
PCA + Leiden
SC3
DEC
DESC
scDeepCluster
scDCC
scNAME
scMAE
AttentionAE-sc
scDSC
scGNN
scCDCG
```

If runtime is limited, first compare:

```text
PCA + KMeans
PCA + Leiden
scMAE
scNAME
AttentionAE-sc
scCDCG
CAAM-scMAE variants
```

### CAAM-scMAE Variants

```text
V0: scMAE-compatible baseline
  random mask + gene-wise shuffle + MLP encoder

V1: gene-axis encoder
  random mask + gene-wise shuffle + gene-module attention

V2: bi-axial encoder
  random mask + gene-wise shuffle + gene-axis + cell-axis attention

V3: informative mask selector
  learned selector + MLP encoder

V4: full CAAM-scMAE
  learned selector + bi-axial encoder
```

Do not put all variants in the main benchmark table. The main table should include only:

```text
CAAM-scMAE full
```

Internal variants belong in ablation tables.

## 3. Metrics

### Clustering Metrics

```text
ACC
NMI
ARI
macro-F1 after best mapping
```

### Biological Interpretation Metrics

```text
marker-overlap annotation accuracy
rare cell recovery
top marker gene overlap
cluster-cell type Sankey consistency
```

### Embedding Quality Metrics

```text
silhouette score
cell-type ASW
embedding cosine similarity distribution
collapse score
UMAP visualization
```

### Mask Diagnostics

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

### Scalability Metrics

```text
runtime_seconds
peak_gpu_memory_mb
parameter_count
failed_by_oom
failed_by_nan
```

## 4. Stage 1: scMAE-Compatible Reproduction

### Goal

Verify that our training pipeline can reproduce a strong masked autoencoder baseline before adding CAAM components.

### Natural-Language Algorithm

For each mini-batch, randomly select a fixed fraction of expression entries. For each gene, shuffle expression values across cells. Replace masked entries with shuffled values, train the encoder to reconstruct the original expression matrix, and train the mask head to predict which entries were corrupted.

### LaTeX

$$
X'_{ij}=X_{\pi_j(i),j}
$$

$$
\widetilde{X}=(1-M)\odot X+M\odot X'
$$

$$
\mathcal{L}_{base}
=
\frac{1}{NG}\sum_{i,j}(1+\lambda M_{ij})(\widehat{X}_{ij}-X_{ij})^2
+\gamma\mathrm{BCEWithLogits}(\widehat{M},M)
$$

### Practical Code Interface

```python
result = train_variant(
    variant="scmae_compatible",
    corruption="gene_wise_shuffle",
    mask_policy="random",
    encoder="mlp",
    mask_ratio=0.3,
    reconstruction_loss="weighted_mse",
    mask_loss="bce",
)
```

### Pass Criteria

```text
training loss decreases
mask prediction loss is finite
embedding_final.npy exists
metrics.json exists
no NaN
no OOM
actual_mask_ratio within 5% of target
```

## 5. Stage 2: Encoder Ablation

### Goal

Test whether explicit two-dimensional context improves cell embeddings beyond an MLP encoder.

### Variants

```text
MLP encoder
gene-axis module encoder
bi-axial context encoder
```

### Key Question

```text
Does bi-axial context improve clustering metrics and embedding distinguishability without inducing representation collapse?
```

### Required Outputs

```text
attention_stats.json
context_selection.json
embedding_similarity.json
```

## 6. Stage 3: Mask Strategy Ablation

### Goal

Test whether learned mask selection produces more informative reconstruction tasks than random masking.

### Variants

```text
random mask
observed-only random mask
coverage-regularized learned mask
entropy-regularized learned mask
full constrained selector
```

### Key Question

```text
Does learned masking improve downstream clustering, or only make reconstruction harder?
```

### Failure Cases

```text
mask collapses to a small set of genes
mask selector always selects high-expression genes
training becomes unstable
student loss increases but clustering metrics do not improve
zero_to_zero_fraction remains too high
effective_changed_fraction remains too low
```

## 7. Stage 4: Formal Benchmark

### Goal

Compare full CAAM-scMAE against established baselines under a standardized multi-dataset protocol.

### Seeds

Use at least:

```text
3 seeds for development
5 seeds for formal benchmark if compute allows
```

### Main Table

```text
rows: methods
columns: datasets
values: ACC / NMI / ARI mean ± std
```

### Ablation Table

```text
V0: scMAE-compatible
V1: gene-axis
V2: bi-axial
V3: selector only
V4: full
```

### Diagnostic Table

```text
mask_entropy
effective_changed_fraction
zero_to_zero_fraction
runtime
memory
parameter_count
```

## 8. Required Artifacts Per Run

Every run must save:

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

For bi-axial variants:

```text
context_indices.npy
context_selection.json
attention_stats.json
gene_module_assignment.npz
```

For learned-mask variants:

```text
selector_stats.json
gene_mask_frequency.npy
mask_epoch_overlap.json
```

## 9. Manuscript Figures

### Figure 1: Method Overview

```text
X -> mask selector -> corrupted X -> bi-axial encoder -> embedding -> reconstruction + mask prediction
```

### Figure 2: Overall Benchmark

```text
ACC/NMI/ARI across datasets
average rank
statistical test
```

### Figure 3: Ablation

```text
V0 vs V1 vs V2 vs V3 vs V4
```

### Figure 4: Embedding and Rare Cell Analysis

```text
UMAP
rare cell recovery
marker-overlap annotation
```

### Figure 5: Mask Diagnostics

```text
mask entropy
gene coverage
zero-to-zero fraction
effective changed fraction
```

## 10. Rules Before Claiming Success

Do not claim CAAM-scMAE is better than scMAE unless:

```text
1. V0 scMAE-compatible baseline runs correctly.
2. Full model improves clustering metrics on multiple datasets.
3. Improvement is not explained only by parameter count.
4. Mask diagnostics show non-degenerate behavior.
5. Runtime and memory are acceptable.
6. Rare cell / marker-overlap analysis does not contradict clustering metrics.
```

If adversarial masking fails but bi-axial context works, revise method name:

```text
from: Context-aware Adversarial Axial Masked Autoencoder
to: Context-aware Axial Masked Autoencoder
```

If bi-axial context fails but mask selection works, revise method name:

```text
from: CAAM-scMAE
to: Adaptive Masked scMAE
```
