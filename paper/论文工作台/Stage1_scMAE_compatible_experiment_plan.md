# Stage 1 Experiment Plan: scMAE-Compatible Baseline

Updated: 2026-06-26

## 1. Purpose

Stage 1 tests whether our pipeline can train a strong scMAE-compatible masked autoencoder on the local plant scRNA-seq dataset before adding bi-axial attention or constrained adversarial masking.

This is the necessary first gate for the CAAM-scMAE paper. If this baseline fails, later improvements cannot be interpreted cleanly.

## 2. Dataset

Use:

```text
数据文件/SRP182008.h5ad
```

Dataset audit:

```text
13,514 cells × 53,678 genes
sparsity: 97.60%
Celltype labels: 15 categories, evaluation only
Seurat_clusters: 24 clusters, evaluation only
source: Zhang et al. 2019, Molecular Plant
```

BibTeX:

```text
zhang2019arabidopsisroot
he2024scplantdb
```

## 3. Algorithm: Natural-Language Description

For each mini-batch, we randomly choose a fixed fraction of gene expression entries to corrupt. For each gene, we shuffle expression values across cells in the same mini-batch, preserving the marginal distribution of that gene. Masked entries are replaced with their shuffled values, while unmasked entries remain unchanged.

The model receives the corrupted expression vector and produces:

```text
1. a cell embedding z_i;
2. a reconstructed expression vector x_hat_i;
3. a predicted mask vector m_hat_i.
```

Training optimizes a weighted reconstruction loss and a mask prediction loss. The final embedding is evaluated with clustering metrics, but cell labels are never used during training.

## 4. Algorithm: LaTeX Formulation

Given an expression matrix:

$$
X\in\mathbb{R}^{N\times G},
$$

sample a binary mask:

$$
M_{ij}\sim \mathrm{Bernoulli}(\rho).
$$

For each gene \(j\), sample a cell permutation:

$$
\pi_j:\{1,\dots,N\}\rightarrow\{1,\dots,N\}.
$$

The shuffled replacement is:

$$
X'_{ij}=X_{\pi_j(i),j}.
$$

The corrupted input is:

$$
\widetilde{X}=(1-M)\odot X+M\odot X'.
$$

The encoder and heads are:

$$
Z=f_\theta(\widetilde{X}),
$$

$$
\widehat{X}=g_\theta(Z),
$$

$$
\widehat{M}=h_\theta(Z).
$$

The weighted reconstruction loss is:

$$
\mathcal{L}_{rec}
=
\frac{1}{NG}
\sum_{i=1}^{N}
\sum_{j=1}^{G}
(1+\lambda M_{ij})
(\widehat{X}_{ij}-X_{ij})^2.
$$

The mask prediction loss is:

$$
\mathcal{L}_{mask}
=
\mathrm{BCEWithLogits}(\widehat{M},M).
$$

The total objective is:

$$
\mathcal{L}
=
\mathcal{L}_{rec}
+\gamma\mathcal{L}_{mask}.
$$

## 5. Algorithm: Practical Code Interface

Minimal PyTorch-style interface:

```python
result = train_variant(
    data_path="数据文件/SRP182008.h5ad",
    variant="scmae_compatible",
    input_mode="log1p_hvg",
    n_top_genes=2000,
    encoder="mlp",
    corruption="gene_wise_shuffle",
    mask_policy="random",
    mask_ratio=0.3,
    reconstruction_loss="weighted_mse",
    mask_loss="bce",
    labels_for_eval=["Celltype", "Seurat_clusters"],
    labels_for_train=None,
    save_dir="results/stage1_srp182008/scmae_compatible_seed42",
    seed=42,
)
```

Core corruption function:

```python
def gene_wise_shuffle_corruption(x, mask):
    batch_size, n_genes = x.shape
    x_prime = torch.empty_like(x)
    for j in range(n_genes):
        perm = torch.randperm(batch_size, device=x.device)
        x_prime[:, j] = x[perm, j]
    return x * (1.0 - mask) + x_prime * mask
```

Core loss:

```python
def scmae_compatible_loss(x, x_hat, mask, mask_logits, lambda_masked=4.0, gamma=1.0):
    rec_weight = 1.0 + lambda_masked * mask
    loss_rec = (rec_weight * (x_hat - x).pow(2)).mean()
    loss_mask = torch.nn.functional.binary_cross_entropy_with_logits(mask_logits, mask)
    return loss_rec + gamma * loss_mask
```

## 6. Experimental Grid

### Minimal Smoke Test

```text
genes: HVG 2000
mask_ratio: 0.3
seed: 42
epochs: 2
batch_size: 128
```

Pass criteria:

```text
no crash
no NaN
embedding_final.npy saved
metrics.json saved
mask_diagnostics.json saved
```

### Development Run

```text
genes: HVG 2000
mask_ratio: 0.2, 0.3, 0.4
seeds: 42, 43, 44
epochs: 50 or early stopping
batch_size: 128
```

### Full-Gene Feasibility Check

Run only after HVG version is stable:

```text
genes: all 53,678
mask_ratio: 0.3
seed: 42
epochs: 2
batch_size: memory-dependent
```

The full-gene check is for scalability diagnostics, not the first main result.

## 7. Metrics

Evaluate embeddings with:

```text
KMeans known-K using Celltype K=15
KMeans known-K using Seurat_clusters K=24
fixed-resolution Leiden
ACC
NMI
ARI
macro-F1 after Hungarian matching
silhouette / ASW
```

Note:

```text
Known-K KMeans is an oracle-style diagnostic.
Fixed-resolution Leiden is closer to unknown-K use.
Both should be reported separately.
```

## 8. Mask Diagnostics

Required diagnostics:

```text
actual_mask_ratio_global
actual_mask_ratio_observed
zero_to_zero_fraction
effective_changed_fraction
gene_mask_frequency
mask_entropy
top_masked_gene_concentration
```

Why this matters:

```text
SRP182008 has 97.60% sparsity. Random masks may land mostly on zeros, so a nominal 30% mask ratio may provide much weaker effective corruption.
```

## 9. Required Artifacts

Each run must save:

```text
config.json
metrics.json
embedding_final.npy
cluster_labels_celltype_kmeans.npy
cluster_labels_seurat_kmeans.npy
cluster_labels_leiden_fixed.npy
training_history.json
runtime.json
param_count.json
mask_diagnostics.json
```

## 10. Interpretation Rules

If Stage 1 performs poorly:

```text
Do not proceed directly to full CAAM-scMAE.
First check preprocessing, HVG selection, loss scaling, mask ratio, and zero-to-zero fraction.
```

If Stage 1 performs reasonably:

```text
Proceed to Stage 2 encoder ablation.
```

If Stage 1 performs well but mask diagnostics show mostly zero-to-zero corruption:

```text
This supports the paper's motivation for informative or sparsity-aware masking.
```

If Stage 1 performs well and mask diagnostics are already strong:

```text
The adversarial selector must show additional benefit beyond a strong random mask baseline.
```

## 11. Current Execution Status

Completed:

```text
diagnostics-only smoke
HVG 2000
max cells 2048
mask_ratio 0.3
seed 42

two-epoch MLP training smoke
HVG 2000
max cells 2048
mask_ratio 0.3
latent_dim 64
hidden_dim 512
seed 42
```

Result files:

```text
Stage1_diagnostics_result_2026-06-26.md
Stage1_training_smoke_result_2026-06-26.md
```

Current status:

```text
The scMAE-compatible Stage 1 training entry point runs in the machine Torch environment on CPU.
The smoke test is only a functionality check and should not be treated as manuscript performance evidence.
```

Next implementation step:

```text
Add Hungarian-mapped ACC and macro-F1, save KMeans cluster labels, then run the development baseline on all 13,514 cells.
```
