# Stage 1 Diagnostics Result: SRP182008 HVG-2000 Random Mask

Date: 2026-06-26

## 1. Run Summary

Command:

```powershell
python "实验代码\stage1_scmae_compatible.py" `
  --save-dir "实验结果\stage1_srp182008\diagnostics_hvg2000_cells2048_seed42" `
  --n-top-genes 2000 `
  --max-cells 2048 `
  --mask-ratio 0.3 `
  --seed 42 `
  --diagnostics-only `
  --force
```

Artifacts:

```text
实验结果/stage1_srp182008/diagnostics_hvg2000_cells2048_seed42/config.json
实验结果/stage1_srp182008/diagnostics_hvg2000_cells2048_seed42/data_summary.json
实验结果/stage1_srp182008/diagnostics_hvg2000_cells2048_seed42/mask_diagnostics.json
实验结果/stage1_srp182008/diagnostics_hvg2000_cells2048_seed42/runtime.json
实验结果/stage1_srp182008/diagnostics_hvg2000_cells2048_seed42/hvg_gene_indices.npy
实验结果/stage1_srp182008/diagnostics_hvg2000_cells2048_seed42/selected_cell_indices.npy
实验结果/stage1_srp182008/diagnostics_hvg2000_cells2048_seed42/gene_mask_frequency.npy
```

This was a diagnostics-only run because the current Python environment does not have `torch` installed.

## 2. Natural-Language Interpretation

The run tested random scMAE-style masking on a 2,048-cell subset and 2,000 highly variable genes from SRP182008. The target mask ratio was 30%, and the observed global mask ratio was almost exactly 30%.

However, because the dataset is highly sparse, only about one third of masked entries were originally nonzero. Nearly half of all masked positions were zero-to-zero after gene-wise shuffling. In other words, a large fraction of nominally masked entries did not create a meaningful reconstruction perturbation.

This supports a key motivation for CAAM-scMAE:

```text
In highly sparse scRNA-seq data, random mask ratio is not equivalent to effective corruption ratio.
```

Therefore, later variants should report not only the target mask ratio, but also:

```text
zero_to_zero_fraction
effective_changed_fraction
actual_mask_ratio_observed
```

## 3. Key Numbers

| Quantity | Value |
|---|---:|
| selected cells | 2,048 |
| selected genes | 2,000 |
| total selected entries | 4,096,000 |
| target mask ratio | 0.3000 |
| actual global mask ratio | 0.29998 |
| observed fraction in selected HVG matrix | 0.32587 |
| masked observed fraction among masked entries | 0.32593 |
| zero-to-zero fraction among masked entries | 0.48322 |
| effective changed fraction among masked entries | 0.47869 |
| masked zero fraction among masked entries | 0.67407 |
| nonzero-to-zero fraction among masked entries | 0.19141 |
| zero-to-nonzero fraction among masked entries | 0.19084 |
| normalized mask entropy | 0.99993 |

## 4. LaTeX Interpretation

The nominal mask ratio is:

$$
\rho
=
\frac{\sum_{i,j}M_{ij}}{NG}
\approx
0.300.
$$

The observed-entry mask ratio is:

$$
\rho_{\mathrm{obs}}
=
\frac{\sum_{i,j}M_{ij}\mathbf{1}(X_{ij}>0)}
{\sum_{i,j}\mathbf{1}(X_{ij}>0)}
\approx
0.300.
$$

But the effective changed fraction is:

$$
\rho_{\mathrm{changed}}
=
\frac{\sum_{i,j}M_{ij}\mathbf{1}(X'_{ij}\neq X_{ij})}
{\sum_{i,j}M_{ij}}
\approx
0.479.
$$

And the zero-to-zero fraction is:

$$
\rho_{0\rightarrow 0}
=
\frac{\sum_{i,j}M_{ij}\mathbf{1}(X_{ij}=0)\mathbf{1}(X'_{ij}=0)}
{\sum_{i,j}M_{ij}}
\approx
0.483.
$$

Thus, under this random mask and gene-wise shuffle setting:

$$
\rho_{\mathrm{changed}} < \rho,
$$

meaning that the effective corruption signal is substantially weaker than the nominal mask budget.

## 5. Practical Code Entry

Script:

```text
实验代码/stage1_scmae_compatible.py
```

The script currently supports:

```text
h5ad CSR reading
log1p HVG selection by sparse variance
cell subsampling
random mask generation
gene-wise shuffle corruption diagnostics
artifact writing
```

It does not yet train a neural model because `torch` is unavailable in the current environment.

## 6. Manuscript Implication

This result can support a cautious motivation sentence:

```text
On the SRP182008 Arabidopsis root dataset, a random 30% mask over 2,000 highly variable genes produced an effective changed fraction of only 47.9%, with 48.3% of masked entries remaining zero after gene-wise shuffling. This indicates that nominal mask budgets can substantially overestimate effective corruption in sparse scRNA-seq matrices.
```

Do not use this as evidence that CAAM-scMAE improves clustering. It only supports the need for mask diagnostics and potentially informative or sparsity-aware mask policies.

## 7. Next Step

Implement a Torch-based training runner or install Torch, then run:

```text
HVG 2000
mask_ratio 0.3
epochs 2
seed 42
batch_size 128
```

Required training artifacts:

```text
embedding_final.npy
metrics.json
training_history.json
mask_diagnostics.json
runtime.json
param_count.json
```
