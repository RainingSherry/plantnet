# Stage 1 Training Smoke Result: scMAE-Compatible MLP

Date: 2026-06-26

## 1. Purpose

This run verifies that the Stage 1 scMAE-compatible training entry point works on the cached SRP182008 HVG subset. It is a smoke test, not a performance result for the manuscript.

## 2. Cache Generation

Cache command:

```powershell
python "实验代码\stage1_scmae_compatible.py" `
  --save-dir "实验结果\stage1_srp182008\cache_hvg2000_cells2048_seed42" `
  --n-top-genes 2000 `
  --max-cells 2048 `
  --mask-ratio 0.3 `
  --seed 42 `
  --diagnostics-only `
  --save-cache `
  --force
```

Cache artifacts:

```text
x_hvg_log1p.npy
labels_Celltype.npy
labels_Seurat_clusters.npy
hvg_gene_indices.npy
hvg_gene_names.npy
selected_cell_indices.npy
mask_diagnostics.json
data_summary.json
```

## 3. Training Command

Training command:

```powershell
& "S:\Up\Code\Deeplearning\Miniconda3\envs\machine\python.exe" `
  "实验代码\stage1_train_mlp_from_cache.py" `
  --cache-dir "实验结果\stage1_srp182008\cache_hvg2000_cells2048_seed42" `
  --save-dir "实验结果\stage1_srp182008\train_mlp_hvg2000_cells2048_seed42_epochs2" `
  --epochs 2 `
  --batch-size 128 `
  --mask-ratio 0.3 `
  --latent-dim 64 `
  --hidden-dim 512 `
  --seed 42 `
  --device auto `
  --force
```

## 4. Natural-Language Interpretation

The two-epoch MLP masked-autoencoder smoke test ran successfully on CPU. The training loss decreased from 1.87 to 1.39, showing that the training loop, corruption function, reconstruction head, mask head, and artifact saving path are functional.

The clustering metrics are low, which is expected because this run used only 2 epochs and 2,048 cells. These numbers should not be used as manuscript performance claims.

## 5. Training History

| Epoch | Loss | Reconstruction Loss | Mask Loss |
|---:|---:|---:|---:|
| 1 | 1.8677 | 1.1840 | 0.6837 |
| 2 | 1.3912 | 0.7495 | 0.6417 |

## 6. Smoke Metrics

| Evaluation target | Classes | ARI | NMI | Macro-F1 unmapped |
|---|---:|---:|---:|---:|
| Celltype | 15 | 0.0446 | 0.1586 | 0.0533 |
| Seurat_clusters | 24 | 0.0838 | 0.2432 | 0.0409 |

Embedding:

```text
embedding_final.npy shape: 2048 × 64
```

Model size:

```text
parameters: 1,318,368
trainable parameters: 1,318,368
```

Runtime:

```text
elapsed_seconds: 7.71
device: CPU
cuda_available in machine env: false
```

## 7. Mask Diagnostics During Training

| Quantity | Value |
|---|---:|
| actual global mask ratio | 0.30000 |
| actual observed mask ratio | 0.29974 |
| masked observed fraction among masked | 0.32559 |
| zero-to-zero fraction among masked | 0.48444 |
| effective changed fraction among masked | 0.47550 |
| masked zero fraction among masked | 0.67441 |
| nonzero-to-zero fraction among masked | 0.18994 |
| zero-to-nonzero fraction among masked | 0.18997 |
| normalized mask entropy | 0.99996 |

These diagnostics are consistent with the diagnostics-only run: random gene-wise shuffle masking in this sparse HVG matrix still produces nearly half zero-to-zero masked entries.

## 8. Artifacts

Saved under:

```text
实验结果/stage1_srp182008/train_mlp_hvg2000_cells2048_seed42_epochs2
```

Files:

```text
config.json
embedding_final.npy
mask_diagnostics.json
metrics.json
param_count.json
runtime.json
training_history.json
```

## 9. Manuscript Use

This result can be used to say:

```text
We established a working scMAE-compatible Stage 1 training pipeline and verified that the loss decreases in a short smoke test.
```

Do not use it to say:

```text
CAAM-scMAE improves clustering.
scMAE-compatible baseline is strong.
The model reaches publishable performance.
```

## 10. Next Step

Run a proper development baseline:

```text
HVG 2000
all 13,514 cells
mask_ratio 0.2 / 0.3 / 0.4
seeds 42 / 43 / 44
epochs 50 or early stopping
```

Before that, improve evaluation:

```text
1. add Hungarian-mapped ACC and macro-F1;
2. save cluster labels;
3. optionally run PCA + KMeans baseline on the same cached subset;
4. test machine or pytorch environment with CUDA if available.
```
