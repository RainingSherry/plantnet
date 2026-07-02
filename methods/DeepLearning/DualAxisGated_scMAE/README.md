# DualAxisGated-scMAE

A dual-axis attention encoder with a reliability-gated cell axis, built on the
rank13 DEC scMAE backbone. This is an **experimental development method**, not a
verified formal-benchmark method, and is **not** registered in any method
manifest.

## Core idea

Replace the plain MLP encoder of the scMAE backbone with a **dual-axis
encoder**, so cells can inform each other during encoding while rare/boundary
cells stay protected.

- **Gene axis (always on).** Expression `x[B,G]` is pooled into `M` gene modules
  through a hard, label-free `[G, M]` assignment (KMeans on gene profiles —
  "shared gene-row weights"), tokenized, and passed through a few
  self-attention layers. This is the per-cell representation base and is
  **never gated**: it plays the role of the self-anchored scMAE backbone.
- **Cell axis (reliability-gated).** The gene-axis cell summary cross-attends to
  `K` learnable prototypes (`K << N`, so it scales — no `N×N` attention). The
  resulting context is multiplied by the per-cell reliability `r_i ∈ [0,1]`, so
  core cells (`r≈1`) absorb cross-cell information while rare/boundary cells
  (`r≈0`) fall back to the pure gene-axis path (= pure scMAE behavior), exactly
  where neighbor smoothing is unsafe.

`r_i = neighbor_agreement · local_density` (see `reliability.py`, reused
verbatim from `GatedNeighborMix_scMAE`). It is recomputed from the current
embedding every few epochs after warmup.

The encoder output is a `latent [B, hidden_size]` vector — a drop-in
replacement for the MLP encoder — so every downstream head (mask predictor,
decoder, DEC cluster centers, the reliability-gated NeighborMix pseudo term) is
**unchanged** from the proven rank13 / GatedNeighborMix pipeline.

## Why plant data

Across the whole scMAE improvement search, every added mechanism (APA
representation alignment, CAAM adversarial mask, GatedNeighborMix,
AdaptiveGranularity anchor/adapt) failed to beat the backbone on **animal**
datasets, because scMAE is already saturated there (Quake ARI 0.957,
Limb_Muscle 0.990). This method instead targets **plant** single-cell data,
where the problem is open: classic Leiden (ACC ~0.61–0.69) often ties or beats
all deep methods, most deep methods were never run, and APA-scMAE ranks last
(0.16–0.43).

**Target to beat (ACC), from the plant benchmark:**

| Dataset   | Leiden | scMAE | best deep baseline |
| --------- | -----: | ----: | -----------------: |
| SRP171040 |  0.690 | 0.679 |       scDCC 0.527  |
| SRP182008 |  0.612 | 0.577 |  scGPT/scDCC ~0.68 |
| SRP235541 |  0.657 | 0.557 |       scDCC 0.567  |

First milestone: on ≥2 plant datasets, no collapse **and** ACC/NMI ≥ scMAE.

## Files

- `model.py` — `DualAxisEncoder` + `DualAxisGatedScMAE` (backbone with the
  encoder swapped; all other heads identical to the rank13 DEC scMAE)
- `reliability.py` — per-cell reliability field `r_i` (verbatim reuse)
- `loss.py` — gated scMAE + confidence-gated DEC KL + gated NeighborMix pseudo
  (verbatim reuse)
- `run.py` — data loading, PCA-KNN graph, gene-module construction,
  train/eval loop, smoke/screen stages

## Key implementation notes

- **Scaled encoder input, unscaled reconstruction target.** Inherited from
  rank13 and kept intact: reconstructing scaled→scaled drops the scMAE loss ~6×
  and lets DEC dominate → collapse. Do not change.
- **Gate is a training-time regularizer by default.** `extract_all` uses the
  full encoder (`r=None`) at evaluation, so the gate shapes parameters during
  training but the final embedding uses full model capacity. Pass
  `--eval_gate true` to apply the gate at eval as well (ablation).
- **Memory.** Gene axis runs over `M≈64` module tokens (not `G` genes); cell
  axis over `K≈16` prototypes (not `N` cells). No `N×N` or `G×G` attention.

## Usage

Environment: `/data/luolie/conda/envs/scssl_bench_py310/bin/python`
(the active `plantnet` conda env has no numpy).

Smoke (CPU, 3 epochs):

```bash
/data/luolie/conda/envs/scssl_bench_py310/bin/python \
  methods/DeepLearning/DualAxisGated_scMAE/run.py \
  --data_path /data/luolie/biopipeline/scCluBench/data/SRP224648.h5ad \
  --save_dir methods/DeepLearning/DualAxisGated_scMAE/runs/smoke/SRP224648/seed42 \
  --dataset_name SRP224648 --label_key auto --n_clusters 4 --smoke --no_cuda
```

Formal (GPU 1–6):

```bash
/data/luolie/conda/envs/scssl_bench_py310/bin/python \
  methods/DeepLearning/DualAxisGated_scMAE/run.py \
  --data_path /data/luolie/biopipeline/scCluBench/data/SRP171040.h5ad \
  --save_dir methods/DeepLearning/DualAxisGated_scMAE/runs/plant/SRP171040/seed42 \
  --dataset_name SRP171040 --n_clusters 12 --gpu 1
```

## Dual-axis / gate hyperparameters

| Flag                | Default | Meaning |
| ------------------- | ------: | ------- |
| `--n_gene_modules`  |      64 | gene modules `M` (gene-axis token count) |
| `--token_dim`       |      48 | per-module token dimension |
| `--n_prototypes`    |      16 | prototypes `K` for the cell axis |
| `--gene_layers`     |       2 | gene-axis self-attention layers |
| `--attn_heads`      |       4 | attention heads |
| `--attn_dropout`    |     0.1 | attention dropout |
| `--eval_gate`       |   false | apply reliability gate at eval too |
| `--pseudo_weight`   |     0.3 | NeighborMix pseudo term weight (0 = dual-axis alone) |
| `--alpha_min`       |     0.6 | min self-weight for fully-reliable cells |

## Built-in ablations

- **Dual-axis alone:** `--pseudo_weight 0` (turns off the NeighborMix branch)
- **Gate off:** set reliability to all-ones (no cell-axis retreat) to test
  whether the gate prevents collapse on rare cells
- **Module granularity:** sweep `--n_gene_modules 32 / 64 / 128`

## Status

- Smoke on SRP224648 (n_clusters=4, CPU, 3 epochs): passes end-to-end, no
  collapse, reliability spans 0–1 with core fraction ~0.86.
- Not yet run at full length / on GPU across the plant suite.
