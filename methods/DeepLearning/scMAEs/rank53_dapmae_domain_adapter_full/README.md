# rank53_dapmae_domain_adapter_full

Independent-full scMAE candidate based on DAP-MAE, "Domain-Adaptive Point Cloud Masked Autoencoder for Effective Cross-Domain Learning".

## Theory basis

The PDF and GitHub code emphasize two transferable components: a heterogeneous domain adapter (HDA) with domain-specific branches during pretraining, and a domain feature generator (DFG) guided by contrastive/domain objectives. I inspected the public repository at `https://github.com/CVI-SZU/DAP-MAE`, especially `models/HDAdapter.py`, `models/DAP_MAE.py`, and `ContrastiveLoss.py`.

## scMAE integration

This candidate fills scMAE's target/domain-adapter gap:

- It builds unsupervised pseudo-domains from log-expression SVD plus KMeans.
- `DAPMAEScMAE` uses DAP-style heterogeneous domain adapter branches and domain tokens.
- The original scMAE mask prediction and masked expression reconstruction are retained.
- A DFG head predicts pseudo-domain prototypes, with domain CE and within-batch contrastive guidance.

## Expression semantics

`scaled_expr` may be used only as encoder input. Reconstruction targets and pseudo-domain SVD features are derived from unscaled log-expression. No scaled expression is treated as count data.

## NeighborMix

NeighborMix is not used. The relationship is independent and complementary; `mixed_cell_fraction=0.0`.

## Outputs

Each run writes `embedding_final.npy`, `labels.npy`, `training_history.json`, `diagnostics.json`, `summary.json`, `args.json`, and optional fixed-k metrics. Smoke/screen rows are written only to `新模型独立快筛单次结果.csv` and summarized in `新模型独立快筛汇总结果.csv`.
