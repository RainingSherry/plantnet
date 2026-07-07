# rank58_fuzzy_rough_core_boundary_full

Independent-full scMAE candidate based on fuzzy rough iterative computation for scRNA-seq gene selection.

## Theory basis

The source paper defines fuzzy symmetric relations between cells and fuzzy rough lower/upper approximations to handle uncertainty and noise in high-dimensional scRNA-seq data. The improvement report recommends mapping this family of methods to core/boundary clustering constraints rather than directly forcing hard labels.

## scMAE integration

This candidate fills the boundary / clustering-head gap:

- scMAE mask prediction and masked expression reconstruction remain the primary loss.
- After a warmup, KMeans initializes fuzzy cluster centers in latent space.
- High-membership cells act as fuzzy lower-approximation core cells and receive a KL self-training loss.
- Low-confidence boundary cells receive no hard pull, which is the rough-set veto against over-confident mixing or clustering.

## Expression semantics

`scaled_expr` may be used only as encoder input. Reconstruction uses unscaled log-expression. No count distribution is used in this candidate.

## NeighborMix

NeighborMix is not used. The relationship is independent and complementary; `mixed_cell_fraction=0.0`.

## Outputs

Each run writes `embedding_final.npy`, `labels.npy`, `fuzzy_membership.npy`, `training_history.json`, `diagnostics.json`, `summary.json`, `args.json`, and optional fixed-k metrics. Smoke/screen rows are written only to `新模型独立快筛单次结果.csv` and summarized in `新模型独立快筛汇总结果.csv`.
