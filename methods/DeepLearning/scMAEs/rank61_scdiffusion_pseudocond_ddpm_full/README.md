# rank61_scdiffusion_pseudocond_ddpm_full

Independent-full scMAE candidate based on scDiffusion.

## Theory basis

scDiffusion combines an autoencoder latent space with a diffusion backbone and classifier guidance for controlled single-cell generation. The improvement report recommends using diffusion scheduling and denoising targets as a lightweight auxiliary for scMAE, without training a heavy generator or evaluating generated pseudo-cells as real cells.

## scMAE integration

This candidate fills the robust denoising / conditional target gap:

- scMAE mask prediction and masked expression reconstruction remain the primary loss.
- A linear DDPM beta schedule is used in the scMAE latent space.
- Pseudo conditions are obtained by unsupervised KMeans on log-expression SVD features, not by ground-truth labels.
- A time- and pseudo-condition-conditioned network predicts the injected latent noise.

## Expression semantics

`scaled_expr` may be used only as encoder input. Pseudo-condition discovery and reconstruction targets use unscaled log-expression. No generated cells enter KMeans evaluation.

## NeighborMix

NeighborMix is not used. The relationship is independent and complementary; this implementation never mixes cells and reports `mixed_cell_fraction=0.0`.

## Source note

The GitHub repository `https://github.com/EperLuo/scDiffusion` was cloned and inspected. The relevant source files are `guided_diffusion/gaussian_diffusion.py`, `guided_diffusion/script_util.py`, `guided_diffusion/train_util.py`, `guided_diffusion/cell_model.py`, and the README training workflow.

## Outputs

Each run writes `embedding_final.npy`, `labels.npy`, `pseudo_conditions.npy`, `training_history.json`, `diagnostics.json`, `summary.json`, `args.json`, and optional fixed-k metrics. Smoke/screen rows are written only to `新模型独立快筛单次结果.csv` and summarized in `新模型独立快筛汇总结果.csv`.
