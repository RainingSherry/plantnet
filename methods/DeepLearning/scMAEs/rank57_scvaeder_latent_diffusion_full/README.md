# rank57_scvaeder_latent_diffusion_full

Independent-full scMAE candidate based on scVAEDer.

## Theory basis

scVAEDer trains a VAE on scRNA-seq expression, then trains a diffusion model over the VAE latent embeddings so generated latent samples better match the real latent distribution than a simple Gaussian prior. The GitHub URL in the index required credentials in this environment, so this implementation is reconstructed from the PDF and the improvement report.

## scMAE integration

This candidate fills scMAE's robust latent-prior gap:

- The scMAE encoder, mask predictor, and masked expression reconstruction remain the main path.
- A `mu/logvar` VAE adapter provides latent uncertainty.
- A time-conditioned latent denoiser predicts diffusion noise added to the latent embedding.
- ZINB likelihood is enabled only when raw counts can be aligned to the selected HVG genes.

## Expression semantics

`scaled_expr` may be used only as encoder input. Masked expression reconstruction uses unscaled log-expression. ZINB uses aligned raw counts and size factors only; if raw counts cannot be aligned, ZINB is disabled and recorded in `preprocess_config.json` and `summary.json`.

## NeighborMix

NeighborMix is not used. The relationship is independent and complementary; `mixed_cell_fraction=0.0`.

## Outputs

Each run writes `embedding_final.npy`, `labels.npy`, `posterior_variance.npy`, `training_history.json`, `diagnostics.json`, `summary.json`, `args.json`, and optional fixed-k metrics. Smoke/screen rows are written only to `新模型独立快筛单次结果.csv` and summarized in `新模型独立快筛汇总结果.csv`.
