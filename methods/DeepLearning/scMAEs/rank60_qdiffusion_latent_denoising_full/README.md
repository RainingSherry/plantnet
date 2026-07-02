# rank60_qdiffusion_latent_denoising_full

Independent-full scMAE candidate based on q-Diffusion.

## Theory basis

q-Diffusion defines q-deformed kernels with adaptive nearest-neighbor bandwidths to retain high-dimensional gene coexpression geometry. Its Julia implementation exposes `DeformedKernel`, `q_exp`, `q_sum`, and `AdaptiveKernel`; the improvement report recommends using this as a lightweight denoising auxiliary rather than training a heavy generative model.

## scMAE integration

This candidate fills the robust target / semantic denoising gap:

- scMAE mask prediction and masked expression reconstruction remain the primary loss.
- A q-deformed adaptive KNN kernel is built from unscaled log-expression.
- The model adds a time-conditioned latent denoiser.
- Training corrupts the latent with a q-shaped noise schedule and denoises toward the clean scMAE latent, with a weak q-neighbor context consistency term.

## Expression semantics

`scaled_expr` may be used only as encoder input. q-kernel construction and reconstruction targets use unscaled log-expression. No NB/ZINB count likelihood is used.

## NeighborMix

NeighborMix is not used. The relationship is independent and complementary; q-kernel context could later define confidence gates, but this implementation never mixes cells and reports `mixed_cell_fraction=0.0`.

## Source note

The improvement report listed `https://github.com/zsteve/qDiffusion`, which was not anonymously readable in this environment. The index lists `https://github.com/marmarelis/QDiffusion.jl`, which was cloned and inspected. This Python candidate ports only the kernel/noise-schedule idea needed for scMAE adaptation.

## Outputs

Each run writes `embedding_final.npy`, `labels.npy`, `qdiffusion_neighbors.npy`, `qdiffusion_weights.npy`, `training_history.json`, `diagnostics.json`, `summary.json`, `args.json`, and optional fixed-k metrics. Smoke/screen rows are written only to `新模型独立快筛单次结果.csv` and summarized in `新模型独立快筛汇总结果.csv`.
