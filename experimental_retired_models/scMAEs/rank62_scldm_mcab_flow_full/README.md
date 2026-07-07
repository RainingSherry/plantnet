# rank62_scldm_mcab_flow_full

Independent-full scMAE candidate based on scLDM.

## Theory basis

scLDM proposes exchangeable gene-expression modeling with a Multi-head Cross-Attention Block (MCAB) for permutation-invariant latent pooling and permutation-equivariant decoding, plus latent flow/diffusion modeling with linear interpolants and Diffusion Transformers. The local improvement report recommends reusing the lightweight denoising/latent regularization idea rather than training a full generator.

## scMAE integration

This candidate fills the semantic target / exchangeable gene representation gap:

- scMAE mask prediction and masked expression reconstruction remain the primary loss.
- A small MCAB-style gene-token cross-attention adapter pools expression tokens into the scMAE latent.
- A linear-interpolant flow head predicts the vector from Gaussian noise to the current latent.
- No generated cells are used for evaluation.

## Expression semantics

`scaled_expr` may be used only as encoder input. Reconstruction targets use unscaled log-expression. No count likelihood is used in this candidate.

## NeighborMix

NeighborMix is not used. The relationship is independent and complementary; `mixed_cell_fraction=0.0`.

## Source note

The index lists `https://github.com/czi-ai/scldm`, which was cloned and inspected. The report URL `https://github.com/OmicsML/scLDM` required credentials in this environment. Relevant inspected files include `README.md`, `src/scldm/layers.py`, `src/scldm/nnets.py`, `src/scldm/vae.py`, and `src/scldm/transport/path.py`.

## Outputs

Each run writes `embedding_final.npy`, `labels.npy`, `training_history.json`, `diagnostics.json`, `summary.json`, `args.json`, and optional fixed-k metrics. Smoke/screen rows are written only to `新模型独立快筛单次结果.csv` and summarized in `新模型独立快筛汇总结果.csv`.
