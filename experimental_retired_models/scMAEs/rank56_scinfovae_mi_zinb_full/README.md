# rank56_scinfovae_mi_zinb_full

Independent-full scMAE candidate based on ScInfoVAE.

## Theory basis

ScInfoVAE combines InfoVAE mutual-information regularization with a ZINB observation model for scRNA-seq representation learning. The paper argues that ordinary VAEs can ignore latent variables when the decoder is too flexible, while InfoVAE/MMD encourages useful latent representations. No GitHub URL is provided in the index, so this implementation is reconstructed from the PDF and the scMAE improvement report.

## scMAE integration

This candidate fills scMAE's robust loss / uncertainty gap:

- The original scMAE encoder shape is retained.
- A lightweight `mu/logvar` adapter turns the hidden embedding into an InfoVAE latent variable.
- Mask prediction and masked expression reconstruction remain the primary objective.
- MMD and tiny KL weights regularize the latent distribution.
- ZINB likelihood is enabled only when raw counts can be aligned to the selected HVG genes.

## Expression semantics

`scaled_expr` may be used only as encoder input. Masked expression reconstruction uses unscaled log-expression. ZINB uses aligned raw counts and size factors only; if raw counts cannot be aligned, ZINB is disabled and recorded in `preprocess_config.json` and `summary.json`.

## NeighborMix

NeighborMix is not used. The relationship is independent and complementary; `mixed_cell_fraction=0.0`.

## Outputs

Each run writes `embedding_final.npy`, `labels.npy`, `posterior_variance.npy`, `training_history.json`, `diagnostics.json`, `summary.json`, `args.json`, and optional fixed-k metrics. Smoke/screen rows are written only to `新模型独立快筛单次结果.csv` and summarized in `新模型独立快筛汇总结果.csv`.
