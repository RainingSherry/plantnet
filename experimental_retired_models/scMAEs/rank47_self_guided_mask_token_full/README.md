# rank47_self_guided_mask_token_full

Independent-full scMAE candidate adapted from **Self-Guided Masked Autoencoder**.

## Method Basis

The local report recommends self-guided MAE as a masked modeling / discriminator-target method and suggests self-guided mask, gene-specific rank tokens, and curriculum-style adaptive mask. The paper finds that MAE learns token relations early and uses its internal progress to generate informed masks.

The GitHub URLs in the index/report were tested but were not anonymously readable in this environment, so this implementation is reconstructed from the PDF and local report.

## scMAE Gap Addressed

This candidate targets the **mask / semantic target** gap:

- scMAE mask prediction is retained.
- masked expression reconstruction is retained.
- Gene-specific log-expression quantile bins define rank-token targets.
- An EMA of per-gene reconstruction difficulty guides future mask probability.
- A difficulty head predicts gene difficulty as an internal self-guidance signal.

## Data Semantics

- `scaled_expr` is used only as encoder input when `--scale_input true`.
- `log_expr` is used as masked expression reconstruction target and quantile-token source.
- No count likelihood, NB/ZINB, or generated-cell evaluation is used.

## NeighborMix Relationship

NeighborMix is not used here. The method is independent and does not mix cell expressions. Therefore `mixed_cell_fraction=0.0`.

## Differences From Original Self-Guided MAE

- Image patches are replaced by genes.
- Patch clustering is mapped to gene-specific rank-token targets and reconstruction-difficulty-guided masks.
- The scMAE reconstruction and mask-prediction objectives remain the main training path.

## Screen Caveat

Smoke and screen results are candidate evidence only. They must not be appended to `全benchmark结果.csv` and are not formal performance claims.
