# rank48_mae_asymmetric_decoder_full

Independent-full scMAE candidate adapted from **Masked Autoencoders Are Scalable Vision Learners**.

## Method Basis

The local report recommends MAE-style masked modeling only when mapped to gene-expression structure. The official MAE repository was inspected; the relevant ideas are per-sample random masking, asymmetric encoder/decoder, a decoder embedding, mask tokens, and reconstruction loss on masked content.

## scMAE Gap Addressed

This candidate targets the **mask / decoder target** gap:

- scMAE random cell-swap corruption is retained.
- scMAE mask prediction is retained.
- masked expression reconstruction is retained.
- The decoder is changed into an asymmetric latent-to-decoder bottleneck with a learned mask-token bias.

## Data Semantics

- `scaled_expr` is used only as encoder input when `--scale_input true`.
- `log_expr` is used as reconstruction target.
- No count likelihood, NB/ZINB, tokenized counts, or generated-cell evaluation is used.

## NeighborMix Relationship

NeighborMix is not used here. The method is independent and does not mix cells beyond the original scMAE random swap corruption. Therefore `mixed_cell_fraction=0.0`.

## Differences From Original MAE

- Image patches are replaced by gene-expression vectors.
- The implementation avoids full ViT patch attention for screen-time feasibility.
- The asymmetry is implemented as a stronger decoder bottleneck around the original scMAE MLP encoder.

## Screen Caveat

Smoke and screen results are candidate evidence only. They must not be appended to `全benchmark结果.csv` and are not formal performance claims.
