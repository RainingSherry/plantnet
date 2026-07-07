# rank50_cmt_collaborative_mask_target_full

Independent-full scMAE candidate adapted from **The Dynamic Duo of Collaborative Masking and Target for Advanced Masked Autoencoder Learning**.

## Method Basis

The local report recommends mapping this method to gene-specific tokens, replaced-expression detection, and curriculum masking while preserving scMAE. The PDF proposes CMT-MAE, which linearly combines teacher and student attention maps for masking and uses teacher/student features as collaborative reconstruction targets.

No GitHub URL is listed in the index, so this implementation is reconstructed from the PDF and local report.

## scMAE Gap Addressed

This candidate targets the **mask / semantic target / teacher** gap:

- scMAE mask prediction is retained.
- masked expression reconstruction is retained.
- Fixed SVD features from log expression act as teacher targets.
- EMA student encoder features act as student targets.
- Teacher gene variance and EMA student reconstruction error jointly guide mask probabilities.

## Data Semantics

- `scaled_expr` is used only as encoder input when `--scale_input true`.
- `log_expr` is used as reconstruction target and fixed teacher feature source.
- No count likelihood, NB/ZINB, or generated-cell evaluation is used.

## NeighborMix Relationship

NeighborMix is not used here. The method is independent and does not mix cells beyond masked replacement corruption. Therefore `mixed_cell_fraction=0.0`.

## Differences From Original CMT-MAE

- Image-patch attentions are replaced by gene-level teacher/student saliency.
- CLIP teacher features are replaced by fixed SVD features from unlabeled log-expression data.
- Student targets use an EMA encoder to avoid label leakage.

## Screen Caveat

Smoke and screen results are candidate evidence only. They must not be appended to `全benchmark结果.csv` and are not formal performance claims.
