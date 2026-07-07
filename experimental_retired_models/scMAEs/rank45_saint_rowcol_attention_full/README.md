# rank45_saint_rowcol_attention_full

Independent-full scMAE candidate adapted from **SAINT: Improved Neural Networks for Tabular Data via Row Attention and Contrastive Pre-training**.

## Method Basis

The local report recommends SAINT for feature tokenizer and row-neighbor attention without replacing the biological mask path. The paper uses continuous feature embeddings, column self-attention, intersample row attention, and contrastive/denoising pretraining. The inspected official repository implements `RowColTransformer`, continuous-feature MLP embeddings, row/column attention modes, CutMix/mixup augmentations, and contrastive/denoising pretraining.

This candidate ports the architecture core to scRNA expression tables.

## scMAE Gap Addressed

This candidate targets the **feature tokenizer / row attention** gap:

- scMAE mask prediction is retained.
- masked expression reconstruction is retained.
- Gene expression is projected into gene-group feature tokens.
- Column attention models interactions among feature tokens.
- Batch row attention implements SAINT-style intersample attention.
- A second masked view supplies a contrastive row-consistency objective.

## Data Semantics

- `scaled_expr` is used only as encoder input when `--scale_input true`.
- `log_expr` is used as masked expression reconstruction target.
- No count likelihood, NB/ZINB, tokenized counts, or generated-cell evaluation is used.

## NeighborMix Relationship

NeighborMix is not used here. The method is independent and does not mix cell expressions. Therefore `mixed_cell_fraction=0.0`.

## Differences From Original SAINT

- Categorical/continuous tabular columns are replaced by gene-group expression tokens.
- CutMix and mixup are not used to avoid forced cell mixing.
- The model uses one lightweight row-column attention block for screen-time feasibility.

## Screen Caveat

Smoke and screen results are candidate evidence only. They must not be appended to `全benchmark结果.csv` and are not formal performance claims.
