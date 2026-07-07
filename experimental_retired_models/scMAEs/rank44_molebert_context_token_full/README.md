# rank44_molebert_context_token_full

Independent-full scMAE candidate adapted from **Mole-BERT: Rethinking Pre-training Graph Neural Networks for Molecules**.

## Method Basis

The local report categorizes Mole-BERT as a graph / neighbor reliability method. The paper argues that masked node modeling benefits from a context-aware discrete tokenizer rather than a small, imbalanced raw vocabulary. The inspected official repository includes `VectorQuantizer`, `DiscreteGNN`, `MaskAtom`, masked atom modeling, masked edge prediction, and triplet masked contrastive learning.

This candidate ports the context-aware tokenizer idea to single-cell data: each cell receives an unsupervised discrete context token from SVD cell features and KNN-neighbor context, then scMAE predicts that token while preserving masked expression reconstruction.

## scMAE Gap Addressed

This candidate targets the **semantic target / graph-neighbor context** gap:

- scMAE mask prediction is retained.
- masked expression reconstruction is retained.
- KNN context supplies local cell-neighborhood semantics.
- MiniBatchKMeans over cell-plus-context features forms discrete context tokens.
- A context-gated encoder predicts these tokens as a Mole-BERT-style masked node modeling objective.

## Data Semantics

- `scaled_expr` is used only as encoder input when `--scale_input true`.
- `log_expr` is used as masked expression reconstruction target.
- No NB/ZINB or count likelihood is used.
- Context tokens are unsupervised graph-context labels, not class labels.

## NeighborMix Relationship

NeighborMix is not used here. The method is independent and complementary: it learns context tokens and edge confidence, but does not mix cell expressions. Therefore `mixed_cell_fraction=0.0`.

## Differences From Original Mole-BERT

- Molecular atom tokens are replaced by unsupervised cell-context tokens.
- The full GNN/VQ-VAE pretraining stack is not imported; this is a lightweight scMAE-specific reconstruction from the paper and code ideas.
- TMCL is not used in this candidate to avoid stacking multiple paper mechanisms.

## Screen Caveat

Smoke and screen results are candidate evidence only. They must not be appended to `全benchmark结果.csv` and are not formal performance claims.
