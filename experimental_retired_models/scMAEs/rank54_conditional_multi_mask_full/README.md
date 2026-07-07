# rank54_conditional_multi_mask_full

Independent-full scMAE candidate based on Conditional MAE / multiple masking.

## Theory basis

The PDF studies multiple masking stages in MAE and describes Conditional MAE, where later masking is conditioned on representations from earlier visible/unmasked content. No GitHub URL is provided in the index, so this implementation is reconstructed from the PDF and the scMAE improvement report.

## scMAE integration

This candidate fills scMAE's mask/robustness gap:

- First-stage mask corrupts expression input and retains scMAE mask prediction.
- Second-stage mask is applied to the intermediate encoder representation through a learned latent mask token.
- The decoder still reconstructs masked log-expression targets.
- A light condition head predicts the clean stage-1 representation to keep the latent masking stage anchored.

## Expression semantics

`scaled_expr` may be used only as encoder input. Masked expression reconstruction targets are unscaled log-expression. No count objective is computed from scaled expression.

## NeighborMix

NeighborMix is not used. The relation is independent and complementary; `mixed_cell_fraction=0.0`.

## Outputs

Each run writes `embedding_final.npy`, `labels.npy`, `training_history.json`, `diagnostics.json`, `summary.json`, `args.json`, and optional fixed-k metrics. Smoke/screen rows are written only to `新模型独立快筛单次结果.csv` and summarized in `新模型独立快筛汇总结果.csv`.
