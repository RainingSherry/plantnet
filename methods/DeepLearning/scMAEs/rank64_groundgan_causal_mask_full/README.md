# rank64_groundgan_causal_mask_full

Independent-full scMAE candidate based on GRouNdGAN.

## Theory basis

GRouNdGAN imposes a user-defined gene regulatory network inside a causal generator through sparse masks. Its generator first produces TF expressions, then uses a masked causal generator so each target gene is generated only from its regulating TFs and gene-specific noise. The implementation uses `MaskedLinear` so masked weights and masked gradients remain zero.

## scMAE integration

This candidate fills the graph / target dependency gap:

- scMAE mask prediction and masked expression reconstruction remain the primary loss.
- A pseudo TF pool is chosen from high-variance genes in unscaled log-expression.
- Each target gene receives top-k pseudo regulators by absolute coexpression.
- A causal dependency adapter predicts target expression from only its masked regulator set and the scMAE cell embedding.

## Expression semantics

`scaled_expr` may be used only as encoder input. Causal graph construction, causal targets, and reconstruction targets use unscaled log-expression. This candidate does not use generated pseudo-cells in evaluation.

## NeighborMix

NeighborMix is not used. The relationship is independent and complementary; `mixed_cell_fraction=0.0`.

## Source note

The GitHub repository `https://github.com/Emad-COMBINE-lab/GRouNdGAN` was cloned and inspected. Relevant source files include `src/networks/masked_causal_generator.py`, `src/layers/masked_linear.py`, and `src/preprocessing/grn_creation.py`.

## Outputs

Each run writes `embedding_final.npy`, `labels.npy`, `pseudo_tf_pool.npy`, `causal_regulators.npy`, `causal_regulator_weights.npy`, `training_history.json`, `diagnostics.json`, `summary.json`, `args.json`, and optional fixed-k metrics. Smoke/screen rows are written only to `新模型独立快筛单次结果.csv` and summarized in `新模型独立快筛汇总结果.csv`.
