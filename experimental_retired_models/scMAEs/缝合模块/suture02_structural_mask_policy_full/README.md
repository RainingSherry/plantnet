# suture02_structural_mask_policy_full

This independent-full candidate adapts the Structural-Aware Multi-Scale Masking Module idea into a one-dimensional gene-level mask policy for scMAE.

## Theory

The reference SMMM module uses structural saliency masks to filter redundant image skip features. This implementation does not use image convolutions and does not reshape gene vectors into 2D images. Instead, it estimates gene saliency from log-expression mean, variance, and dropout statistics, then masks likely redundant genes more often while protecting high-risk marker-like genes.

## scMAE Gap

The module targets the `mask` gap. Vanilla random masking can damage rare-cell or boundary marker genes. Here the mask distribution is learned and regularized so high variance / high dropout / low mean genes are less likely to be destructively zero-masked.

## Data Semantics

`scaled_expr` is used only as encoder input when `--scale_input true`. `log_expr` is the reconstruction target and the source of gene statistics. No count, token, NB/ZINB, or diffusion target is derived from scaled expression.

## NeighborMix

This candidate is independent of NeighborMix. It performs no cell mixing, so `mixed_cell_fraction` is always `0.0`.

## Notes

The module can be disabled with `--policy_weight 0`, which falls back to uniform scMAE-style masking. Screen results are written only to quick-screen CSV files and must not be treated as formal benchmark evidence.
