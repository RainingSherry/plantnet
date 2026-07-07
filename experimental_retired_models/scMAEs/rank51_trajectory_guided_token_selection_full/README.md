# rank51_trajectory_guided_token_selection_full

Independent-full scMAE candidate based on rank51, "Reinforcement Learning meets Masked Video Modeling: Trajectory-Guided Adaptive Token Selection".

## Theory basis

The reference report classifies this paper as a policy-search/adaptive token selection method. The PDF describes a warmup phase with random masking followed by trajectory-guided token sampling, where token scores are derived from temporal motion/trajectory cues rather than a fixed random mask. No GitHub URL is provided in `02_整理索引.csv`, so this implementation is reconstructed from the report and PDF.

## scMAE integration

This candidate fills scMAE's mask/target gap. For scRNA-seq, the video token trajectory is translated into per-gene reconstruction dynamics:

- `level`: EMA of per-gene masked reconstruction difficulty.
- `velocity`: EMA of change in difficulty.
- `variance`: fixed log-expression variance.
- `zero_fraction`: fixed expression dropout/zero fraction.

After warmup, these features are passed through this directory's own `TrajectoryTokenSampler` to produce gene-token mask probabilities. The main model remains an scMAE-style encoder with mask prediction and masked expression reconstruction.

## Expression semantics

`scaled_expr` may be used only as encoder input. `log_expr` is used as the reconstruction target and as the source of gene trajectory statistics. No count or token objective is built from scaled expression.

## NeighborMix

NeighborMix is not used in this candidate. The relation is independent and complementary. `mixed_cell_fraction` is always `0.0`.

## Outputs

Each run saves `embedding_final.npy`, `labels.npy`, `training_history.json`, `diagnostics.json`, `summary.json`, `args.json`, and optional fixed-k metrics. Smoke/screen rows are appended only to `新模型独立快筛单次结果.csv` and summarized in `新模型独立快筛汇总结果.csv`; this candidate must not update `全benchmark结果.csv`.
