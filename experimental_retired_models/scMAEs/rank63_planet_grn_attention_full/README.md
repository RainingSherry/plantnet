# rank63_planet_grn_attention_full

Independent-full scMAE candidate based on Planet, an attention-guided probabilistic diffusion model for cell-type-specific GRN generation.

## Theory basis

The source paper frames gene regulatory network inference as a diffusion-style edge denoising problem conditioned on gene expression profiles. Planet uses a Triple Hybrid-Attention Transformer to combine graph attention, cross-attention, self-attention, and time-step embeddings so high-confidence regulatory edges are preserved while noisy edges are pruned.

## scMAE integration

This candidate fills the graph / semantic target gap:

- scMAE mask prediction and masked expression reconstruction remain the primary loss.
- A sparse gene-gene regulatory prior is reconstructed from unscaled log-expression coexpression.
- A time-guided gene-pair attention adapter predicts clean regulatory edge labels from scMAE cell embeddings.
- The GRN module is an auxiliary denoising target; generated GRN edges are not used as clustering labels.

## Expression semantics

`scaled_expr` may be used only as encoder input. Reconstruction and GRN prior construction use unscaled log-expression. No NB/ZINB count likelihood or generated pseudo-cell evaluation is used.

## NeighborMix

NeighborMix is not used. The relationship is independent and complementary; future work could use high-confidence GRN edges to guide gene masking, but this implementation does not mix cells and reports `mixed_cell_fraction=0.0`.

## Source note

The index provides no GitHub URL for rank63. This implementation is reconstructed from the PDF and the scMAE improvement report.

## Outputs

Each run writes `embedding_final.npy`, `labels.npy`, `grn_edges.npy`, `grn_edge_targets.npy`, `training_history.json`, `diagnostics.json`, `summary.json`, `args.json`, and optional fixed-k metrics. Smoke/screen rows are written only to `新模型独立快筛单次结果.csv` and summarized in `新模型独立快筛汇总结果.csv`.
