# rank59_soft_graph_clustering_full

Independent-full scMAE candidate based on Soft Graph Clustering for single-cell RNA Sequencing Data.

## Theory basis

The source paper argues that hard binary cell graphs can amplify wrong inter-cluster edges in noisy scRNA-seq data. It proposes continuous soft edge weights from complementary similarity views and graph-aware clustering objectives. The local improvement report recommends using this family as a shallow edge-confidence and graph-reliability module rather than stacking deep GNN layers.

## scMAE integration

This candidate fills the graph / neighbor reliability gap:

- scMAE mask prediction and masked expression reconstruction remain the primary loss.
- A soft KNN graph is built on unscaled log-expression targets with Euclidean and cosine channels.
- The model adds a shallow pairwise edge adapter over scMAE embeddings.
- Training adds soft edge reconstruction, DropEdge robustness, and residual neighbor consistency.

## Expression semantics

`scaled_expr` may be used only as encoder input. Graph construction and reconstruction targets use unscaled log-expression. No count distribution is used in this candidate.

## NeighborMix

NeighborMix is not used. The relationship is independent and complementary: the learned edge confidence could later gate NeighborMix, but this implementation does not mix cells and reports `mixed_cell_fraction=0.0`.

## Implementation note

The index lists `https://github.com/seandavi/awesome-single-cell`, which is a general resource collection and does not contain a runnable scSGC implementation. This directory therefore reconstructs the usable scMAE adapter from the paper and the improvement report.

## Outputs

Each run writes `embedding_final.npy`, `labels.npy`, `soft_graph_neighbors.npy`, `soft_graph_weights.npy`, `edge_confidence.npy`, `training_history.json`, `diagnostics.json`, `summary.json`, `args.json`, and optional fixed-k metrics. Smoke/screen rows are written only to `新模型独立快筛单次结果.csv` and summarized in `新模型独立快筛汇总结果.csv`.
