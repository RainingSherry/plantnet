# rank27_graph_barlow_twins_full

Independent-full scMAE candidate adapted from **Graph Barlow Twins: A self-supervised representation learning framework for graphs**.

## Method Basis

Graph Barlow Twins trains one symmetric graph encoder on two distorted graph views. Its loss normalizes both embedding batches, computes their cross-correlation matrix, pushes diagonal entries toward one, and pushes off-diagonal entries toward zero. The official GitHub implementation confirms feature masking, edge dropping, and the Barlow Twins loss without negative samples, EMA target encoders, or predictor asymmetry.

## scMAE Gap Addressed

This candidate targets the **graph / neighbor reliability / robust loss** gap:

- scMAE mask prediction is retained.
- masked expression reconstruction is retained.
- PCA-KNN supplies a shallow cell graph.
- Two graph views use feature masking and DropEdge.
- A symmetric shallow graph adapter avoids deep GNN oversmoothing.
- The Graph Barlow Twins loss adds invariance plus redundancy reduction.
- A light edge-confidence head records neighbor reliability diagnostics.

## Data Semantics

- `scaled_expr` is used only as encoder input when `--scale_input true`.
- `log_expr` is used as masked expression reconstruction target.
- No count likelihood, NB/ZINB, token, or generated-cell objective is used.

## NeighborMix Relationship

NeighborMix is not used. This candidate is independent and potentially complementary: it estimates edge confidence and boundary risk, but never mixes cell expressions. `mixed_cell_fraction=0.0`.

## Screen Caveat

Smoke and screen results are candidate evidence only. They must not be appended to `全benchmark结果.csv` and are not formal performance claims.
