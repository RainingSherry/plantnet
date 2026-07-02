# rank42_maskgae_edge_reconstruct_full

Independent-full scMAE candidate adapted from **MaskGAE: Masked Graph Modeling Meets Graph Autoencoders**.

## Method Basis

The local report recommends MaskGAE as a graph and neighbor-reliability module for scMAE. The paper frames masked graph modeling as masking graph structure and reconstructing missing edges. The inspected official repository implements edge masking, path masking, an encoder, an edge decoder, negative edge sampling, and a degree decoder.

This candidate uses the MaskGAE edge-wise variant because it maps cleanly to cell KNN graphs without adding a deep GNN stack.

## scMAE Gap Addressed

This candidate targets the **graph / neighbor reliability** gap:

- scMAE mask prediction is retained.
- masked expression reconstruction is retained.
- PCA-KNN over cells supplies local graph context.
- A shallow residual graph encoder aggregates only visible neighbors.
- Masked KNN edges are reconstructed with a pairwise edge decoder.
- A degree decoder regularizes how many local edges were hidden.
- A neighbor latent consistency term encourages local smoothness without forced cell mixing.

## Data Semantics

- `scaled_expr` is used only as encoder input when `--scale_input true`.
- `log_expr` is used as masked expression reconstruction target.
- No count likelihood, NB/ZINB, ZINB, diffusion count target, or generated-cell evaluation is used.
- Edge targets are graph structure labels from the PCA-KNN graph, not expression-count labels.

## NeighborMix Relationship

NeighborMix is not used here. The method is independent and potentially complementary: it estimates local edge reliability and boundary risk, but does not mix cell expressions. Therefore `mixed_cell_fraction=0.0`.

## Differences From Original MaskGAE

- The original paper targets graph representation learning; this implementation wraps the idea around scMAE.
- The encoder is intentionally shallow and residual to reduce oversmoothing risk on single-cell KNN graphs.
- Path-wise masking is not used in this candidate; only edge-wise masking is used.
- Masked expression reconstruction remains the main scMAE objective.

## Screen Caveat

Smoke and screen results are candidate evidence only. They must not be appended to `全benchmark结果.csv` and are not formal performance claims.
