# rank43_graphgps_local_global_full

Independent-full scMAE candidate adapted from **GraphGPS: General Powerful Scalable Graph Transformers**.

## Method Basis

The local report lists GraphGPS as a graph / neighbor reliability method. The paper proposes a three-part recipe for graph Transformers: positional or structural encoding, local message passing, and global attention. The inspected official repository implements `GPSLayer` as a residual combination of local MPNN and Transformer-style global attention, with optional LapPE/RWSE-style encodings.

This candidate ports that recipe to single-cell KNN graphs without importing the full PyG GraphGPS stack.

## scMAE Gap Addressed

This candidate targets the **graph / neighbor reliability** gap:

- scMAE mask prediction is retained.
- masked expression reconstruction is retained.
- PCA-KNN over cells supplies local graph context.
- Lightweight structural encodings combine SVD coordinates, local density, and KNN in-degree proxy.
- A shallow GraphGPS adapter combines local neighbor messages with batch global multi-head attention.
- A drop-neighbor consistency term encourages stable graph-aware embeddings.

## Data Semantics

- `scaled_expr` is used only as encoder input when `--scale_input true`.
- `log_expr` is used as masked expression reconstruction target.
- No count likelihood, NB/ZINB, diffusion count target, or generated-cell evaluation is used.
- Graph targets and diagnostics are derived from KNN structure, not from scaled expression counts.

## NeighborMix Relationship

NeighborMix is not used here. This candidate is independent and potentially complementary: it estimates edge confidence and boundary risk, but does not mix cell expressions. Therefore `mixed_cell_fraction=0.0`.

## Differences From Original GraphGPS

- The original GraphGPS framework uses PyG and graph-level datasets; this implementation is a single-cell scMAE adapter.
- Full Laplacian eigenvector encodings are replaced by scalable SVD/local-density/KNN-degree structural encodings.
- The local/global stack is intentionally shallow and residual to reduce oversmoothing on cell KNN graphs.
- Global attention is mini-batch attention for scalability.

## Screen Caveat

Smoke and screen results are candidate evidence only. They must not be appended to `全benchmark结果.csv` and are not formal performance claims.
