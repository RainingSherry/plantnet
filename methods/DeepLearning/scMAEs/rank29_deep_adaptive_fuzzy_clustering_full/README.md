# rank29_deep_adaptive_fuzzy_clustering_full

Independent-full scMAE candidate adapted from **Deep Adaptive Fuzzy Clustering for Evolutionary Unsupervised Representation Learning**.

## Method Basis

The paper proposes joint reconstruction and fuzzy clustering in a bottleneck representation space. It emphasizes fuzzy membership, weighted adaptive entropy, and iterative representation refinement. No GitHub URL is provided in the local index, so this implementation is reconstructed from the PDF and the local scMAE improvement report.

## scMAE Gap Addressed

This candidate targets the **boundary / clustering head / robust loss** gap:

- scMAE mask prediction is retained.
- masked expression reconstruction is retained.
- a fuzzy Student-t membership head assigns each cell to known-k cluster centers;
- cluster centers are initialized from SVD anchor KMeans;
- fuzzy KL is activated only after warmup and only for high-confidence core cells;
- low-confidence boundary cells are allowed higher membership entropy;
- balance and center-separation terms reduce cluster collapse.

## Data Semantics

- `scaled_expr` is used only as encoder input when `--scale_input true`.
- `log_expr` is used as masked expression reconstruction target.
- SVD anchor is computed from encoder input only as stabilizing bottleneck context.
- No count likelihood, NB/ZINB, token objective, or generated-cell evaluation is used.

## NeighborMix Relationship

NeighborMix is not used. The method is independent and potentially complementary. Since no cell mixing is performed, `mixed_cell_fraction=0.0`.

## Screen Caveat

Smoke and screen results are candidate evidence only. They must not be appended to `全benchmark结果.csv` and are not formal performance claims.
