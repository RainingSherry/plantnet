# rank26_bgrl_graph_bootstrap_full

Independent-full scMAE candidate adapted from **BGRL: Large-Scale Representation Learning on Graphs via Bootstrapping**.

## Method Basis

BGRL learns graph node representations by predicting a target encoder's representation of a second augmented graph view. The target encoder is updated by EMA, and the objective uses cosine prediction without negative examples. The official GitHub implementation confirms:

- separate online and target encoders;
- a predictor on the online branch;
- target encoder parameters updated by momentum;
- simple feature masking and edge masking augmentations.

## scMAE Gap Addressed

This candidate targets the **graph / neighbor reliability / teacher** gap:

- scMAE mask prediction is retained.
- masked expression reconstruction is retained.
- PCA-KNN over cells supplies a shallow graph context.
- A residual graph adapter uses one-hop neighbor summaries only, avoiding deep GNN oversmoothing.
- Two graph views use feature dropout and DropEdge.
- A BGRL EMA target branch adds graph bootstrap consistency.
- An edge-confidence head provides edge reliability diagnostics and a light positive-vs-random edge loss.

## Data Semantics

- `scaled_expr` is used only as encoder input when `--scale_input true`.
- `log_expr` is used as masked expression reconstruction target.
- No count likelihood, NB/ZINB, or generated-cell evaluation is used.

## NeighborMix Relationship

NeighborMix is not used here. The method is independent and potentially complementary: it estimates `edge_confidence` and boundary risk, but does not mix cell expressions. Therefore `mixed_cell_fraction=0.0`.

## Screen Caveat

Smoke and screen results are candidate evidence only. They must not be appended to `全benchmark结果.csv` and are not formal performance claims.
