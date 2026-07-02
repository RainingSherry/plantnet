# rank28_graphormer_structural_bias_full

Independent-full scMAE candidate adapted from **Graphormer**.

## Method Basis

Graphormer shows that Transformers can model graphs when structural information is encoded into the attention mechanism. The paper proposes degree-based centrality encoding, shortest-path spatial encoding as attention bias, and edge encoding as an additional attention bias. The official GitHub implementation confirms that `attn_bias` is added directly to attention logits and that preprocessing provides `spatial_pos`, degree, and edge inputs.

## scMAE Gap Addressed

This candidate targets the **graph / neighbor reliability** gap:

- scMAE mask prediction is retained.
- masked expression reconstruction is retained.
- each cell forms a local KNN graph with itself as the center token and neighbors as local tokens;
- centrality encoding distinguishes the center cell and high-rank neighbors;
- KNN-rank spatial encoding and edge/drop-edge type encoding are used as attention bias;
- a light edge-confidence head provides neighbor reliability diagnostics;
- a scaled-expression SVD anchor is used only as encoder-side stabilizing context, not as count or token target.

## Data Semantics

- `scaled_expr` is used only as encoder input when `--scale_input true`.
- `log_expr` is used as masked expression reconstruction target.
- No NB/ZINB count likelihood or generated-cell evaluation is used.

## NeighborMix Relationship

NeighborMix is not used. This candidate is independent and potentially complementary: it estimates edge confidence and boundary risk but never mixes cell expressions. `mixed_cell_fraction=0.0`.

## Screen Caveat

Smoke and screen results are candidate evidence only. They must not be appended to `全benchmark结果.csv` and are not formal performance claims.
