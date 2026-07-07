# rank30_fuzzy_rough_quantification_full

Independent-full scMAE candidate adapted from **Fuzzy Rough Sets Based on Fuzzy Quantification**.

## Method Basis

The paper introduces fuzzy quantifier-based fuzzy rough sets, where lower approximation captures certain/core membership and upper approximation captures possible/boundary membership. The local index provides no GitHub URL, so this implementation is reconstructed from the PDF and the local scMAE improvement report.

## scMAE Gap Addressed

This candidate targets the **boundary / clustering head / robust loss** gap:

- scMAE mask prediction is retained.
- masked expression reconstruction is retained.
- fuzzy Student-t membership defines cluster membership;
- rough lower approximation strengthens only high-confidence core cells;
- rough upper approximation leaves uncertain boundary cells entropy-tolerant;
- rough width, balance, and center-separation terms reduce collapse.

## Data Semantics

- `scaled_expr` is used only as encoder input when `--scale_input true`.
- `log_expr` is used as masked expression reconstruction target.
- SVD anchor is computed from encoder input only as stabilizing bottleneck context.
- No count likelihood, NB/ZINB, token objective, or generated-cell evaluation is used.

## NeighborMix Relationship

NeighborMix is not used. This method is independent and potentially complementary. Since no cell mixing is performed, `mixed_cell_fraction=0.0`.

## Screen Caveat

Smoke and screen results are candidate evidence only. They must not be appended to `全benchmark结果.csv` and are not formal performance claims.
