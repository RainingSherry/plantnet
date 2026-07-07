# rank39_ts2vec_gene_module_full

Independent-full scMAE candidate based on **TS2Vec: Towards Universal Representation of Time Series**.

## Theory basis

The local report marks TS2Vec as a usable signal-modeling reference, but warns that single-cell genes have no natural timestamp order. This implementation therefore does not run TS2Vec over raw HVG index order. It first builds a data-driven gene axis by sampling log-expression cells, standardizing genes, taking SVD loadings, and sorting genes by the first-two loading angle. Contiguous blocks on that axis become gene modules.

The borrowed TS2Vec mechanism is hierarchical contextual consistency: two overlapping module-context views are encoded after latent timestamp masking, then contrasted at multiple pooling scales using instance-wise and module-position-wise contrast.

## scMAE connection

The scMAE body remains the main training objective:

- masked expression reconstruction on log-expression targets
- mask prediction BCE

TS2Vec adds one core mechanism only:

- gene-module sequence adapter with latent timestamp masking
- module mean target regression
- hierarchical contextual consistency over overlapping module views

This addresses scMAE's **target / semantic signal** gap by making the encoder respect co-expression-derived module context instead of only independent gene reconstruction.

## Data semantics

- `scaled_expr`: optional encoder input only.
- `log_expr`: masked expression reconstruction target and module target source.
- raw counts are not used as NB/ZINB targets.
- gene modules are derived from log-expression statistics, not labels.

## NeighborMix relationship

NeighborMix is not used here. The relationship is independent and potentially complementary. `mixed_cell_fraction` is always `0.0`.

## Source notes

The GitHub repository `https://github.com/yuezhihan/ts2vec` was readable. The local implementation follows the paper/code ideas for latent timestamp masking, overlapping context views, and hierarchical contrastive loss, but rewrites them for gene-module expression data and keeps an independent scMAE training loop.

Screen results are candidate evidence only and must not be appended to `全benchmark结果.csv`.
