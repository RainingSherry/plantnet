# rank35_videomae_gene_module_full

Independent-full scMAE candidate adapted from **VideoMAE: Masked Autoencoders are Data-Efficient Learners for Self-Supervised Video Pre-Training**.

## Method Basis

The local report recommends VideoMAE as a gene sequence / signal modeling candidate, with the warning that genes must not be treated as an arbitrary random sequence. The paper and official implementation emphasize three mechanisms used here: high-ratio tube masking, visible-token-only encoder, and a mask-token decoder that reconstructs masked targets.

The official GitHub implementation was inspected for `TubeMaskingGenerator`, `modeling_pretrain.py`, and `engine_for_pretraining.py`. It masks the same spatial map across the temporal axis, sends only visible tokens into the encoder, appends mask tokens in the decoder, and computes reconstruction loss on masked positions.

## scMAE Gap Addressed

This candidate targets the **mask / signal target** gap:

- scMAE mask prediction is retained.
- masked expression reconstruction is retained.
- genes are ordered by a deterministic greedy walk over absolute gene-gene correlations in `log_expr`;
- ordered genes are chunked into gene modules;
- VideoMAE-style tube/block masks hide correlated module positions across pseudo-frames;
- an asymmetric encoder-decoder reconstructs masked ordered gene modules.

## Data Semantics

- `scaled_expr` is used only as encoder input when `--scale_input true`.
- `log_expr` is used as masked expression reconstruction target.
- `gene_order.npy` records the data-derived ordering.
- No count likelihood, NB/ZINB, or token target is produced from scaled expression.

## NeighborMix Relationship

NeighborMix is not used. This method is independent and potentially complementary: it changes the masking geometry but never mixes cells. `mixed_cell_fraction=0.0`.

## Difference From Original VideoMAE

The original VideoMAE operates on video cube/tube tokens. This scRNA-seq implementation replaces video space-time patches with ordered gene modules. The pseudo-temporal structure is explicitly data-derived from gene correlations to avoid imposing arbitrary HVG index order.

## Screen Caveat

Smoke and screen results are candidate evidence only. They must not be appended to `全benchmark结果.csv` and are not formal performance claims.
