# rank41_nonstationary_gene_attention_full

Independent-full scMAE candidate based on **Non-stationary Transformers: Exploring the Stationarity in Time Series Forecasting**.

## Theory basis

The local report marks this method as a usable signal-modeling reference, but requires constructing a gene module/order first because genes have no natural timestamp order. This implementation uses sampled log-expression SVD loadings to order genes, then treats contiguous blocks as gene modules.

The borrowed mechanism is:

- Series Stationarization on each cell's gene-module sequence
- learned de-stationary factors `tau` and `delta` from the raw module sequence and its mean/std
- self-attention scores are rescaled as `QK * tau + delta`

## scMAE connection

The main objective remains scMAE:

- masked expression reconstruction on log-expression targets
- mask prediction BCE

The new core mechanism is a non-stationary gene-module attention adapter plus:

- module mean target regression
- latent-mask view consistency
- light regularization of `tau/delta` to avoid runaway attention factors

This addresses scMAE's **semantic target / non-stationary module-context** gap.

## Data semantics

- `scaled_expr`: optional encoder input only.
- `log_expr`: reconstruction and module-target source.
- raw counts are not used as NB/ZINB targets.
- labels are not used for training.

## NeighborMix relationship

NeighborMix is not used. The relationship is independent and potentially complementary. `mixed_cell_fraction=0.0`.

## Source notes

The indexed GitHub repository `https://github.com/thuml/Nonstationary_Transformers` was readable. The implementation inspected:

- `/tmp/Nonstationary_Transformers_repo/ns_models/ns_Transformer.py`
- `/tmp/Nonstationary_Transformers_repo/ns_layers/SelfAttention_Family.py`
- `/tmp/Nonstationary_Transformers_repo/ns_layers/Transformer_EncDec.py`

This directory rewrites those ideas for single-cell gene modules and keeps an independent model/loss/training loop.

Screen results are candidate evidence only and must not be appended to `全benchmark结果.csv`.
