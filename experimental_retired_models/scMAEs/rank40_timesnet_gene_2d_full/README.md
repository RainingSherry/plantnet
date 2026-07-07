# rank40_timesnet_gene_2d_full

Independent-full scMAE candidate based on **TimesNet: Temporal 2D-Variation Modeling for General Time Series Analysis**.

## Theory basis

The local report marks TimesNet as a usable signal-modeling reference, with the same caveat as TS2Vec: genes have no natural timestamp order. This implementation therefore first builds a data-driven gene-module axis from log-expression SVD loadings, then applies TimesNet-style period discovery and 2D variation modeling on that module sequence.

TimesNet's transferable mechanism is:

- FFT-based top-period discovery on the module sequence
- reshape 1D module sequence into 2D tensors by period
- lightweight inception-style 2D kernels for intra-period and inter-period variation
- adaptive aggregation of multiple period views

## scMAE connection

The main objective remains scMAE:

- masked expression reconstruction on log-expression targets
- mask prediction BCE

TimesNet adds one core mechanism only:

- gene-module 2D variation adapter
- module mean target regression
- two latent-mask module views with SmoothL1 consistency

This addresses scMAE's **semantic target / long-range module context** gap.

## Data semantics

- `scaled_expr`: optional encoder input only.
- `log_expr`: reconstruction and module-target source.
- raw counts are not used as NB/ZINB targets.
- gene modules are data-driven and unsupervised; labels are used only for evaluation.

## NeighborMix relationship

NeighborMix is not used. The relationship is independent and potentially complementary. `mixed_cell_fraction=0.0`.

## Source notes

The indexed GitHub URL `https://github.com/thuml/TimesNet` currently contains README/images and states that complete code moved to `https://github.com/thuml/Time-Series-Library`. The implementation inspected:

- `/tmp/TimesNet_repo/README.md`
- `/tmp/Time-Series-Library_repo/models/TimesNet.py`
- `/tmp/Time-Series-Library_repo/layers/Conv_Blocks.py`

This directory rewrites those ideas for single-cell gene modules and keeps an independent model/loss/training loop.

Screen results are candidate evidence only and must not be appended to `全benchmark结果.csv`.
