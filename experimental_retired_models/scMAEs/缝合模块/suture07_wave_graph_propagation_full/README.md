# suture07_wave_graph_propagation_full

This independent-full candidate adapts Wave Propagation Operator into a graph wave auxiliary branch for scMAE.

The reference module solves a wave equation in 2D image frequency space. This implementation does not use image DCT/FFT and does not reshape gene vectors into images. Instead, it builds a KNN graph on unscaled log-expression, computes first- and second-order graph propagation, and feeds the resulting wave context into a weak latent adapter.

## scMAE integration

The primary objective remains scMAE mask prediction and masked expression reconstruction. The wave branch is an auxiliary context and can be disabled with `--wave_weight 0`.

`scaled_expr` is used only as encoder input. KNN graph construction and reconstruction targets use `log_expr`. No cell mixing is used.

## Defaults

- `wave_k=15`
- `wave_alpha=0.55`
- `wave_damping=0.25`
- `wave_weight=0.12`

## NeighborMix

Independent and complementary. No cells are mixed, so `mixed_cell_fraction=0.0`.
