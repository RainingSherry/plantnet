# suture10_spatial_channel_gene_modulator_full

This independent-full candidate rewrites Spatial-Channel Feature Modulation for scRNA latent vectors.

## Mechanism

The image module uses spatial and channel attention maps. This implementation does not reshape gene vectors into 2D images. It uses:

- a latent channel gate, analogous to channel attention;
- a cell/sample gate, analogous to spatial saliency;
- a small residual latent modulation controlled by `--modulator_weight`.

The scMAE backbone keeps masked expression reconstruction and mask prediction. Setting `--modulator_weight 0` returns the model to a plain scMAE-style latent path.

## Gap Addressed

This candidate targets local latent calibration without graph propagation. It is intended as a safer alternative after prototype and graph-wave candidates damaged Melanoma/Macosko geometry.

## NeighborMix

NeighborMix is not used. There is no cell mixing, and diagnostics report `mixed_cell_fraction=0.0`.

## Source

Reference module:

`/home/luolie/biopipeline/dimension-reduction/plantnet/缝合模块/即插即用/227 Spatial-Channel Feature Modulator (CVPR 2026).py`

Implementation source: rewritten from the mechanism description for vector latent scRNA data.
