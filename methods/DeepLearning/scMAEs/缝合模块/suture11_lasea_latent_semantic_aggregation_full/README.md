# suture11_lasea_latent_semantic_aggregation_full

This independent-full candidate rewrites LaSEA for scRNA latent representations.

## Mechanism

The original image module uses multi-dilation convolutions and random pooling channel attention. This implementation does not reshape gene vectors into images. It applies multi-dilation 1D convolution across the learned latent vector, followed by latent attention and a small residual semantic update.

During training, the masked-view semantic summary is weakly aligned to a detached clean-view semantic summary. This preserves the scMAE masked reconstruction task while adding a semantic consistency target.

## Gap Addressed

This candidate targets semantic target and rare/boundary stability. It avoids graph propagation, prototype memory, and NeighborMix after those mechanisms showed dataset-specific instability.

## NeighborMix

NeighborMix is not used. There is no cell mixing, and diagnostics report `mixed_cell_fraction=0.0`.

## Source

Reference module:

`/home/luolie/biopipeline/dimension-reduction/plantnet/缝合模块/即插即用/237 Latent-Aware Semantic Extraction and Aggregation(2026 一区TOP).py`

Implementation source: rewritten from the mechanism description for vector latent scRNA data.
