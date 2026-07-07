# suture08_parameter_free_edge_structure_full

This independent-full candidate adapts the PFESA idea from parameter-free edge/structure attention into a single-cell graph setting.

## Mechanism

- Input branch separation:
  - `scaled_expr`: encoder input only.
  - `log_expr`: reconstruction target and graph source.
- Graph low-frequency structure:
  - local KNN smoothing on log-expression for small/medium datasets;
  - anchor KNN fallback for large datasets such as Macosko.
- Graph high-frequency edge residual:
  - `edge = log_expr - neighbor_structure`.
- Parameter-free reliability:
  - computed from low-frequency structure energy versus high-frequency residual energy;
  - no learned attention parameters are used for the reliability score.

The reliability score controls a lightweight latent context adapter and a weak structure consistency loss. It does not replace the scMAE backbone, and `--adapter_weight 0` recovers a plain scMAE-style path.

## Relation To scMAE

The model keeps both required scMAE objectives:

- mask prediction;
- masked expression reconstruction.

The added mechanism targets the cluster-geometry and boundary-stability gap: reliable low-frequency cells can receive small graph-context denoising, while high-edge/boundary cells keep more direct reconstruction pressure.

## Relation To NeighborMix

NeighborMix is not used. This candidate performs no cell mixing, and diagnostics always report `mixed_cell_fraction=0.0`.

## Source

Reference module:

`/home/luolie/biopipeline/dimension-reduction/plantnet/缝合模块/即插即用/223 Parameter-Free Edge-Structure Attention(MICCAI 2026).py`

The original module is FFT-based for medical image segmentation. This implementation does not reshape gene vectors into images and does not use image FFT. It reconstructs the edge/structure idea on the log-expression KNN graph.
