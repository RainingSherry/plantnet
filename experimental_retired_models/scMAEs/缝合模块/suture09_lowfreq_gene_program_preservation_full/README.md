# suture09_lowfreq_gene_program_preservation_full

This independent-full candidate rewrites the Low-frequency Preservation idea for scRNA representation learning.

## Mechanism

The original reference module preserves low-frequency image content. This version does not use image FFT, wavelets, or any 2D reshaping of gene vectors.

Instead, it builds a coarse log-expression gene program:

1. sort genes by mean log expression;
2. split the sorted genes into fixed bins;
3. use each bin mean as a low-frequency gene-program target.

The model keeps the scMAE backbone with masked expression reconstruction and mask prediction. The new branch adds a weak gene-program prediction head and a very small gated latent adapter. Setting `--adapter_weight 0` and `--program_weight 0` recovers a plain scMAE-style path.

## Gap Addressed

This candidate targets semantic target and low-frequency preservation. It is meant to avoid the failure mode of overly strong graph/prototype injection by only asking the latent representation to retain coarse expression structure.

## NeighborMix

NeighborMix is not used. There is no cell mixing, and diagnostics report `mixed_cell_fraction=0.0`.

## Source

Reference module:

`/home/luolie/biopipeline/dimension-reduction/plantnet/缝合模块/即插即用/224 Low-frequency Preservation(TPAMI 2026).py`

Implementation source: rewritten from the mechanism description for scRNA data.
