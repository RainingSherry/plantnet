# Rank 13 mask-sc single-cell clustering

This folder independently adapts **Masked Modeling for Single-cell Clustering of scRNA-seq Data** to the scMAE benchmark protocol.

## Source and reproduction scope

- Paper/index item: Rank 13 `Masked Modeling for Single-cell Clustering of scRNA-seq Data`.
- PDF: `../参考文献/01_PDF论文_按推荐程度排序/013_高_Masked_Modeling_for_Single-cell_Clustering_of_scRNA-seq_Data.pdf`.
- GitHub URL: none listed in `02_整理索引.csv`.
- PDF SHA-256 from the local index: `07A1AC9FF48A3521C627E9093309CC2820BEB443701EDB9B514E3C3699D06D53`.

## Implemented paper mechanisms

- Input expression vectors are padded into a square 2D expression matrix.
- A 2D convolutional patch embedding layer creates matrix patch tokens.
- Random masking removes a fixed ratio of patch tokens from the encoder input.
- The encoder is a self-attention Transformer over the visible expression matrix patches.
- A sequence-guided decoder receives encoded visible tokens plus learned mask tokens and predicts masked patch target features.
- The target features come from a sequence-level encoder that is first trained with a contrastive-sc style self-supervised objective on short expression patches, then frozen before mask-sc training.
- Inference uses mean pooling over encoded patch-level features, matching the paper's clustering representation.

## Benchmark adaptation

- The benchmark still uses the shared fair preprocessing: `n_top_genes`, `target_sum`, and `scale_input` are controlled by `run.py` arguments.
- No labels are used during target pretraining or mask-sc training.
- KMeans evaluation follows the same known-`k` protocol as the other independent variants.
- Mask semantics are fixed as `1 = expression matrix patch removed from encoder and predicted by sequence-guided decoder`.

## Deliberate exclusions

- The original paper used a much larger Transformer configuration and 1600 epochs. This implementation keeps the core architecture while exposing depth/width/epoch arguments for the local benchmark.
- No external target encoder weights are loaded, because no official repository or checkpoint is listed in the index.
