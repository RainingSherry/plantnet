# Rank 14 CICL iterative contrastive scMAE

This folder independently adapts **CICL: scRNA-seq Data Clustering by Cluster-aware Iterative Contrastive Learning** to the scMAE benchmark protocol.

## Source and reproduction scope

- GitHub source: `https://github.com/Alunethy/CIRCLE`.
- Local source snapshot: `../external_sources/rank14_circle`, commit `aba15bc81dc0b4999f56c7a82dd5f13bf109f27c`.
- PDF: `../参考文献/01_PDF论文_按推荐程度排序/014_高_CICL_scRNA-seq_Data_Clustering_by_Cluster-aware_Iterative_Contrastive_Learning.pdf`.
- Referenced source files:
  - `training.py`
  - `trans_model.py`
  - `learner/contrastive_utils.py`
  - `learner/cluster_utils.py`

## Implemented paper mechanisms

- Transformer encoder over gene-expression patches.
- Two Gaussian/dropout augmented expression views.
- KMeans centroids updated iteratively from current embeddings.
- Student-t clustering head for pseudo-label probabilities.
- Projection head for contrastive learning.
- Instance-wise InfoNCE loss.
- Cluster-aware contrastive loss using same-pseudo-label positives and different-pseudo-label negatives.
- Small DEC-style target distribution KL term to stabilize the clustering head.

## scMAE-family adaptation

CICL itself is not a masked autoencoder. To keep a fair scMAE-family protocol, this variant adds an auxiliary masked reconstruction branch with explicit mask semantics:

`1 = expression gene replaced for the auxiliary scMAE reconstruction branch`.

The clustering and contrastive mechanisms remain the primary CICL contribution.

## Deliberate exclusions

- The original repository uses distributed training and dataset-specific loaders; this folder uses the shared scMAE data loader and single-process benchmark runner.
- No labels are used during training. Labels are used only by the shared KMeans evaluation routine after training.
