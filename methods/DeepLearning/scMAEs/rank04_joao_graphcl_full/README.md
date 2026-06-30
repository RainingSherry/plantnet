# Rank 04: JOAO GraphCL scMAE

Source paper: **JOAO: Automated Data Augmentations for Graph Contrastive Learning**
(ICML 2021).

This independent variant adapts JOAO to scMAE by treating each training mini-batch
as a cell graph. Nodes are cells, node features are preprocessed expression
vectors, and edges are batch-local KNN similarities. A GCN encoder learns two
augmented graph views with NT-Xent contrastive loss, while a masked-expression
decoder preserves the scMAE reconstruction objective.

JOAO-specific parts implemented here:

- five graph/expression augmentations: edge drop, gene mask, expression noise,
  cell feature drop, and subgraph-style cell keep/drop;
- augmentation pair sampling from a learned `aug_prob`;
- periodic per-augmentation loss evaluation;
- JOAO probability update `P <- ProjSimplex(P + beta * (loss_aug - gamma *
  (P - uniform)))`, following the official code structure;
- graph contrastive objective over cell embeddings.

Mask semantics: `1 = expression gene masked in scMAE reconstruction branch`.
The reconstruction denominator is the number of masked gene entries in the batch,
clamped only for degenerate smoke tests.

Fair protocol shared with the other independent variants:

- repository `scMAE_family` preprocessing and fixed KMeans known-k evaluation;
- default `n_top_genes=1000`, `target_sum=10000`, `scale_input=True`;
- no dependence on old `scMAEs/common/model.py`;
- no GPU 0 or 7 should be used by the benchmark runner.

Not reproduced:

- original molecular/TU dataset loaders and PyG GIN implementations are not used
  directly because the input object here is a single-cell expression matrix, not a
  dataset of small molecular graphs;
- the core JOAO algorithmic ingredients retained are GraphCL contrastive learning,
  graph augmentations, augmentation probability search, and simplex projection.
