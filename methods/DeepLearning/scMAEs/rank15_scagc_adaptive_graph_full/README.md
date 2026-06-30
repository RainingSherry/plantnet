# Rank 15 scAGC adaptive graph scMAE

This folder independently adapts **scAGC: Learning Adaptive Cell Graphs with Contrastive Guidance for Single-Cell Clustering** to the scMAE benchmark protocol.

## Source and reproduction scope

- Paper/index item: Rank 15 `scAGC`.
- PDF: `../参考文献/01_PDF论文_按推荐程度排序/015_高_scAGC_Learning_Adaptive_Cell_Graphs_with_Contrastive_Guidance_for_Single-Cell_Clustering.pdf`.
- GitHub URL in index: none.
- PDF SHA-256 from the local index: `E63386A5EC7B85E74FE71C2102A6318FFFA9A08A0218C78BE9B161AABDAF5D44`.

## Implemented paper mechanisms

- Initial cell graph construction with KNN.
- Adaptive graph refinement with an RBF similarity matrix, Gumbel noise, temperature softmax, top-K straight-through estimator, and symmetrization.
- TAGCN-style graph encoder with 0..K hop kernels.
- Inner-product adjacency decoder and graph reconstruction loss.
- ZINB decoder and negative log-likelihood for nonnegative expression reconstruction.
- Temporal graph contrastive guidance between embeddings from the previous and adaptive graph.
- Student-t clustering distribution, sharpened target distribution, and KL clustering objective.

## scMAE-family adaptation

scAGC is a graph clustering autoencoder, not originally a masked autoencoder. To keep the scMAE-family protocol, this variant adds an auxiliary masked reconstruction branch with fixed semantics:

`1 = expression gene replaced for the auxiliary scMAE reconstruction branch`.

The encoder/evaluation input uses the same scaled HVG matrix as the rest of the benchmark. The ZINB target uses the same selected HVG genes with `scale_input=False` and nonnegative log-normalized expression, because the ZINB likelihood is defined on nonnegative counts/expression values.

## Scalability adaptation

The paper describes a graph over all cells. To keep the same runner usable on larger benchmark datasets, this implementation constructs and refines adaptive cell graphs inside each mini-batch. The core differentiable graph-learning mechanism is preserved, while avoiding a full `N x N` graph in memory.

## Deliberate exclusions

- No official repository or checkpoint is listed in the local index.
- No labels are used during training. Labels are used only by the shared KMeans evaluation routine after training.
