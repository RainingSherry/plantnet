# PlantSPADE-LGCL Method Note

PlantSPADE-LGCL treats plant single-cell expression as sparse support geometry, not as another dense reconstruction task.

## Objects

Let `X in R^{C x G}` be the cell-gene expression matrix.

`X = M o A`

- `M`: support matrix, where `M_cg = 1` if expression of gene `g` is observed as non-zero in cell `c`, otherwise `0`.
- `A`: amplitude matrix, the non-zero expression magnitude after shared normalization and `log1p`.
- `S_c = {g | M_cg = 1}`: the observed support set of cell `c`.

The canonical data rule is:

- `M` comes from raw counts whenever possible: `layers["counts"]`, then `adata.raw.X`, then raw-looking `adata.X`.
- `A` comes from `normalize_total(target_sum=1e4) + log1p`.
- HVG selection is shared across methods, with `n_top_genes=2000` in the main table.

## Views

Local view: a cell-gene bipartite support graph built from `M`.

Global view: TF-IDF on `A`, followed by truncated SVD. This gives a low-rank cell view with gene-frequency correction.

## Model

The local encoder uses LightGCN-style normalized bipartite propagation:

- cell embeddings and gene embeddings are initialized as trainable lookup tables;
- normalized cell-gene support adjacency propagates embeddings across observed support edges;
- layer outputs are averaged to produce local cell and gene embeddings.

Training uses:

- BPR ranking loss on positive support edges versus sampled zero-expression genes;
- InfoNCE alignment between local support-graph cell embeddings and global TF-IDF-SVD cell embeddings;
- optional sparse module regularization.

## Negative Sampling

`negative_sampling.py` implements three modes:

- `random_zero`: uniformly sample genes absent from the current cell.
- `idf_weighted_zero`: preferentially sample globally common genes that are absent in the current cell. These are high-information zeros because many cells express them.
- `neighbor_conflict_zero`: preferentially sample genes expressed by similar cells but absent in the current cell.

The main experiment uses `random_zero`; the other modes are ablations.

## SupportGeneAttention

SupportGeneAttention is a sparse support-set readout. For each cell, attention is computed only over `S_c`; no dense `cell x gene` attention matrix is constructed.

The post-hoc version is the default because it is stable and does not change the main LGCL training objective. A lightweight trainable refiner is also available; only scalar coefficients such as `beta`, `gamma`, and `eta` participate in training.

Attention terms:

- cell-gene embedding similarity;
- optional amplitude term controlled by `beta`;
- optional gene IDF term controlled by `gamma`;
- residual readout strength controlled by `eta`.

Required attention ablations:

- `support_attention`
- `attention_no_idf`
- `attention_no_amplitude`
- `attention_topk_64`
- `attention_topk_128`
- `attention_topk_256`

`attention_no_idf` is a key mechanism ablation. The global view already applies TF-IDF before SVD, so adding IDF again in local support attention can double-count gene-frequency correction. This ablation tests whether local attention should rely on sparse support and amplitude without a second IDF term.

## Interpretation

Attention top genes are candidate explanation genes: genes assigned high weight by the sparse support-set readout. They should not be claimed as marker genes by themselves. The protocol reports overlap with Wilcoxon DEG markers so attention genes and classical differential-expression markers can be compared and cross-validated.
