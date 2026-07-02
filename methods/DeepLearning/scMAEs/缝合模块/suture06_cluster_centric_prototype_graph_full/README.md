# suture06_cluster_centric_prototype_graph_full

This independent-full candidate adapts Cluster-Centric Scanning Module into a scRNA prototype graph adapter.

The original reference scans image cluster centers instead of pixels. This implementation does not use images, Mamba, or pixel scanning. It maps cell latents to trainable prototype centers, propagates information through prototype-prototype cosine similarity, then injects a gated prototype context back into each cell latent.

## scMAE integration

The model keeps scMAE mask prediction and masked expression reconstruction as the primary training objective. The prototype branch is an auxiliary cluster-geometry stabilizer.

`scaled_expr` is used only as encoder input. `log_expr` is the reconstruction target. No count likelihood, NB/ZINB target, diffusion target, or cell mixing is used.

## Recommended defaults

- `n_prototypes = max(2 * n_clusters, 16)`
- `warmup_epochs = 10`
- `proto_weight = 0.15`
- `proto_temperature = 0.25`

Prototype centers are initialized from warmup base embeddings with KMeans, then trained as a memory bank.

## NeighborMix

Independent and complementary. This candidate does not mix cells, so `mixed_cell_fraction=0.0`.
