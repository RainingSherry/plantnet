# PlantDiffCluster

## GRN-Support Graph Diffusion for Plant Single-Cell RNA-Seq Clustering

**PlantDiffCluster** is an interpretable deep learning model for dimensionality reduction and clustering of plant single-cell RNA-seq data. It combines ideas from DOLORIS (GRN-conditional networks + diffusion) and PhytoCluster (VAE + GMM for plant scRNA-seq) with a novel **Gene Support Graph** architecture.

### Core Architecture

```
Input: scRNA-seq gene expression X (n_cells, n_genes)
   │
   ├─→ HVG Selection + Gene Graph Construction
   │     ├─ Co-expression Graph (Pearson correlation on HVG)
   │     ├─ Marker Graph (Plant marker-boosted)
   │     └─ Random Graph (negative control)
   │
   ├─→ GeneGATEncoder (GAT on gene graph)
   │     h_g = GAT(gene_embedding, A)
   │
   ├─→ SupportPooling (gene → cell)
   │     z_c = Pool_{g ∈ S_c}( w_{cg} · h_g )
   │     S_c = {g | X_{cg} > 0}  ← cell's expression support set
   │
   ├─→ MaskDiffusionRefiner (DDPM/DDIM denoising)
   │     z_0^raw → [add noise] → z_t ──→ denoise → z_0^refined
   │
   └─→ ClusterHead (GMM / DEC / Contrastive)
         q(k|z) = softmax(...), L_cluster

Output: cell embedding z (n_cells, D)
        cluster assignment (n_cells,)
        sparsity mask M (n_cells, n_hvg)
```

### Key Innovation

| Traditional VAE | PlantDiffCluster |
|---|---|
| All 20,000 genes go through MLP | Only expressed genes (support set) go through GAT |
| Dense MLP encoding | Sparse graph-based message passing |
| Black-box latent dimensions | Interpretable: each dimension traces to specific genes/GRN modules |

**Core expression:**
```
h_g = GAT_GRN(e_g)                    ← gene contextual embedding
z_c = Pool_{g ∈ S_c}( w_{cg} · h_g )  ← cell embedding from expressed genes
```

Where `S_c = {g | X_{cg} > 0}` is the **expression support set** — the genes actually expressed in cell c.

### Installation

```bash
pip install -r requirements.txt
```

### Quick Start

```bash
# Training (Co-expression graph, 15 clusters, 100 epochs)
python train/train_plantdiffcluster.py \
    --data_path ../../../data/SRP182008.h5ad \
    --graph_type coexpression \
    --n_clusters 15 \
    --epochs 100 \
    --save_dir ./results/srp182008_coexp

# Evaluation
python eval/evaluate.py \
    --model_path ./results/srp182008_coexp/best_model.pt \
    --data_path ../../../data/SRP182008.h5ad \
    --output_dir ./eval_results
```

### Ablation Experiments

Three gene graph types for ablation:

```bash
# 1. HVG Co-expression Graph (recommended)
--graph_type coexpression

# 2. Marker-boosted Graph
--graph_type marker

# 3. Random Graph (negative control)
--graph_type random
```

### Support Set Weighting Strategies

Four strategies for `w_{cg}` in the support pooling:

| Strategy | Formula | When to Use |
|---|---|---|
| `log1p` (default) | `log(1 + X_{cg})` | General purpose, recommended start |
| `rank` | Rank-normalized weights | When expression values are unreliable |
| `tfidf` | TF × IDF | Highlight rare but important genes |
| `norm` | L2-normalized count | When total counts vary wildly |

### Supported Datasets

- **SRP182008**: 13,514 Arabidopsis root tip cells, 15 cell types (PlantNet data)

### Model Components

| Module | File | Description |
|---|---|---|
| Gene Graph Builder | `graphs/build_gene_graph.py` | Co-expression / Marker / Random graphs |
| Support Graph Builder | `graphs/build_cell_gene_graph.py` | Cell-gene bipartite support graph |
| Gene GAT Encoder | `models/gene_gat_encoder.py` | GAT on gene regulatory graph |
| Support Pooling | `models/support_pooling.py` | Attention-weighted gene aggregation |
| Mask Diffusion Refiner | `models/mask_diffusion_refiner.py` | DDPM/DDIM denoising of cell embeddings |
| Cluster Head | `models/cluster_head.py` | GMM / DEC / Contrastive clustering |
| Main Model | `models/plantdiffcluster.py` | End-to-end PlantDiffCluster |

### Citation

If this work is helpful for your research, please cite:

```bibtex
# (Placeholder - update with your paper when available)
```
