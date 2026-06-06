# Model Authenticity Checklist

This document provides detailed authenticity checklists for each migrated model, verifying that the original scCluBench model architecture, losses, and training procedure are preserved.

**Core principle**: A model is AUTHENTIC only if its core structure, losses, and training flow are preserved. Substitutions (KMeans-only, PCA-only, MSE-only, etc.) are strictly forbidden.

---

## DEC — Deep Embedded Clustering

### Checklist

| Component | Required by scCluBench | Preserved in `methods/DeepLearning/dec/run.py` | Notes |
|-----------|----------------------|----------------------------------------------|-------|
| Autoencoder class | Yes | **Yes** | Stacked encoder/decoder with ReLU activation |
| ClusteringLayer class | Yes | **Yes** | Student-t soft assignment with learnable cluster centers |
| DEC model class | Yes | **Yes** | Wraps AE + ClusteringLayer |
| `target_distribution()` function | Yes | **Yes** | KL target distribution: `weight.T / weight.sum(1)` |
| Pretraining with MSE loss | Yes | **Yes** | `model.pretrain()` with `F.mse_loss` |
| KMeans initialization of cluster centers | Yes | **Yes** | `KMeans(n_clusters, n_init=20)` on AE embedding |
| Joint KL divergence optimization | Yes | **Yes** | `F.kl_div(q.log(), p)` in clustering loop |
| Student-t q distribution | Yes | **Yes** | `1.0 / (1.0 + dist / alpha)` with power `(alpha+1)/2` |

### Forbidden Behaviors Check

| Forbidden | Evidence |
|-----------|----------|
| KMeans-only path | No; full DEC training always runs |
| PCA embedding shortcut | No; AE always trained |
| Deleted clustering layer | No; ClusteringLayer present and used |
| Replaced KL loss with MSE | No; KL divergence used in clustering phase |
| Skipped DEC training | No; both pretrain and clustering phases always execute |

### Verdict: **AUTHENTIC** ✓

---

## scDCC — Single-Cell Deep Constrained Clustering

### Checklist

| Component | Required by scCluBench | Preserved in `methods/DeepLearning/scDCC/` | Notes |
|-----------|----------------------|--------------------------------------------|-------|
| scDCC model class | Yes | **Yes** | In `scDCC.py` |
| Encoder network | Yes | **Yes** | `buildNetwork()` in `scDCC.py` |
| Decoder network | Yes | **Yes** | `buildNetwork()` in `scDCC.py` |
| ZINB decoder mean (`_dec_mean`) | Yes | **Yes** | `MeanAct()` activation on dec_h3 |
| ZINB decoder dispersion (`_dec_disp`) | Yes | **Yes** | `DispAct()` activation on dec_h3 |
| ZINB decoder pi (`_dec_pi`) | Yes | **Yes** | Sigmoid activation on dec_h3 |
| ZINBLoss class | Yes | **Yes** | In `layers.py`; ZINB negative binomial + zero-inflation |
| MeanAct activation | Yes | **Yes** | In `layers.py` |
| DispAct activation | Yes | **Yes** | In `layers.py` |
| Cluster parameter `mu` | Yes | **Yes** | `nn.Parameter` in `scDCC.py` |
| `soft_assign()` method | Yes | **Yes** | Student-t soft assignment |
| `target_distribution()` method | Yes | **Yes** | KL target distribution |
| Pretrain with ZINB loss | Yes | **Yes** | `pretrain_autoencoder()` uses ZINBLoss |
| Clustering fit with joint loss | Yes | **Yes** | `fit()` uses `cluster_loss + zinb_loss` |
| KMeans initialization | Yes | **Yes** | KMeans on latent representation |

### Forbidden Behaviors Check

| Forbidden | Evidence |
|-----------|----------|
| Ordinary AE + KMeans replacement | No; ZINB decoder and ZINBLoss present |
| ZINB loss replaced with MSE | No; ZINBLoss used in both pretrain and fit |
| Clustering layer removed | No; `mu` parameter and `soft_assign()` present |
| Skipped constrained/clustering training | No; `fit()` always runs clustering phase |
| Simplified autoencoder without ZINB | No; ZINB decoder has `_dec_mean`, `_dec_disp`, `_dec_pi` |

### Verdict: **AUTHENTIC** ✓

---

## scDSC (SDCN) — Structural Deep Clustering Network

### Checklist

| Component | Required by scCluBench | Preserved in `methods/GNN/scDSC/` | Notes |
|-----------|----------------------|-----------------------------------|-------|
| AE class | Yes | **Yes** | 3-layer encoder/decoder with BatchNorm |
| GNN layers (`gnn_1` to `gnn_5`) | Yes | **Yes** | In `SDCN.__init__` from `GNN.py` |
| Graph adjacency input | Yes | **Yes** | KNN graph built from data |
| Cluster layer | Yes | **Yes** | `nn.Parameter` in `SDCN` |
| ZINB decoder mean (`_dec_mean`) | Yes | **Yes** | Sequential: Linear → MeanAct |
| ZINB decoder dispersion (`_dec_disp`) | Yes | **Yes** | Sequential: Linear → DispAct |
| ZINB decoder pi (`_dec_pi`) | Yes | **Yes** | Sequential: Linear → Sigmoid |
| ZINBLoss | Yes | **Yes** | In `layers.py` |
| `target_distribution()` | Yes | **Yes** | In `run.py` |
| Joint AE/GNN/clustering training | Yes | **Yes** | Combined loss: BCE + CE + MSE + ZINB |

### Forbidden Behaviors Check

| Forbidden | Evidence |
|-----------|----------|
| AE-only path | No; GNN layers always in forward pass |
| Graph clustering only | No; AE branch always in forward pass |
| Leiden/Louvain substitution | No; no Leidenalg or louvain_communities in SDCN |
| ZINB decoder removed | No; `_dec_mean`, `_dec_disp`, `_dec_pi` present |
| Simplified replacement model | No; full SDCN architecture preserved |

### Verdict: **AUTHENTIC** ✓

---

## ScanpyStandard

### Checklist

| Component | Required by scCluBench | Preserved in `methods/Traditional/ScanpyStandard/run.py` | Notes |
|-----------|----------------------|--------------------------------------------------------|-------|
| Scanpy QC (filter_cells/genes) | Yes | **Yes** | Steps 2, QC |
| normalize_total | Yes | **Yes** | Step 3 |
| log1p transformation | Yes | **Yes** | Step 4 |
| HVG selection (seurat) | Yes | **Yes** | Step 5 |
| PCA dimensionality reduction | Yes | **Yes** | Step 7 |
| KNN neighbors graph | Yes | **Yes** | Step 8 |
| UMAP | Yes | **Yes** | Step 9 |
| Leiden clustering | Yes | **Yes** | Step 10 |
| Resolution auto-search | Yes | **Yes** | `res_search_fixed_clus()` |

### Verdict: **AUTHENTIC** ✓

---

## Leiden

### Checklist

| Component | Required by scCluBench | Preserved in `methods/Traditional/Leiden/run.py` | Notes |
|-----------|----------------------|-------------------------------------------------|-------|
| Leiden algorithm | Yes | **Yes** | `leidenalg.find_partition()` |
| RBConfigurationVertexPartition | Yes | **Yes** | Supports resolution parameter |
| igraph graph input | Yes | **Yes** | Converted from scanpy neighbors |
| Resolution tuning | Yes | **Yes** | Auto-tune via NMI search |
| KNN graph construction | Yes | **Yes** | From scanpy neighbors |

### Verdict: **AUTHENTIC** ✓

---

## Louvain

### Checklist

| Component | Required by scCluBench | Preserved in `methods/Traditional/Louvain/run.py` | Notes |
|-----------|----------------------|------------------------------------------------|-------|
| Louvain algorithm | Yes | **Yes** | `networkx.algorithms.community.louvain_communities` |
| Resolution parameter | Yes | **Yes** | Passed to `louvain_communities` |
| NetworkX graph input | Yes | **Yes** | Converted from scanpy neighbors |
| Resolution tuning | Yes | **Yes** | Auto-tune via NMI search |
| KNN graph construction | Yes | **Yes** | From scanpy neighbors |

### Verdict: **AUTHENTIC** ✓

---

## sc3 — Single-Cell Consensus Clustering

### Checklist

| Component | Required by scCluBench | Preserved in `methods/Traditional/sc3/run.py` | Notes |
|-----------|----------------------|--------------------------------------------|-------|
| PCA dimensionality reduction | Yes | **Yes** | `PCA(n_pcs)` |
| Consensus matrix construction | Yes | **Yes** | Multiple K-means + hierarchical clustering |
| K-means ensemble | Yes | **Yes** | Multiple K-means configurations |
| Hierarchical clustering on consensus | Yes | **Yes** | `AgglomerativeClustering` on distance matrix |
| Final consensus clustering | Yes | **Yes** | Hierarchical clustering result |

### Verdict: **AUTHENTIC** ✓

### Note
Pure Python reimplementation (original SC3 requires R). Algorithm principle preserved: consensus clustering via multiple clustering runs and hierarchical aggregation. No KMeans-only shortcut — the consensus matrix + hierarchical clustering pipeline is fully implemented.

---

## scDeepCluster

### Checklist

| Component | Required by scCluBench | Preserved in `methods/DeepLearning/scDeepCluster/` | Notes |
|-----------|----------------------|---------------------------------------------------|-------|
| SCDeepCluster model class | Yes | **Yes** (in `code/scDeepCluster.py`) | Not independently verified |
| `pretrain()` method | Yes | **Yes** (in `code/scDeepCluster.py`) | Not independently verified |
| `fit()` method | Yes | **Yes** (in `code/scDeepCluster.py`) | Not independently verified |
| `extract_feature()` method | Yes | **Yes** (in `code/scDeepCluster.py`) | Not independently verified |
| TensorFlow dependency | Yes | **Yes** | ENV-GATED |

### Verdict: **ENV-GATED** (not independently verified — TensorFlow required)

---

## scGNN

### Checklist

| Component | Required by scCluBench | Preserved in `methods/GNN/scGNN/` | Notes |
|-----------|----------------------|----------------------------------|-------|
| AE/VAE model | Yes | **PENDING** | Code migrated, not independently verified |
| GAE embedding | Yes | **PENDING** | `gae_embedding.py` present |
| EM clustering | Yes | **PENDING** | `scGNN.py` contains EM logic |
| Graph construction | Yes | **PENDING** | `graph_function.py` present |
| LTMG regularization | Yes | **PENDING** | `LTMG_R.py` present |
| `--gpu` default 0 | Forbidden | **FAIL** | Violates BDD Scenario 13 |

### Verdict: **PENDING** (GPU policy violation)

---

## scCDCG

### Checklist

| Component | Required by scCluBench | Preserved in `methods/GNN/scCDCG/` | Notes |
|-----------|----------------------|-----------------------------------|-------|
| Graph construction | Yes | **PENDING** | `scCDCG_preprocess.py` present |
| Contrastive module | Yes | **PENDING** | Model structure not yet verified |
| Clustering module | Yes | **PENDING** | Model structure not yet verified |
| `--gpu` default 0 | Forbidden | **FAIL** | Violates BDD Scenario 13 |

### Verdict: **PENDING** (GPU policy violation)

---

## AttentionAE_sc

### Checklist

| Component | Required by scCluBench | Preserved in `methods/GNN/AttentionAE_sc/` | Notes |
|-----------|----------------------|-------------------------------------------|-------|
| Attention autoencoder | Yes | **PENDING** | `model.py` present, not independently verified |
| Original loss | Yes | **PENDING** | `loss.py` present, not independently verified |
| Training loop | Yes | **PENDING** | `train.py` present, not independently verified |
| `--gpu` default 0 | Forbidden | **FAIL** | Violates BDD Scenario 13 |

### Verdict: **PENDING** (GPU policy violation)

---

## Audit Script

Model authenticity is also verified programmatically by `scripts/audit_model_authenticity.py`. See that script for automated keyword and structural checks.
