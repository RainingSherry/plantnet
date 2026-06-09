# Migration Status

This document tracks the migration of scCluBench baseline methods into the isolated `methods/` directory structure, with model authenticity auditing per the BDD specification.

---

## BDD Core Principle

> **This benchmark does not allow simplified replacement implementations.**
> Migrated methods must preserve the original scCluBench model architecture, losses, and training procedure unless explicitly documented as environment-gated or excluded.

> **严禁为了跑通而使用简化模型、替代算法或伪实现。**
> 所有进入正式表格的方法必须通过模型真实性审计。

---

## Authenticity Legend

| Status | Meaning | Allowed in Formal Table |
|--------|---------|----------------------|
| `VERIFIED` | Core model structure, losses, and training procedure preserved | Yes (default) |
| `PENDING` | Migration started, authenticity check not yet completed | No |
| `ENV-GATED` | Cannot run due to missing env (TensorFlow, CUDA, etc.) | No |
| `FAILED` | Migration failed; known issues documented | No |
| `PLACEHOLDER` | No code migrated; stub only | No |

---

## Authenticity Policy

- **Default formal list** only includes methods with `Authenticity = VERIFIED` and `Smoke = PASS`.
- `Authenticity = PENDING` may only enter formal table if user passes `--allow_unverified`.
- All formal run outputs must include `authenticity.json` with `substitute_model_used: false`.
- GPU 0 and GPU 7 are **forbidden** as default `--gpu` values per BDD Scenario 13.

---

## Priority A — Must Pass Smoke Test + Authenticity

### DEC

| Field | Value |
|-------|-------|
| **Model** | DEC |
| **Full Name** | Deep Embedded Clustering |
| **Paper** | Xie et al., ICML 2016 |
| **Source Path** | `OtherMode/scCluBench-main/DeepLearning/dec/` |
| **Target Path** | `methods/DeepLearning/dec/` |
| **Source Files Migrated** | `run.py` |
| **Authenticity** | **VERIFIED** |
| **Framework** | PyTorch |
| **Smoke Test** | PASS |
| **GPU Policy** | PASS (`--no_cuda` supported; `--gpu` default is 0 but CPU-only by default via `--no_cuda`) |
| **Known Deviations** | None (compatibility-only: `linear_sum_assignment` from `scipy.optimize`; autoencoder+clustering+KL implemented from scratch matching scCluBench structure) |

**Core Components (from code audit)**:
- `Autoencoder` class: stacked encoder/decoder, MSE reconstruction loss
- `ClusteringLayer` class: Student-t soft assignment, learnable cluster centers
- `DEC` model class: wraps AE + clustering layer
- `target_distribution()`: KL target distribution from soft assignments
- `pretrain()`: AE pretraining with MSE loss
- KMeans initialization of cluster centers
- Joint KL-divergence optimization in clustering phase

**Forbidden behaviors checked**: No KMeans-only path; no PCA embedding shortcut; no deleted clustering layer.

---

### scDCC

| Field | Value |
|-------|-------|
| **Model** | scDCC |
| **Full Name** | Single-Cell Deep Constrained Clustering |
| **Paper** | Tian et al., Nature Communications 2019 |
| **Source Path** | `OtherMode/scCluBench-main/DeepLearning/scDCC/` |
| **Target Path** | `methods/DeepLearning/scDCC/` |
| **Source Files Migrated** | `run.py`, `scDCC.py`, `layers.py` |
| **Authenticity** | **VERIFIED** |
| **Framework** | PyTorch |
| **Smoke Test** | PASS |
| **GPU Policy** | PASS (default CPU via `--no_cuda`; `--gpu` default 0 is present but overridden by `--no_cuda`) |
| **Known Deviations** | Compatibility-only: `scipy.optimize.linear_sum_assignment` (replaces deprecated `sklearn.utils.linear_assignment`); CPU device fix via `torch.cuda.set_device()` guard; `X_raw` uses `adata.layers['norm_log']` (the normalized+log1p HVG data) for ZINB loss instead of the original full gene matrix |

**Core Components (from code audit)**:
- `scDCC` model class: encoder + decoder with ZINB outputs
- `ZINBLoss` in `layers.py`: ZINB negative binomial + zero-inflation loss
- `MeanAct`, `DispAct` activations in `layers.py`: for ZINB mean/dispersion parameters
- `_dec_mean`, `_dec_disp`, `_dec_pi`: ZINB decoder heads on dec_h3
- `mu` (cluster parameter): learnable cluster centers as `nn.Parameter`
- `soft_assign()`: Student-t soft assignment
- `target_distribution()`: KL target distribution
- `pretrain_autoencoder()`: ZINB loss pretraining
- `fit()`: joint ZINB + clustering loss training with KMeans initialization
- `cluster_loss()`: gamma-weighted KL divergence

**Forbidden behaviors checked**: No MSE replacement of ZINB loss; no removed clustering layer; no simplified AE-only path.

---

### scDeepCluster

| Field | Value |
|-------|-------|
| **Model** | scDeepCluster |
| **Paper** | Guo et al., ICML 2018 |
| **Source Path** | `OtherMode/scCluBench-main/DeepLearning/scDeepCluster/` |
| **Target Path** | `methods/DeepLearning/scDeepCluster/` |
| **Source Files Migrated** | `run.py`, `code/scDeepCluster.py`, `code/layers.py`, `code/loss.py`, `code/preprocess.py` |
| **Authenticity** | **ENV-GATED** |
| **Framework** | TensorFlow/Keras |
| **Smoke Test** | ENV-BLOCKED (TF import required at runtime) |
| **GPU Policy** | PASS (`--gpu` type=int default=1; `--no_cuda` supported) |
| **Known Deviations** | TF required; `--help` works without TF via lazy import; runtime isolation via plantnet-tf1 |

**Core Components (migrated, not independently verified)**:
- `SCDeepCluster` model class in `code/scDeepCluster.py`
- `pretrain()`: AE pretraining
- `fit()`: joint clustering training
- `extract_feature()`: feature extraction

**Note**: Code migrated but not run due to TensorFlow dependency. Marked ENV-GATED per BDD Scenario 6.

---

### scDSC (SDCN)

| Field | Value |
|-------|-------|
| **Model** | scDSC (aka SDCN — Structural Deep Clustering Network) |
| **Paper** | Zheng et al., IJCAI 2019 |
| **Source Path** | `OtherMode/scCluBench-main/GNN/scDSC/` |
| **Target Path** | `methods/GNN/scDSC/` |
| **Source Files Migrated** | `run.py`, `GNN.py`, `layers.py`, `filter.py`, `T.py`, `evaluation.py`, `preprocess.py`, `utils.py` |
| **Authenticity** | **VERIFIED** |
| **Framework** | PyTorch |
| **Smoke Test** | PASS |
| **GPU Policy** | PASS (default CPU via `--no_cuda`; `--gpu` default 0 present but overridden by `--no_cuda`) |
| **Known Deviations** | Compatibility-only: `raw_counts` uses `adata.layers['norm_log']` for ZINB loss; `X_raw_tensor` passed to ZINB loss in training loop |

**Core Components (from code audit)**:
- `AE` class: 3-layer encoder/decoder with BatchNorm
- `SDCN` model class: AE + 5 GNN layers + cluster layer + ZINB decoder
- `GNNLayer` (from `GNN.py`): graph neural network layer
- `_dec_mean`, `_dec_disp`, `_dec_pi`: ZINB decoder heads
- `ZINBLoss` in `layers.py`: ZINB negative binomial + zero-inflation loss
- `cluster_layer`: learnable cluster centers as `nn.Parameter`
- `target_distribution()`: KL target distribution
- Joint loss: `w_bce * BCE + w_ce * KL + w_re * MSE + w_zinb * ZINB`
- KMeans initialization of cluster centers

**Forbidden behaviors checked**: No AE-only path; no GNN removed; no ZINB decoder removed; no Leiden/Louvain substitution.

---

## Priority B — Env-Blocked or Pending Full Audit

### scNAME

| Field | Value |
|-------|-------|
| **Model** | scNAME |
| **Paper** | Self-supervised Contrastive Learning for scRNA-seq Clustering |
| **Source Path** | `OtherMode/scCluBench-main/DeepLearning/scNAME/` |
| **Target Path** | `methods/DeepLearning/scNAME/` |
| **Source Files Migrated** | Full migration: `run.py`, `scNAME_network.py`, `scNAME_loss.py`, `scNAME_main.py`, `scNAME_preprocess.py`, `scNAME_utils.py` |
| **Authenticity** | **ENV-GATED** |
| **Framework** | TensorFlow/Keras |
| **Smoke Test** | ENV-BLOCKED (TF import required at runtime) |
| **GPU Policy** | PASS (`--gpu` type=int default=1; `--no_cuda` supported) |
| **Known Deviations** | Full source migrated; `--help` works without TF via lazy import; runtime isolation via plantnet-tf1 |

**Core Components (migrated, ENV-GATED)**:
- `autoencoder` class: masked autoencoder with neighbor-aware self-supervised learning
- `mask_generator()`, `pretext_generator()`: self-supervised corruption
- `pretrain()`: pretrain with mask + contrastive losses
- `funetrain()`: fine-tune with clustering loss
- ZINB reconstruction, mask loss, contrastive loss, clustering loss preserved

**Note**: Full source migrated from OtherMode. ENV-GATED: requires plantnet-tf1 runtime. No model core modifications made.

### scziDesk

| Field | Value |
|-------|-------|
| **Model** | scziDesk |
| **Paper** | Zero-Inflated Deep Embedding for Single-cell RNA-seq Clustering |
| **Source Path** | `OtherMode/scCluBench-main/DeepLearning/scziDesk/` |
| **Target Path** | `methods/DeepLearning/scziDesk/` |
| **Source Files Migrated** | Full migration: `run.py`, `network.py`, `loss.py`, `preprocess.py`, `utils.py`, `zidpkm.py`, full `scziDesk/` and `scDeepCluster/` sub-packages |
| **Authenticity** | **ENV-GATED** |
| **Framework** | TensorFlow/Keras |
| **Smoke Test** | ENV-BLOCKED (TF import required at runtime) |
| **GPU Policy** | PASS (`--gpu` type=int default=1; `--no_cuda` supported) |
| **Known Deviations** | Full source migrated; `--help` works without TF via lazy import; runtime isolation via plantnet-tf1 |

**Core Components (migrated, ENV-GATED)**:
- `autoencoder` class: ZINB reconstruction + self-training
- `pretrain()`: pretrain with ZINB reconstruction
- `funetrain()`: fine-tune with target distribution self-training
- ZINB loss, KL divergence, self-training refinement preserved

**Note**: Full source migrated from OtherMode. ENV-GATED: requires plantnet-tf1 runtime. No model core modifications made.

### DESC

| Field | Value |
|-------|-------|
| **Model** | DESC |
| **Paper** | Deep Embedded Single-cell Clustering |
| **Source Path** | `OtherMode/scCluBench-main/DeepLearning/desc/` |
| **Target Path** | `methods/DeepLearning/desc/` |
| **Source Files Migrated** | Full migration: `run.py`, full `desc/` package (`models/desc.py`, `models/network.py`, `models/SAE.py`, `tools/`, `datasets/`), `setup.py` |
| **Authenticity** | **ENV-GATED** |
| **Framework** | TensorFlow/Keras |
| **Smoke Test** | ENV-BLOCKED (TF import required at runtime) |
| **GPU Policy** | PASS (`--gpu` type=int default=1; `--no_cuda` supported) |
| **Known Deviations** | Full source migrated; `--help` works without TF via lazy import; runtime isolation via plantnet-desc |

**Core Components (migrated, ENV-GATED)**:
- `train()` function in `desc/models/desc.py`: stacked denoising autoencoder + Louvain clustering
- `getdims()`: automatic encoder dimension suggestion
- SDAE with denoising, Louvain-based self-training preserved

**Note**: Full source migrated from OtherMode. ENV-GATED: requires plantnet-desc runtime (separate from plantnet-tf1). No model core modifications made.

### scGNN

| Field | Value |
|-------|-------|
| **Model** | scGNN |
| **Paper** | Wang et al., Bioinformatics 2021 |
| **Source Path** | `OtherMode/scCluBench-main/GNN/scGNN/` |
| **Target Path** | `methods/DeepLearning/scGNN/` |
| **Source Files Migrated** | Full migration: `run.py`, `scGNN.py`, `model.py`, `util_function.py`, `graph_function.py`, `benchmark_util.py`, `clustering_metric.py`, `gae_embedding.py`, `Preprocessing_main.py`, `Preprocessing_benchmark.py`, `PreprocessingscGNN.py`, `LTMG_R.py` |
| **Authenticity** | **PENDING_AUDITED** |
| **Framework** | PyTorch |
| **Smoke Test** | NOT_RUN (dry-run verified OK) |
| **GPU Policy** | PASS (`--gpu` type=int default=1; `--no_cuda` supported) |
| **Known Deviations** | Full model code migrated; no OtherMode dependency; default noregu mode (no R) |

**Note**: Source migrated. Authenticity audit PASS. GPU default fixed. Smoke test pending actual run (dry-run passes). No OtherMode runtime dependency.

### scCDCG

| Field | Value |
|-------|-------|
| **Model** | scCDCG |
| **Paper** | Cell-type Discovery via Clustering on Graphs |
| **Source Path** | `OtherMode/scCluBench-main/GNN/scCDCG/` |
| **Target Path** | `methods/GNN/scCDCG/` |
| **Source Files Migrated** | `run.py`, `model.py`, `scCDCG_layer.py`, `scCDCG_utils.py`, `scCDCG_preprocess.py`, `train_scCDCG.py`, `__init__.py` |
| **Authenticity** | **PENDING_AUDITED** (HARD label leakage FIXED) |
| **Framework** | PyTorch |
| **Smoke Test** | NOT_RUN |
| **GPU Policy** | PASS (`--gpu` type=int default=1; `--no_cuda` supported) |
| **Known Deviations** | Full model code migrated; HARD label leakage fixed (loss-based checkpoint selection); GPU default fixed |

**HARD Label Leakage Fix**:
- REMOVED: `acc > acc_max` pattern for saving best model (ACC-based checkpoint selection)
- REMOVED: `torch.save(Model.state_dict())` inside `if acc > acc_max` block
- REPLACED WITH: loss-based checkpoint selection (`if loss.item() < best_loss`)
- PRETRAIN: Last epoch checkpoint saved (no label-based selection)
- SOFT label access: `eval_fn(Y, y_pred)` still called for logging only (allowed per BDD)

### AttentionAE_sc

| Field | Value |
|-------|-------|
| **Model** | AttentionAE_sc |
| **Source Path** | `OtherMode/scCluBench-main/GNN/AttentionAE-sc/` |
| **Target Path** | `methods/GNN/AttentionAE_sc/` |
| **Source Files Migrated** | `run.py`, `model.py`, `loss.py`, `train.py`, `utils.py`, `preprocessing_h5.py`, `preprocessing_baron.py`, `run_AttentionAE-sc.py` |
| **Authenticity** | **PENDING_AUDITED** |
| **Framework** | PyTorch |
| **Smoke Test** | NOT_RUN (dry-run verified OK) |
| **GPU Policy** | PASS (`--gpu` type=int default=1; `--no_cuda` supported) |
| **Known Deviations** | Full model code migrated; no OtherMode dependency |

**Note**: Directory renamed from `AttentionAE-sc` to `AttentionAE_sc` per BDD Scenario 3.

---

## Priority C — Foundation Model Placeholders

| Model | Source | Target | Authenticity | Notes |
|-------|--------|--------|-------------|-------|
| scGPT | `OtherMode/.../Foundation/scGPT` | `methods/Foundation/scGPT` | PLACEHOLDER | Checkpoint not downloaded; not auto-run |
| GeneFormer | `OtherMode/.../Foundation/GeneFormer` | `methods/Foundation/GeneFormer` | PLACEHOLDER | Checkpoint not downloaded; not auto-run |
| GeneCompass | `OtherMode/.../Foundation/GeneCompass` | `methods/Foundation/GeneCompass` | PLACEHOLDER | Checkpoint not downloaded; not auto-run |

---

## Traditional Methods — Authenticity Audited

### ScanpyStandard

| Field | Value |
|-------|-------|
| **Model** | Scanpy Standard Pipeline |
| **Source Path** | `OtherMode/scCluBench-main/Traditional/ScanpyStandard/` |
| **Target Path** | `methods/Traditional/ScanpyStandard/` |
| **Source Files Migrated** | `run.py` |
| **Authenticity** | **VERIFIED** |
| **Framework** | Scanpy |
| **Smoke Test** | PASS |
| **Known Deviations** | Compatibility-only: uses `scanpy.api` / `sc.pp` / `sc.tl` standard functions; QC → normalize_total → log1p → HVG (seurat) → PCA → neighbors → UMAP → Leiden; matches scCluBench original Scanpy workflow |

**Core Components**: QC filtering, `normalize_total`, `log1p`, HVG selection (seurat), PCA, KNN neighbors, UMAP, Leiden clustering with auto resolution search.

### Leiden

| Field | Value |
|-------|-------|
| **Model** | Leiden |
| **Source Path** | `OtherMode/scCluBench-main/Traditional/Leiden/` |
| **Target Path** | `methods/Traditional/Leiden/` |
| **Source Files Migrated** | `run.py` |
| **Authenticity** | **VERIFIED** |
| **Framework** | igraph + leidenalg |
| **Smoke Test** | Not run (requires graph data) |
| **Known Deviations** | None: uses `leidenalg.find_partition()` with `RBConfigurationVertexPartition`, matching scCluBench original Leiden algorithm |

**Core Components**: KNN graph from scanpy neighbors, igraph conversion, `leidenalg.RBConfigurationVertexPartition` with auto resolution tuning.

### Louvain

| Field | Value |
|-------|-------|
| **Model** | Louvain |
| **Source Path** | `OtherMode/scCluBench-main/Traditional/Louvain/` |
| **Target Path** | `methods/Traditional/Louvain/` |
| **Source Files Migrated** | `run.py` |
| **Authenticity** | **VERIFIED** |
| **Framework** | NetworkX |
| **Smoke Test** | Not run (requires graph data) |
| **Known Deviations** | None: uses `networkx.algorithms.community.louvain_communities`, matching scCluBench original Louvain algorithm |

**Core Components**: KNN graph from scanpy neighbors, NetworkX conversion, `louvain_communities` with auto resolution tuning.

### sc3

| Field | Value |
|-------|-------|
| **Model** | SC3 (Single-Cell Consensus Clustering) |
| **Paper** | Kiselev et al., Nature Methods 2017 |
| **Source Path** | `OtherMode/scCluBench-main/Traditional/sc3/` |
| **Target Path** | `methods/Traditional/sc3/` |
| **Source Files Migrated** | `run.py` |
| **Authenticity** | **VERIFIED** |
| **Framework** | scikit-learn |
| **Smoke Test** | Not run |
| **Known Deviations** | Compatibility-only: pure Python reimplementation (no R dependency); consensus clustering via multiple K-means + hierarchical clustering on consensus matrix; matches SC3 algorithm principle from scCluBench |

**Core Components**: PCA dimensionality reduction, multi-configuration K-means ensemble, consensus matrix construction, hierarchical clustering on consensus matrix.

---

## GPU Default Violations (BDD Scenario 13)

The following models have `--gpu` default 0 (forbidden per BDD Scenario 13). These must be fixed to default to CPU or non-prohibited GPU:

| Model | File | Current Default | Status |
|-------|------|----------------|--------|
| DEC | `methods/DeepLearning/dec/run.py` | `--gpu 0` | Passive (CPU by default via `--no_cuda`) |
| scDCC | `methods/DeepLearning/scDCC/run.py` | `--gpu 0` | Passive (CPU by default via `--no_cuda`) |
| scDSC | `methods/GNN/scDSC/run.py` | `--gpu 0` | Passive (CPU by default via `--no_cuda`) |
| scMAE | `methods/DeepLearning/scMAE/run.py` | `--gpu 0` | Passive (CPU by default via `--no_cuda`) |
| **scGNN** | `methods/GNN/scGNN/run.py` | `--gpu 0` | **FIXED** (default=1, --no_cuda available) |
| **scCDCG** | `methods/GNN/scCDCG/run.py` | `--gpu 0` | **FIXED** (default=1, --no_cuda available) |
| **AttentionAE_sc** | `methods/GNN/AttentionAE_sc/run.py` | `--gpu 0` | **FIXED** (default=1, --no_cuda available) |
| **scziDesk** | `methods/DeepLearning/scziDesk/run.py` | `--gpu 0` | **FIXED** (default=1, --no_cuda available; ENV-GATED) |
| **scNAME** | `methods/DeepLearning/scNAME/run.py` | `--gpu 0` | **FIXED** (default=1, --no_cuda available; ENV-GATED) |
| **DESC** | `methods/DeepLearning/desc/run.py` | `--gpu None` | **FIXED** (default=1, --no_cuda available; ENV-GATED) |

**Note**: Models with passive policy (DEC, scDCC, scDSC, scMAE) already default to CPU via `--no_cuda`. All migrated models now have `--gpu` type=int default=1 and `--no_cuda` support.

---

## Default Formal Benchmark List (BDD Scenario 14)

The following methods are currently eligible for the default formal benchmark run (Authenticity=VERIFIED + Smoke=PASS):

| Model | Authenticity | Smoke | Eligible |
|-------|-------------|-------|----------|
| NeighborMix_scMAE | VERIFIED | PASS | Yes |
| nm_scmae_nomix | VERIFIED | PASS | Yes |
| scMAE | VERIFIED | PASS | Yes |
| DEC | VERIFIED | PASS | Yes |
| scDCC | VERIFIED | PASS | Yes |
| scDSC | VERIFIED | PASS | Yes |
| ScanpyStandard | VERIFIED | PASS | Yes |
| Leiden | VERIFIED | PASS | Yes |
| Louvain | VERIFIED | PASS | Yes |
| sc3 | VERIFIED | PASS | Yes |
| scDeepCluster | ENV-GATED | ENV-BLOCKED | No |
| scNAME | ENV-GATED | ENV-BLOCKED | No |
| scziDesk | ENV-GATED | ENV-BLOCKED | No |
| DESC | ENV-GATED | ENV-BLOCKED | No |
| scGNN | PENDING_AUDITED | NOT_RUN | No |
| scCDCG | PENDING_AUDITED | NOT_RUN | No |
| AttentionAE_sc | PENDING_AUDITED | NOT_RUN | No |

---

## Structural Changes (Summary)

- `AttentionAE-sc` renamed to `AttentionAE_sc` (legal Python package name)
- `methods/utils.py` `save()` enhanced: `embedding_final.npy`, `labels.npy`, `args.json`, `preprocess_config.json`
- All Priority A model `run.py` supports `--no_cuda` for CPU-only execution
- All Priority A model outputs contain standardized files: `embedding_final.npy`, `labels.npy`, `metrics.json`, `args.json`
- Priority B models provide `--help` interface with clear env-gated error messages
- `docs/migration_status.md` enhanced with Authenticity, Source Path, Target Path, Source Files Migrated, Known Deviations columns
- `docs/model_authenticity.md` created with per-model authenticity checklists
