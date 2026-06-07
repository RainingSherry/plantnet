# PlantNet Single-Cell Benchmark

This repository contains a multi-baseline benchmark for single-cell clustering on plant datasets. The benchmark provides an isolated `methods/` directory with standardized data interfaces, evaluation metrics, and unified result formats.

## Project Structure

```
methods/
├── shared_utils.py          # Shared utilities
├── preprocess.py            # Standard data preprocessing (h5ad → normalized/scaled)
├── evaluation.py           # 8 clustering metrics (ACC, NMI, ARI, F1-macro, FMI, V-measure, Homogeneity, Completeness)
├── utils.py               # save() interface and helpers
├── DeepLearning/
│   ├── NeighborMix_scMAE/  # [PROPOSED] NeighborMix-scMAE
│   ├── NeighborMix_MAE/    # NeighborMix-MAE (MAE variant)
│   ├── scMAE/             # scMAE baseline
│   ├── dec/               # DEC (Deep Embedded Clustering)
│   ├── scDCC/             # scDCC (Deep Constrained Clustering)
│   ├── scDeepCluster/     # scDeepCluster [TensorFlow, env-blocked if TF not installed]
│   ├── scNAME/             # scNAME [env-blocked: requires TensorFlow]
│   ├── scziDesk/           # scziDesk [env-blocked: requires TensorFlow]
│   ├── desc/               # DESC [env-blocked: requires TensorFlow]
│   ├── scVI/              # scVI baseline
│   ├── PhytoCluster/       # PhytoCluster
│   └── PlantSPADE_LGCL/    # [LEGACY] Retained for reference
├── GNN/
│   ├── scDSC/             # scDSC (Structural Deep Clustering Network)
│   ├── scGNN/             # scGNN [env-blocked]
│   ├── scCDCG/            # scCDCG [env-blocked]
│   └── AttentionAE_sc/    # AttentionAE [env-blocked]
├── Foundation/
│   ├── scGPT/             # scGPT [placeholder: requires checkpoint download]
│   ├── GeneFormer/        # GeneFormer [placeholder]
│   ├── GeneCompass/       # GeneCompass [placeholder]
│   └── scPlantLLM/        # scPlantLLM [placeholder]
└── Traditional/
    ├── ScanpyStandard/    # Scanpy standard pipeline
    ├── Leiden/            # Leiden clustering
    ├── Louvain/           # Louvain clustering
    └── sc3/              # SC3 clustering

OtherMode/
└── scCluBench-main/      # [REFERENCE ONLY] Original scCluBench repository.
                           # Kept as-is for historical reference. Not imported by runnable methods.
```

## Method Roles

| Role | Methods | Notes |
|---|---|---|
| **Proposed** | NeighborMix-scMAE | Main method in this benchmark |
| **Ablation** | NM-scMAE-noMix | NeighborMix-scMAE without neighbor mixing |
| **First baselines** | DEC, scDCC, scDeepCluster, scDSC, scMAE | Priority A PyTorch baselines, smoke-tested |
| **TF baselines** | scNAME, scziDesk, DESC, scDeepCluster | TensorFlow/Keras; env-blocked if TF not installed |
| **GNN baselines** | scGNN, scCDCG, AttentionAE_sc | GNN methods; env-blocked pending dependency check |
| **Foundation** | scGPT, GeneFormer, GeneCompass | Large model placeholders; not auto-run |
| **Reference** | OtherMode/scCluBench-main | Original scCluBench repository; read-only reference |

## Standard Outputs

Every runnable model emits the following files to its `--save_dir`:

| File | Description |
|---|---|
| `embedding_final.npy` | Final embedding (n_cells × n_z), float32 |
| `labels.npy` | Ground truth labels, int64 |
| `metrics.json` | 8 evaluation metrics |
| `args.json` | Running arguments |
| `embedding_{epoch}.npy` | Per-epoch embeddings |
| `metrics_{epoch}.json` | Per-epoch metrics |
| `embedding.h5` | HDF5 bundle of embeddings and predictions |
| `types_{epoch}_pred.csv` | Ground truth vs predicted labels |

## Quick Start

### Smoke Test (1 epoch, CPU)

```bash
# DEC
python methods/DeepLearning/dec/run.py \
  --data_path data/SRP182008.h5ad \
  --save_dir results/smoke_dec \
  --n_clusters 15 --pretrain_epochs 1 --epochs 1 --no_cuda

# scDCC
python methods/DeepLearning/scDCC/run.py \
  --data_path data/SRP182008.h5ad \
  --save_dir results/smoke_scdcc \
  --n_clusters 15 --pretrain_epochs 1 --epochs 1 --no_cuda

# scDSC (use smaller dataset for CPU)
python methods/GNN/scDSC/run.py \
  --data_path data/subsample_2k.h5ad \
  --save_dir results/smoke_scdsc \
  --n_clusters 7 --pretrain_epochs 1 --epochs 1 --no_cuda
```

### Standard Preprocessing

All models use `methods/preprocess.py::prepare_data_for_model()`:

```python
from preprocess import prepare_data_for_model
X, Y, sf, adata = prepare_data_for_model(
    'data/my_data.h5ad',
    size_factors=True,
    filter_min_counts=True,
    logtrans_input=True,
    normalize_input=True
)
# X: (n_cells, n_hvg) — scaled and normalized
# Y: cell type labels
# sf: size factors
# adata: full AnnData with layers['norm_log'] for ZINB loss
```

## Model Authenticity Policy

This benchmark **does not allow simplified replacement implementations**. Migrated methods must preserve the original scCluBench model architecture, losses, and training procedure unless explicitly documented as environment-gated or excluded.

### 本项目严禁为了跑通而使用简化模型、替代算法或伪实现。所有进入正式表格的方法必须通过模型真实性审计。

**Key principles:**
- No replacing deep models with KMeans, PCA, or plain AE
- No removing core losses (ZINB, KL divergence, clustering loss)
- No skipping training phases (pretrain, clustering)
- All formal benchmark results include `authenticity.json` with `substitute_model_used: false`

### Authenticity Status

| Status | Meaning | Eligible for Formal Table |
|--------|---------|------------------------|
| `VERIFIED` | Core architecture, losses, and training preserved | Yes (default) |
| `ENV-GATED` | Blocked by missing env (TensorFlow, etc.) | No |
| `PENDING` | Code migrated, audit not yet complete | No |
| `PLACEHOLDER` | Stub only | No |
| `FAILED` | Known issues | No |

Run `python scripts/audit_model_authenticity.py` to check all models.

See [docs/migration_status.md](docs/migration_status.md) and [docs/model_authenticity.md](docs/model_authenticity.md) for detailed checklists.

## Formal Benchmark

### Formal Method List

The default formal benchmark includes 10 methods:

| Role | Method | Key |
|------|--------|-----|
| **Proposed** | NeighborMix-scMAE | `neighbormix_scmae` |
| **Ablation** | NM-scMAE-noMix | `nm_scmae_nomix` |
| **External baseline** | scMAE | `scmae` |
| **Deep baseline** | DEC | `dec` |
| **Deep baseline** | scDCC | `scdcc` |
| **Deep baseline** | scDSC | `scdsc` |
| **Traditional** | ScanpyStandard | `scanpy_standard` |
| **Traditional** | Leiden | `leiden` |
| **Traditional** | Louvain | `louvain` |
| **Traditional** | SC3 | `sc3` |

### Quick Run

```bash
# Run all 10 formal methods (GPU enabled, GPU 1)
python scripts/run_formal_benchmark.py \
  --data_path data/subsample_2k.h5ad \
  --dataset_name subsample_2k \
  --out_dir results/formal \
  --n_clusters 7 \
  --seeds 42 \
  --gpu 1

# CPU-only run
python scripts/run_formal_benchmark.py \
  --data_path data/subsample_2k.h5ad \
  --dataset_name subsample_2k \
  --out_dir results/formal \
  --n_clusters 7 \
  --seeds 42 \
  --no_cuda

# Preflight / dry run (validate without executing)
python scripts/run_formal_benchmark.py \
  --data_path data/subsample_2k.h5ad \
  --dataset_name subsample_2k \
  --out_dir results/formal \
  --n_clusters 7 \
  --seeds 42 \
  --gpu 1 \
  --dry_run

# Run specific methods
python scripts/run_formal_benchmark.py \
  --data_path data/subsample_2k.h5ad \
  --methods neighbormix_scmae dec scdcc \
  --n_clusters 7 --seeds 42 --gpu 1
```

### GPU Policy

**GPU 0 and GPU 7 are FORBIDDEN** (occupied by other users). Default `--gpu` is `1`. Only `--gpu 1-6` or `--no_cuda` are allowed. The runner validates GPU policy at startup and will exit if a forbidden GPU is detected.

### Formal Results

Output per method run:
- `embedding_final.npy`, `labels.npy`, `metrics.json`, `args.json`
- `authenticity.json` (authenticity audit record, `substitute_model_used: false`)
- `status.json` (execution metadata: status, runtime, GPU, commit SHA, branch)
- `command.txt` (exact command executed)
- `run.log` (stdout/stderr capture)

Summary tables in output root:
- `benchmark_summary.csv` — full per-run traceability (seed, status, runtime, GPU, commit, command, error)
- `benchmark_summary_mean_std.csv` — mean ± std over seeds (only `status=success` + `authenticity=VERIFIED` + `substitute_model_used=false` rows included)
- `benchmark_summary_mean_std.md` — markdown version of the aggregated table

**Only runs where `status=success`, `authenticity=VERIFIED`, and `substitute_model_used=false` are eligible for the main paper table.**

## Migration Status

See [docs/migration_status.md](docs/migration_status.md) for detailed migration status of all baselines.

## Dependencies

| Package | Required By |
|---|---|
| PyTorch | DEC, scDCC, scDSC, NeighborMix_scMAE, scMAE |
| TensorFlow/Keras | scDeepCluster, scNAME, scziDesk, DESC |
| scanpy | All models (preprocessing) |
| scipy, scikit-learn | All models (evaluation, clustering) |

## Notes

- `OtherMode/scCluBench-main/` is reference-only. No runnable code in `methods/` imports from it.
- Foundation models (scGPT, GeneFormer, GeneCompass) are placeholders and are not included in default benchmark runs.
- GNN models (scDSC, scGNN, scCDCG) require sufficient memory for the k-NN adjacency matrix.
