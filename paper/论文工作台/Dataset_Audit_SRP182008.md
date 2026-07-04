# Dataset Audit: SRP182008

Updated: 2026-06-26

## 1. Local File

```text
数据文件/SRP182008.h5ad
```

File size:

```text
94,738,395 bytes
```

## 2. Source Identification

The local dataset is labeled as:

```text
Dataset: SRP182008
Orig.ident: SRX5290443
Organ: Root
Tissue: Root tip
ACE: 10 days old seedling
Condition: Normal
Genotype: Col-0
Libraries: 10x Genomics
```

This matches the Arabidopsis root single-cell RNA-seq dataset from:

```text
Zhang et al. 2019, Molecular Plant
A Single-Cell RNA Sequencing Profiles the Developmental Landscape of Arabidopsis Root
PMID: 31004836
DOI: 10.1016/j.molp.2019.04.004
```

BibTeX key added:

```text
zhang2019arabidopsisroot
```

Plant database context:

```text
scPlantDB provides a plant single-cell atlas database and includes plant atlas resources across species.
BibTeX key: he2024scplantdb
```

## 3. Matrix Structure

The expression matrix is stored as a CSR sparse matrix:

```text
X encoding-type: csr_matrix
X shape: 13,514 cells × 53,678 genes
X dtype: float64
X nnz: 17,378,932
X density: 0.0239575743
X sparsity: 0.9760424257
```

The raw matrix has the same shape and nnz:

```text
raw/X shape: 13,514 × 53,678
raw/X nnz: 17,378,932
raw/X density: 0.0239575743
```

Implication:

```text
This is an extremely sparse plant scRNA-seq dataset. It is a good development dataset for testing whether mask strategies collapse into weak zero-to-zero corruption.
```

## 4. Observation Fields

Available `obs` fields:

```text
ACE
Celltype
Condition
Dataset
Genotype
Libraries
Organ
Orig.ident
Percent.mt
Seurat_clusters
Tissue
nCount_RNA
nFeature_RNA
```

Potential evaluation labels:

```text
Celltype
Seurat_clusters
```

Do not use these fields during training.

Potential technical covariates:

```text
nCount_RNA
nFeature_RNA
Percent.mt
```

Fields with only one value in this file:

```text
Orig.ident = SRX5290443
ACE = 10 days old seedling
Condition = Normal
Genotype = Col-0
Libraries = 10x Genomics
Organ = Root
Tissue = Root tip
Dataset = SRP182008
```

Implication:

```text
This local file does not contain multiple batches or conditions. It is suitable for Stage 1 smoke testing and plant single-dataset development, but it cannot test batch robustness by itself.
```

## 5. Cell Type Distribution

`Celltype` has 15 categories:

| Celltype | Count | Fraction |
|---|---:|---:|
| Root stele | 2569 | 0.1901 |
| Root cortex | 1751 | 0.1296 |
| Phloem/Pericycle | 1524 | 0.1128 |
| Root hair | 1153 | 0.0853 |
| Columella root cap | 1082 | 0.0801 |
| Lateral root cap | 948 | 0.0701 |
| Root endodermis | 854 | 0.0632 |
| Unknow | 805 | 0.0596 |
| Companion cell | 690 | 0.0511 |
| S phase | 478 | 0.0354 |
| Protoxylem | 396 | 0.0293 |
| G2/M phase | 354 | 0.0262 |
| Sieve element | 344 | 0.0255 |
| Metaxylem | 300 | 0.0222 |
| Root epidermis | 266 | 0.0197 |

Notes:

```text
The label "Unknow" appears in the data and should be preserved as-is unless a preprocessing decision explicitly remaps it.
The smallest labeled cell type is Root epidermis with 266 cells, about 1.97%.
```

## 6. Seurat Cluster Distribution

`Seurat_clusters` has 24 clusters.

Largest clusters:

| Cluster | Count | Fraction |
|---|---:|---:|
| 0 | 1524 | 0.1128 |
| 1 | 1480 | 0.1095 |
| 2 | 1477 | 0.1093 |
| 3 | 1092 | 0.0808 |
| 4 | 1082 | 0.0801 |
| 5 | 809 | 0.0599 |
| 6 | 690 | 0.0511 |
| 7 | 571 | 0.0423 |
| 8 | 511 | 0.0378 |
| 9 | 396 | 0.0293 |

Smallest cluster:

```text
Cluster 23: 71 cells, 0.53%
```

Implication:

```text
If Seurat_clusters is used only as an evaluation target, this dataset contains at least one rare cluster below 1%, useful for rare-cluster recovery diagnostics.
```

## 7. Per-Cell Expression Summary

Per-cell nonzero genes:

```text
min: 234
q25: 646
median: 1007
mean: 1285.99
q75: 1714
max: 3812
```

`nCount_RNA`:

```text
min: 500
q25: 968
median: 1889.5
mean: 3255.99
q75: 4166
max: 39150
```

`nFeature_RNA`:

```text
min: 234
q25: 646
median: 1007
mean: 1285.99
q75: 1714
max: 3812
```

`Percent.mt`:

```text
min: 0
median: 0
mean: 0.000178
max: 0.151745
```

## 8. Embeddings Already Present

The file contains:

```text
obsm/X_umap: 13,514 × 2
obsm/X_tsne: 13,514 × 3
```

These should not be used for model training. They may be useful only for sanity-check visualization.

## 9. Stage 1 Recommendation

This dataset is appropriate for:

```text
scMAE-compatible smoke test
mask diagnostics under extreme sparsity
plant single-cell development experiment
rare Seurat cluster recovery analysis
```

This dataset is not sufficient for:

```text
batch robustness
cross-condition generalization
multi-tissue benchmark claims
formal top-journal performance claims
```

## 10. Immediate Decisions

Recommended Stage 1 setup:

```text
input genes: HVG 2000 first, full genes only after memory check
evaluation labels: Celltype and Seurat_clusters, evaluation only
mask ratios: 0.2, 0.3, 0.4
seeds: 3 development seeds
primary diagnostics: zero_to_zero_fraction, effective_changed_fraction, actual_mask_ratio_observed
```

