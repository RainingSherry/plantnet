# Safe-RDG-PCA v2 Benchmark-Compatible Progress

## Purpose

This v2 run corrects the Safe-RDG-PCA PCA backbone to match the 20260705
benchmark runner:

```text
methods/Traditional/PCA_KMeans/run.py
```

The previous Safe-RDG-PCA run used the correct data files, but its internal
PCA/KMeans protocol was not benchmark-compatible. Specifically, it z-scored
before PCA, z-scored after PCA, and used KMeans `n_init=50`. Therefore the
previous result is exploratory / invalid as a direct PCA benchmark comparison.

## Corrections

- `z_raw` now uses benchmark-compatible PCA:
  `PCA(n_components=dim, random_state=seed).fit_transform(data.astype(np.float64))`.
- No extra z-score is applied before or after benchmark PCA.
- `dim = max(2, min(raw_pca_dim, n_cells - 1, n_genes - 1))`.
- KMeans known-K defaults to `n_init=20`.
- RDG internal standardized PCA remains available only for gene-module
  eigengene/module computations.
- Full v2 results are written under this directory and do not overwrite:
  `experiment_reports/safe_rdg_pca_20260709`.

## Required Gate

Before full rerun, `pca_kmeans` must match 20260705
`pca_kmeans_known_k` for aligned dataset/seed metrics.

Compatibility command:

```bash
python experiment_reports/safe_rdg_pca_20260709_v2_benchmark_compatible/check_pca_compatibility.py \
  --new_run_root experiment_reports/safe_rdg_pca_20260709_v2_benchmark_compatible/compat_smoke_runs \
  --out_dir experiment_reports/safe_rdg_pca_20260709_v2_benchmark_compatible/analysis_compat_smoke \
  --datasets Bach,Pollen,Wang,worm_neuron_cell \
  --seeds 42,3047,3407
```

If compatibility fails, stop and diagnose before running the full suite.

## Compatibility Results

Smoke compatibility passed exactly:

```text
datasets = Bach,Pollen,Wang,worm_neuron_cell
seeds = 42,3047,3407
max_abs_diff(ARI/NMI/ACC) = 0.0
```

The first smoke attempt with the base Python environment failed on
`Bach seed=3407`. Diagnosis showed that PCA embeddings were nearly identical,
but KMeans differed between sklearn 1.8.0 and the 20260705 benchmark
environment. The v2 benchmark launcher therefore defaults to:

```text
/data/luolie/conda/envs/scclubench-main/bin/python
```

Full compatibility against the stored 20260705 result is exact for 45/48 runs.
The only differences are `hrvatin` seeds 42/3047/3407:

```text
strict tolerance = 1e-10
max_abs_ari_diff = 6.8509e-05
max_abs_nmi_diff = 1.2457e-04
max_abs_acc_diff = 4.1437e-05
```

Diagnosis:

- Old vs v2 `hrvatin` PCA embeddings differ by only `5.96e-08` max absolute
  value.
- Re-running the original `methods/Traditional/PCA_KMeans/run.py` now in the
  same `scclubench-main` environment reproduces the v2 `hrvatin` metrics.
- This is a tiny historical KMeans sensitivity/numerical drift in the stored
  20260705 `hrvatin` result, not a data-source error or a Safe-RDG PCA-backbone
  mismatch.

A relaxed diagnostic report with tolerance `1e-3` passes all 48/48:

```text
experiment_reports/safe_rdg_pca_20260709_v2_benchmark_compatible/analysis_full_relaxed_1e-3/pca_compatibility_report.json
```

The strict report is preserved here:

```text
experiment_reports/safe_rdg_pca_20260709_v2_benchmark_compatible/analysis_full/pca_compatibility_report.json
```

## Full Run Result

Completed and verified:

```text
VERIFY OK run_root=experiment_reports/safe_rdg_pca_20260709_v2_benchmark_compatible/runs runs=48 variants=10 analysis_checked=True
```

Analysis files:

- `analysis_full/all_runs.csv`
- `analysis_full/summary_by_variant.csv`
- `analysis_full/summary_by_variant_dataset.csv`
- `analysis_full/oracle_best_runs.csv`
- `analysis_full/summary_report.json`
- `analysis_full/pca_compatibility_report.csv`
- `analysis_full/status_full.csv`

Key full-suite metrics:

| Variant | Mean ARI | Median ARI | Mean regret vs PCA | Negative transfer | Graph activation |
| --- | ---: | ---: | ---: | ---: | ---: |
| `safe_rdg_calibrated_threshold` | 0.7173 | 0.7725 | 0.0110 | 0.0625 | 0.1875 |
| `pca_kmeans` | 0.6968 | 0.7477 | 0.0000 | 0.0000 | 0.0000 |
| `safe_rdg_heuristic` | 0.6635 | 0.6886 | 0.0867 | 0.7083 | 1.0000 |
| `rdg_always_on` | 0.6608 | 0.6792 | 0.0914 | 0.6875 | 1.0000 |
| `rdg_gene_only` | 0.6549 | 0.6729 | 0.0949 | 0.6875 | 1.0000 |
| `rdg_concat_kmeans` | 0.5790 | 0.6321 | 0.1376 | 0.8125 | 0.0000 |

Negative controls are harmful on average:

| Variant | Mean ARI | Mean regret vs PCA | Negative transfer |
| --- | ---: | ---: | ---: |
| `neg_random_cell_graph` | 0.4191 | 0.3355 | 0.7708 |
| `neg_shuffled_gene_cell_graph` | 0.1865 | 0.5416 | 0.9167 |
| `neg_degree_shuffle_graph` | 0.1526 | 0.5776 | 0.9375 |

## Success Criteria

- `safe_rdg_calibrated_threshold` mean/median ARI is higher than
  benchmark-compatible `pca_kmeans`.
- `safe_rdg_calibrated_threshold` has lower regret and lower negative transfer
  than `rdg_always_on`.
- Target rescue passes: `Pollen` or `worm_neuron_cell` best weak-dataset gain
  is approximately `0.416`.
- Graph activation is not zero: calibrated threshold enables graph in 9/48
  runs.
- `oracle_best` still exceeds PCA, but the mean oracle gain is now smaller:
  approximately `0.0756`.
- `safe_rdg_heuristic` fails as a main method: it activates graph in 48/48 runs
  and has high negative transfer.

## Current Judgment

After benchmark-compatible correction, PCA is much stronger than in the
exploratory v1 run. The RDG direction still has a measurable oracle and
calibrated-gate benefit, but always-on graph augmentation and the current
heuristic gate are not safe enough. The only defensible Safe-RDG result at this
stage is the cross-dataset `safe_rdg_calibrated_threshold` variant, with the
strict caveat that it is a meta-calibrated method rather than a fully
unsupervised one.
