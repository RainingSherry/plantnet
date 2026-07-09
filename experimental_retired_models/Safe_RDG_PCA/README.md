# Safe-RDG-PCA

Risk-aware graph augmentation for a strong PCA+KMeans anchor.

This experimental runner is intentionally classical: no GNN, no DEC,
no autoencoder, and no default expression residual smoothing.  It tests
whether reliable cell graphs and bootstrap-stable gene module graphs can
improve PCA+KMeans without negative transfer.

Primary variants:

- `pca_kmeans`
- `pca_spectral_kmeans`
- `rdg_cell_only`
- `rdg_gene_only`
- `rdg_concat_kmeans`
- `rdg_always_on`
- `safe_rdg_heuristic`

Calibrated gate variants are synthesized by
`experiment_reports/safe_rdg_pca_20260709/calibrate_gate.py` from Stage A
outputs, so the expensive graph construction is not repeated.

Optional negative controls are available with
`--include_negative_controls true`:

- `neg_random_cell_graph`
- `neg_degree_shuffle_graph`
- `neg_shuffled_gene_cell_graph`

Data handling:

- The runner uses `methods.DeepLearning.scMAE_family.load_scmae_dataset`,
  matching the current PCA_KMeans/scMAE preprocessing path:
  source selection -> normalize/log1p when raw -> HVG -> optional scaling ->
  label encoding.
- `experiment_reports/safe_rdg_pca_20260709/run_benchmark.py` discovers
  datasets from `data/scMAE` and the existing all-methods catalog.  By default
  it reuses processed h5ad files in
  `result/scmae_all_methods_20260705_full/converted_data`.
- If `--data_preference raw` is used, legacy scMAE `X/Y` h5 files under
  `data/scMAE` are first materialized as per-run `_prepared_input/*.h5ad`,
  then passed through the same `scMAE_family` preprocessing path.
