# Formal VarFloor/PCA Completion, 2026-07-04

This directory records the formal completion run requested for
`results/260629全benchmark结果.csv`.

Methods:

- `varfloor_scmae`: VarFloor-scMAE, implemented at
  `methods/DeepLearning/scMAE_DEC_StdFloor/run.py`
- `pca_kmeans_known_k`: PCA+KMeans known-K, implemented at
  `methods/Traditional/PCA_KMeans/run.py`

Protocol:

- 23 datasets already present in `results/260629全benchmark结果.csv`
- seeds: `42,2024,3407`
- known-K evaluation with labels used only for final scoring
- run artifacts:
  `results/canonical/formal_varfloor_scmae_pca_20260704`

Outputs:

- `task_manifest.json`: full task list
- `formal_completion_runs.csv`: 138 successful dataset/method/seed runs
- `formal_completion_summary.csv`: 46 aggregate rows
- `260629全benchmark结果.with_varfloor_pca.csv`: merged CSV copy

The reference CSV was updated after backing up the previous file to
`results/260629全benchmark结果.csv.bak_20260704_155211`.

