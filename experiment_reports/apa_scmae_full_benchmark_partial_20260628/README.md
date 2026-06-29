# APA-scMAE full benchmark partial report (2026-06-28)

This package summarizes the APA-scMAE full benchmark results available at packaging time. The full benchmark target is 8 datasets x 3 seeds x 80 epochs. At packaging time, 22 of 24 dataset/seed runs had written metrics; CRA007122 seeds 2 and 3 were still pending/running.

## Included

- `tables/apa_scmae_full_metrics_by_run.csv`: current per-run metric table copied from benchmark output.
- `tables/apa_scmae_full_metrics_summary.csv`: current mean/std summary copied from benchmark output.
- `tables/current_dataset_metric_means.csv`: lightweight quicklook means by dataset and clustering method from available runs.
- `tables/completion_status.csv`: 24 expected dataset/seed pairs with completion status at packaging time.
- `manifests/task_manifest.json`: expected task manifest for the full APA-scMAE run.
- `related_desc_scname_fix/`: small CSV/log/script artifacts from the related DESC/scNAME migration fix report; large run artifacts are intentionally excluded.

## Current completion

- Expected dataset/seed runs: 24
- Runs present in metrics table: 22
- Metric rows: 44
- Missing/pending pairs: CRA007122:seed2, CRA007122:seed3

## Quicklook

The values are low across most datasets. This table is intentionally partial and should be replaced/augmented after CRA007122 seeds 2 and 3 finish.

| Dataset | Cluster method | n seeds | ACC mean | NMI mean | ARI mean | F1 macro mean |
|---|---:|---:|---:|---:|---:|---:|
| CRA002977_1 | kmeans_known_k | 3 | 0.2778 | 0.1013 | 0.0641 | 0.1900 |
| CRA002977_1 | leiden_fixed | 3 | 0.1197 | 0.1335 | 0.0246 | 0.1609 |
| CRA007122 | kmeans_known_k | 1 | 0.2815 | 0.1113 | 0.0794 | 0.2196 |
| CRA007122 | leiden_fixed | 1 | 0.1159 | 0.0895 | 0.0267 | 0.1386 |
| SRP145013 | kmeans_known_k | 3 | 0.2476 | 0.0731 | 0.0826 | 0.1632 |
| SRP145013 | leiden_fixed | 3 | 0.1370 | 0.0650 | 0.0238 | 0.1324 |
| SRP171040 | kmeans_known_k | 3 | 0.2170 | 0.1142 | 0.0542 | 0.1951 |
| SRP171040 | leiden_fixed | 3 | 0.2278 | 0.1252 | 0.0405 | 0.1860 |
| SRP182008 | kmeans_known_k | 3 | 0.2048 | 0.1226 | 0.0608 | 0.1703 |
| SRP182008 | leiden_fixed | 3 | 0.1794 | 0.1287 | 0.0432 | 0.1835 |
| SRP224648 | kmeans_known_k | 3 | 0.4318 | 0.1259 | 0.0446 | 0.3476 |
| SRP224648 | leiden_fixed | 3 | 0.1111 | 0.0965 | 0.0137 | 0.2325 |
| SRP235541 | kmeans_known_k | 3 | 0.1614 | 0.1131 | 0.0479 | 0.1304 |
| SRP235541 | leiden_fixed | 3 | 0.1509 | 0.1202 | 0.0393 | 0.1365 |
| SRP309176 | kmeans_known_k | 3 | 0.1633 | 0.0638 | 0.0263 | 0.1252 |
| SRP309176 | leiden_fixed | 3 | 0.1474 | 0.0715 | 0.0188 | 0.1241 |

## Excluded large artifacts

Large benchmark directories, embeddings, canonical h5ad files, model artifacts, and DESC/scNAME run artifacts were not included in this Git package to keep the repository lightweight.
