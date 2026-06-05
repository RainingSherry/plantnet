# Backfill Smoke-Passed Models

Scope: dec, scDCC, scDSC, scGNN, scCDCG, AttentionAE-sc on Mouse_Pancreas_1, SRP171040, SRP182008.

Notes:
- AttentionAE-sc used n_neighbors=5. SRP171040 additionally used n_heads=1 because the default 8-head attention ran out of memory on 33956 cells.
- SRP171040/scDSC was retried successfully with low BLAS thread counts after an OpenBLAS segfault.
- Mouse_Pancreas_1/scGNN was retried successfully after automatic label-column detection was fixed.

## Best Row Per Dataset/Method

| dataset          | method         | cluster_method | acc    | nmi    | ari    | f1_macro | n_pred_clusters |
| ---------------- | -------------- | -------------- | ------ | ------ | ------ | -------- | --------------- |
| Mouse_Pancreas_1 | scDSC          | kmeans         | 0.6013 | 0.6566 | 0.4789 | 0.3712   | 13              |
| Mouse_Pancreas_1 | scDCC          | kmeans         | 0.4745 | 0.6746 | 0.3764 | 0.3897   | 13              |
| Mouse_Pancreas_1 | dec            | leiden         | 0.4931 | 0.6411 | 0.358  | 0.3404   | 11              |
| Mouse_Pancreas_1 | scGNN          | leiden         | 0.4417 | 0.622  | 0.3441 | 0.2928   | 12              |
| Mouse_Pancreas_1 | scCDCG         | kmeans         | 0.3998 | 0.4649 | 0.2492 | 0.1847   | 13              |
| Mouse_Pancreas_1 | AttentionAE-sc | kmeans         | 0.4109 | 0.4351 | 0.1839 | 0.278    | 13              |
| SRP171040        | dec            | kmeans         | 0.5108 | 0.6263 | 0.4385 | 0.4428   | 12              |
| SRP171040        | scDSC          | kmeans         | 0.471  | 0.6003 | 0.3958 | 0.439    | 12              |
| SRP171040        | scGNN          | kmeans         | 0.502  | 0.5278 | 0.3593 | 0.4564   | 12              |
| SRP171040        | scCDCG         | kmeans         | 0.4142 | 0.4689 | 0.3012 | 0.3764   | 12              |
| SRP171040        | AttentionAE-sc | kmeans         | 0.4127 | 0.4804 | 0.2886 | 0.344    | 12              |
| SRP171040        | scDCC          | kmeans         | 0.4031 | 0.4686 | 0.2595 | 0.3686   | 12              |
| SRP182008        | dec            | kmeans         | 0.5367 | 0.5767 | 0.4168 | 0.4654   | 15              |
| SRP182008        | scGNN          | leiden         | 0.5391 | 0.5642 | 0.3692 | 0.5801   | 20              |
| SRP182008        | scDSC          | leiden         | 0.4879 | 0.5874 | 0.3466 | 0.5637   | 26              |
| SRP182008        | AttentionAE-sc | leiden         | 0.345  | 0.4003 | 0.2014 | 0.3234   | 39              |
| SRP182008        | scDCC          | leiden         | 0.3579 | 0.3816 | 0.2007 | 0.3617   | 24              |
| SRP182008        | scCDCG         | leiden         | 0.3139 | 0.3509 | 0.1832 | 0.278    | 14              |

## Files

- benchmarks/unified_protocol/metrics_long.csv
- benchmarks/unified_protocol/backfill_smoke_passed_best.csv
- benchmarks/unified_protocol/backfill_smoke_passed_all.csv
- results/backfill_smoke_passed_best.csv
- results/backfill_smoke_passed_all.csv
