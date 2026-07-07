# Cross-Branch Model Benchmark Audit (2026-07-04)

This report consolidates benchmark evidence found across remote branches without copying raw model weights or full run directories. Ordinary `scMAE` rows are intentionally excluded; `varfloor_scmae` is retained because it corresponds to `scMAE_DEC_StdFloor`.

## Source Branches

- `origin/scmae-refactor`: unified protocol evals, backfill summary/log, and failed `scVI` runs.
- `origin/main`: existing VarFloor/PCA completion report, PCA objective-swap pointer, and `results_master_20260619` smoke record for `scCDCG`.
- `origin/codex/plantspade-lgcl-gated-fusion`: model implementation diff evidence only; no results were copied from this branch.

## Model Coverage

| model | status | datasets | records | branches | notes |
| --- | --- | --- | ---: | --- | --- |
| `PhytoCluster` | success | 3: Mouse_Pancreas_1, SRP171040, SRP182008 | 6 | origin/scmae-refactor | See CSV for per-row metrics and source paths. |
| `scDeepCluster` | NO_RECORD |  | 0 |  | No benchmark result record found in searched benchmark/result paths. |
| `varfloor_scmae` | success | 23: Bach, Baron, CRA002977_1, CRA007122, Guo, Limb_Muscle, Macosko, Melanoma_5K ... | 23 | origin/main | See CSV for per-row metrics and source paths. |
| `scVI` | failed | 3: Mouse_Pancreas_1, SRP171040, SRP182008 | 3 | origin/scmae-refactor | See CSV for per-row metrics and source paths. |
| `GeneCompass` | NO_RECORD |  | 0 |  | No benchmark result record found in searched benchmark/result paths. |
| `GeneFormer` | NO_RECORD |  | 0 |  | No benchmark result record found in searched benchmark/result paths. |
| `scCDCG` | completed_no_status_json:1, success | 4: Bach, Mouse_Pancreas_1, SRP171040, SRP182008 | 7 | origin/main, origin/scmae-refactor | See CSV for per-row metrics and source paths. |
| `scGNN` | success | 3: Mouse_Pancreas_1, SRP171040, SRP182008 | 6 | origin/scmae-refactor | See CSV for per-row metrics and source paths. |
| `PCA_KMeans` | success | 23: Bach, Baron, CRA002977_1, CRA007122, Guo, Limb_Muscle, Macosko, Melanoma_5K ... | 32 | origin/main | See CSV for per-row metrics and source paths. |
| `ScanpyStandard` | NO_RECORD |  | 0 |  | No benchmark result record found in searched benchmark/result paths. |

## Key Comparable Results

- `PhytoCluster` has unified eval JSON records for Mouse_Pancreas_1, SRP171040, and SRP182008 on `origin/scmae-refactor`.
- `scCDCG` has unified/backfill records for Mouse_Pancreas_1, SRP171040, and SRP182008, plus a Bach smoke row in `origin/main` results master.
- `scGNN` has backfill smoke-passed rows for Mouse_Pancreas_1, SRP171040, and SRP182008 on `origin/scmae-refactor`.
- `scVI` has only failed run records for Mouse_Pancreas_1, SRP171040, and SRP182008; no metrics are reported for those failed runs.
- `varfloor_scmae` and `PCA_KMeans` each have 23-dataset aggregate rows in `formal_varfloor_pca_completion_20260704/260629全benchmark结果.with_varfloor_pca.csv`.

## Branch Implementation Comparison

- `origin/codex/plantspade-lgcl-gated-fusion` is the only branch found with substantive target-model implementation changes relative to `origin/main`.
- It modifies `methods/GNN/scGNN/run.py`, `methods/GNN/scGNN/scGNN.py`, and `methods/DeepLearning/scDeepCluster/run.py`, and expands the scDeepCluster source/materials under `methods/DeepLearning/scDeepCluster/`.
- Other manifest-only hits add CAAM-related registration in `methods/method_manifest.yaml`; they do not alter the target models in this audit.

## Files

- `cross_branch_benchmark_records.csv`: normalized benchmark/failed-run evidence.
- `branch_model_diff_summary.csv`: target model branch-diff audit.
- `cross_branch_benchmark_summary.md`: this summary.
