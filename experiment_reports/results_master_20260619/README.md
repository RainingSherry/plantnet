# Results master table 2026-06-19

Generated from `results/` with symlinks followed.

## Files

- `all_results_master_table.csv`: primary table. It contains one row per run artifact plus summary-only rows that are not covered by a run artifact.
- `all_results_evidence_table.csv`: exhaustive evidence table. It keeps every compatible run artifact and every compatible summary CSV row, including duplicate aggregate rows.
- `dataset_model_status_summary.csv`: primary run-level aggregation by experiment group, dataset, and model key.
- `dataset_coverage_matrix.csv`: dataset-level coverage and status counts.
- `build_results_master.py`: reproducible generator.

## Counts

- total evidence rows: 4578
- primary master rows: 1971
- summary CSV evidence rows: 2908
- dataset/model summary rows: 611

## Primary Status Counts

| status | count |
| --- | ---: |
| `completed_no_status_json` | 1166 |
| `success` | 489 |
| `summary_only` | 292 |
| `failed` | 15 |
| `unknown` | 6 |
| `timeout` | 2 |
| `no_score_marker` | 1 |

## Primary Experiment Groups

| group | count |
| --- | ---: |
| `experiments/neighbormix_beta_mechanism_20260617` | 606 |
| `canonical/scmae_11datasets_20260609_12` | 396 |
| `experiments/neighbormix_ra_rg` | 343 |
| `experiments/neighbormix_stochastic_regularization_20260616` | 216 |
| `experiments/rc_nm_checkpoint_v4_1` | 204 |
| `canonical/formal_benchmark_20260607_08` | 117 |
| `experiments/cutaware_neighbormix_20260615` | 80 |
| `experiments/scgpt_plantnet_20260615` | 5 |
| `scratch/smoke` | 4 |

## Evidence Row Kinds

| kind | count |
| --- | ---: |
| `summary_csv` | 2908 |
| `run_artifact` | 1670 |
