# Current Experiment Status

Last reviewed after Diffusion/DOLORIS cleanup.

## Stale Failure Reports

`failure_report.json` means the latest suite-level attempt failed. A stale failure report is an older failure marker in a run directory that later produced successful fixed-protocol outputs such as `*_fixed.csv` and `*_metrics.json`.

Stale reports were renamed to `stale_*failure*.json` so status checks no longer count them as active failures.

Current stale markers:

- `results/PlantSPADE_LGCL_protocol/SRP235541/plantspade_lgcl_support_attention/seed_1/stale_failure_report.json`
- `results/PlantSPADE_LGCL_protocol/SRP235541/plantspade_lgcl_support_attention/seed_1/stale_eval_support_attention_recovery_failure.json`
- `results/PlantSPADE_LGCL_protocol/SRP309176/scmae/seed_1/stale_failure_report.json`
- `results/PlantSPADE_LGCL_protocol/SRP309176/scmae/seed_2/stale_failure_report.json`

There are currently no active `failure_report.json` files under `results/PlantSPADE_LGCL_protocol`.

## Main Progress

Completed fixed-protocol seeds by dataset and method:

| Dataset | traditional_pca | phytocluster | scvi | scmae | lgcl_baseline | lgcl_support_attention | lgcl_gated_fusion |
|---|---:|---:|---:|---:|---:|---:|---:|
| SRP182008 | 5 | 5 | 5 | 5 | 5 | 5 | 3 |
| SRP235541 | 5 | 2 | 4 | 4 | 4 | 3 | 3 |
| SRP171040 | 5 | 2 | 3 | 3 | 5 | 5 | 3 |
| SRP309176 | 2 | 2 | 0 | 2 | 2 | 2 | 0 |
| SRP145013 | 2 | 2 | 0 | 2 | 2 | 2 | 0 |
| CRA002977_1 | 2 | 2 | 0 | 2 | 2 | 2 | 0 |
| SRP224648 | 3 | 2 | 3 | 3 | 3 | 3 | 3 |
| CRA007122 | 2 | 2 | 0 | 2 | 2 | 2 | 0 |

The most complete main dataset is `SRP182008`: all six configured main methods have five fixed-protocol seeds. `SRP171040` has five LGCL seeds but incomplete baseline seeds. `SRP235541` has partial baseline and LGCL completion. Other datasets are currently two- or three-seed progress runs.

## Important Tables

- `results/PlantSPADE_LGCL_protocol/tables/table_main_fixed_protocol.csv`: current aggregated fixed-protocol main table.
- `results/PlantSPADE_LGCL_protocol/tables/all_results_long.csv`: raw per-seed rows for all discovered fixed/oracle/sweep CSVs.
- `results/PlantSPADE_LGCL_protocol/tables/all_results_mean_std.csv`: grouped mean/std over all protocols and variants.
- `results/PlantSPADE_LGCL_protocol/tables/table_attention_ablation.csv`: support-attention and attention-ablation rows.
- `results/PlantSPADE_LGCL_protocol/tables/table_negative_sampling_ablation.csv`: negative-sampling ablation rows; currently empty because those runs are not present.
- `results/PlantSPADE_LGCL_protocol/tables/table_oracle_supplement.csv`: oracle/supplement rows, not for main reporting.
- `results/PlantSPADE_LGCL_protocol/tables/dataset_profiles_summary.csv`: dataset profile summary for all eight configured datasets.
- `results/PlantSPADE_LGCL_protocol/tables/comparison_8datasets_7methods_2seeds_fixed_protocol.csv`: curated two-seed fixed-protocol comparison.
- `results/PlantSPADE_LGCL_protocol/tables/table_main_8x7_summary.csv`: compact curated summary table.
- `results/PlantSPADE_LGCL_protocol/tables/small_main_experiment_4datasets_6methods_3seeds/`: four-dataset, six-method, three-seed subset tables.
- `results/PlantSPADE_LGCL_gated_fusion_comparison/fixed_protocol_summary.csv`: gated-fusion comparison summary.
- `results/PlantSPADE_LGCL_gated_fusion_comparison/gate_summary_by_seed.csv`: gated-fusion gate weight summary.
