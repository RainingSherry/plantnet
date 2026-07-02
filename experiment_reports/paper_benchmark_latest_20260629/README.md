# Latest paper benchmark 2026-06-29

This package consolidates the latest non-ablation benchmark results from `results/` and `experiment_reports/`.

## Included sources

- `results/canonical/formal_benchmark_20260607_08`
- `results/canonical/scmae_11datasets_20260609_12`
- `results/experiments/scgpt_plantnet_20260615`
- `experiment_reports/desc_scname_fix_benchmark_20260628/desc_scname_updated_values_all_runs.csv`
- `experiment_reports/apa_scmae_full_benchmark_partial_20260628/tables/apa_scmae_full_metrics_by_run.csv`

## Exclusion rule

Ablation/mechanism/smoke result groups are excluded, including beta mechanism, stochastic regularization, cut-aware, RA/RG sensitivity, RC checkpoint, and APA v2 ablation outputs.

APA_scMAE `kmeans_known_k` is used as the primary model readout in the paper table. The `leiden_fixed` readout is kept in the evidence table only.

## Counts

- evidence rows: 586
- paper candidate rows: 564
- paper formatter usable scored rows: 531
- evidence datasets: 24
- paper datasets: 24
- paper scored datasets: 23
- paper methods: 15
- paper scored methods: 15

## Files

- `latest_benchmark_evidence_all_rows.csv`: all selected non-ablation evidence rows, including skipped runs and APA secondary readouts.
- `latest_benchmark_paper_runs.csv`: rows used to build the paper benchmark tables.
- `paper_tables/`: numeric, formatted CSV, LaTeX, and XLSX benchmark tables.
- `latest_benchmark_status_by_dataset_method.csv`: status counts.
- `latest_benchmark_coverage_matrix.csv`: dataset/method coverage.
