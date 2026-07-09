# Safe-RDG-PCA Progress

## Benchmark-Compatibility Notice

This run is now marked exploratory / invalid as a direct comparison to the
20260705 `pca_kmeans_known_k` benchmark. The data source was correct
(`result/scmae_all_methods_20260705_full/converted_data`), but the Safe-RDG-PCA
internal PCA/KMeans protocol was not benchmark-compatible: it z-scored before
PCA, z-scored after PCA, and used KMeans `n_init=50`.

The corrected benchmark-compatible rerun is tracked separately under:

```text
experiment_reports/safe_rdg_pca_20260709_v2_benchmark_compatible
```

## Current Implementation

- Model runner: `experimental_retired_models/Safe_RDG_PCA/run.py`
- Benchmark launcher: `experiment_reports/safe_rdg_pca_20260709/run_benchmark.py`
- Gate calibration: `experiment_reports/safe_rdg_pca_20260709/calibrate_gate.py`
- Summary/reporting: `experiment_reports/safe_rdg_pca_20260709/summarize.py`
- Output verifier: `experiment_reports/safe_rdg_pca_20260709/verify_outputs.py`

The runner now uses the same preprocessing entry point as the current
PCA_KMeans/scMAE methods: `methods.DeepLearning.scMAE_family.load_scmae_dataset`.
By default, the benchmark launcher discovers datasets from `data/scMAE` and
reuses processed h5ad files from
`result/scmae_all_methods_20260705_full/converted_data`.

## Implemented Variants

Stage A:

- `pca_kmeans`
- `pca_spectral_kmeans`
- `rdg_cell_only`
- `rdg_gene_only`
- `rdg_concat_kmeans`
- `rdg_always_on`
- `safe_rdg_heuristic`

Optional negative controls:

- `neg_random_cell_graph`
- `neg_degree_shuffle_graph`
- `neg_shuffled_gene_cell_graph`

Negative controls are built lazily. They are not constructed unless
`--include_negative_controls`, `negative_controls_all`, or a specific
negative-control variant is requested.

Calibrated variants are synthesized after Stage A:

- `safe_rdg_calibrated_threshold`
- `safe_rdg_calibrated_logistic_explore`

## Smoke Evidence

Target smoke:

```bash
python experiment_reports/safe_rdg_pca_20260709/run_benchmark.py \
  --datasets Pollen,worm_neuron_cell \
  --seeds 42 \
  --max_workers 1 \
  --gene_bootstrap_B 5 \
  --include_negative_controls \
  --out_root experiment_reports/safe_rdg_pca_20260709/smoke_target_runs \
  --no_resume
```

Analysis output:

- `experiment_reports/safe_rdg_pca_20260709/analysis_target_smoke`

Key result: `safe_rdg_heuristic` improved both target datasets vs PCA.
Mean delta ARI vs PCA was approximately `+0.440`.

Strong-dataset smoke:

```bash
python experiment_reports/safe_rdg_pca_20260709/run_benchmark.py \
  --datasets Wang,Limb_Muscle \
  --seeds 42 \
  --max_workers 1 \
  --gene_bootstrap_B 5 \
  --out_root experiment_reports/safe_rdg_pca_20260709/smoke_strong_runs \
  --no_resume
```

Analysis output:

- `experiment_reports/safe_rdg_pca_20260709/analysis_strong_smoke`

Key result: `safe_rdg_heuristic` and `safe_rdg_calibrated_threshold` were
non-inferior to PCA on these two strong-dataset probes. No ARI drop below
`-0.03` was observed in this smoke subset.

Formal-parameter pilot:

```bash
python experiment_reports/safe_rdg_pca_20260709/run_benchmark.py \
  --datasets Pollen,Wang \
  --seeds 42 \
  --max_workers 1 \
  --gene_bootstrap_B 20 \
  --include_negative_controls \
  --out_root experiment_reports/safe_rdg_pca_20260709/pilot_b20_runs \
  --no_resume
```

Analysis output:

- `experiment_reports/safe_rdg_pca_20260709/analysis_pilot_b20`

Verifier:

```bash
python experiment_reports/safe_rdg_pca_20260709/verify_outputs.py \
  --run_root experiment_reports/safe_rdg_pca_20260709/pilot_b20_runs \
  --analysis_dir experiment_reports/safe_rdg_pca_20260709/analysis_pilot_b20 \
  --include_negative_controls \
  --require_analysis
```

Result: verifier passed for 2 dataset-seed runs and 10 variants. In this
B=20 pilot, `safe_rdg_calibrated_threshold` and `rdg_always_on` had mean
ARI approximately `0.920`, while `safe_rdg_heuristic` had mean ARI
approximately `0.905`. Negative controls were harmful on average.

## Current Interpretation

- The dual graph direction has real rescue potential on `Pollen` and
  `worm_neuron_cell`.
- `rdg_gene_only` is the strongest Stage A branch on the full 16-dataset
  suite, while `rdg_cell_only` and PCA spectral can be harmful. This supports
  keeping the branch-level ablations rather than presenting a single fused
  graph result.
- The heuristic gate activates graphs on every full-suite run. It improves over
  PCA on average, but it does not reduce mean regret relative to always-on.
- The cross-dataset calibrated threshold gate is currently the safer main
  candidate: it activates graphs in 25/48 runs, has positive enabled-graph mean
  gain, and reduces both mean regret and negative transfer relative to
  `rdg_always_on`.
- Negative controls are harmful on average, especially the shuffled
  gene-module graph, which supports that the true gene-module branch is not
  merely benefiting from arbitrary graph smoothing.

## Full Run Command

```bash
python experiment_reports/safe_rdg_pca_20260709/run_benchmark.py \
  --datasets all \
  --seeds 42,3047,3407 \
  --max_workers 1 \
  --gene_bootstrap_B 20 \
  --include_negative_controls \
  --out_root experiment_reports/safe_rdg_pca_20260709/runs \
  --no_resume

python experiment_reports/safe_rdg_pca_20260709/calibrate_gate.py \
  --run_root experiment_reports/safe_rdg_pca_20260709/runs \
  --out_dir experiment_reports/safe_rdg_pca_20260709/analysis_full

python experiment_reports/safe_rdg_pca_20260709/summarize.py \
  --run_root experiment_reports/safe_rdg_pca_20260709/runs \
  --out_dir experiment_reports/safe_rdg_pca_20260709/analysis_full

python experiment_reports/safe_rdg_pca_20260709/verify_outputs.py \
  --run_root experiment_reports/safe_rdg_pca_20260709/runs \
  --analysis_dir experiment_reports/safe_rdg_pca_20260709/analysis_full \
  --datasets all \
  --seeds 42,3047,3407 \
  --include_negative_controls \
  --require_analysis
```

## Full Run Result

Completed and verified:

```text
VERIFY OK run_root=experiment_reports/safe_rdg_pca_20260709/runs runs=48 variants=10 analysis_checked=True
```

Analysis output:

- `experiment_reports/safe_rdg_pca_20260709/analysis_full/all_runs.csv`
- `experiment_reports/safe_rdg_pca_20260709/analysis_full/summary_by_variant.csv`
- `experiment_reports/safe_rdg_pca_20260709/analysis_full/summary_by_variant_dataset.csv`
- `experiment_reports/safe_rdg_pca_20260709/analysis_full/oracle_best_runs.csv`
- `experiment_reports/safe_rdg_pca_20260709/analysis_full/summary_report.json`

Key full-suite metrics:

| Variant | Mean ARI | Median ARI | Mean regret vs PCA | Negative transfer | Graph activation |
| --- | ---: | ---: | ---: | ---: | ---: |
| `pca_kmeans` | 0.4983 | 0.5149 | 0.0000 | 0.0000 | 0.0000 |
| `rdg_gene_only` | 0.6776 | 0.6906 | 0.0215 | 0.2083 | 1.0000 |
| `rdg_always_on` | 0.6363 | 0.6854 | 0.0470 | 0.2917 | 1.0000 |
| `safe_rdg_heuristic` | 0.5892 | 0.6847 | 0.0827 | 0.2708 | 1.0000 |
| `safe_rdg_calibrated_threshold` | 0.6430 | 0.6854 | 0.0041 | 0.0417 | 0.5208 |

Success criteria check:

- Basic success: passed for `safe_rdg_heuristic` and
  `safe_rdg_calibrated_threshold`.
- Non-inferiority on PCA-strong datasets: passed for
  `safe_rdg_calibrated_threshold`; heuristic had 2 drops larger than 0.03.
- Target rescue: passed; the best weak-dataset gain was approximately 0.751.
- Gate contribution: passed for `safe_rdg_calibrated_threshold`; not passed for
  heuristic because mean regret was not lower than always-on.
- Direction decision: continue the RDG line, but use the calibrated threshold
  gate or a redesigned unsupervised gate as the paper-facing safety mechanism.
