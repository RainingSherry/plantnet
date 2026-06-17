# NeighborMix Beta Mechanism Package

Package id:

```text
neighbormix_beta_mechanism_20260617
```

Local full result root:

```text
results/experiments/neighbormix_beta_mechanism_20260617
```

This package archives the staged stochastic beta NeighborMix mechanism experiment. It stores code snapshots, summary tables, interpretation notes, and lightweight per-run evidence. Large runtime artifacts are intentionally excluded.

## Completed Matrix

Completed runs:

```text
Stage 1: 6 datasets x 3 seeds x 8 variants  = 144 runs
Stage 2: 6 datasets x 3 seeds x 5 variants  = 90 runs
Stage 3: 6 datasets x 3 seeds x 4 variants  = 72 runs
Stage 4: 3 datasets x 3 seeds x 12 variants = 108 runs
Full:    8 datasets x 3 seeds x 8 variants  = 192 runs
Total: 606 runs
Failed/running markers: 0
```

GPU policy:

```text
GPU 0 and GPU 7 were not used.
Runs were scheduled on GPU 1-6.
```

## Contents

```text
code/
  methods/DeepLearning/NeighborMix_scMAE/run_beta_mechanism.py
  scripts/run_beta_mechanism.py
  scripts/summarize_beta_mechanism.py

results/
  summaries/
    main_results.csv
    group_summary.csv
    stage1_beta_mean_vs_randomness.csv
    stage2_beta_variance.csv
    stage3_local_mix_mechanism.csv
    stage4_bad_edge_robustness.csv
    full_benchmark_summary.csv
    neighbor_diagnostics.csv
    interpretation.md
  run_artifacts/
    <stage>/<dataset>/<variant>/seed<seed>/
      args.json
      eval_fixed.csv
      eval_metrics.json
      metrics.json
      summary.json
      neighbor_diagnostics*.json/csv
      training_history.json
      per_class_neighbor_purity.csv
      dataset_profile.json
      preprocess_config.json
      command.txt

notes/
  experiment_readme.md
  paper_claim_tracker.md
```

Excluded by design:

```text
*.npy
*.h5
*.h5ad
*.pt
*.pth
stdout.log
stderr.log
```

## Main Decision

The experiment supports stochastic beta NeighborMix as an average-performance mechanism, but not as an unconditional robust-improvement claim.

Key evidence:

```text
Stage 1 same-mean control:
random_beta_uniform_0.1 - fixed_beta_0.05 = +0.0638 mean ARI, +0.0257 macro-F1.

Stage 2 distribution control:
truncated_normal_beta_mean0.05_std0.02 was best at mean beta 0.05.

Full benchmark:
truncated_normal_beta_mean0.05_std0.02 mean ARI = 0.7270.
random_beta_uniform_0.1 mean ARI = 0.7234.
fixed_beta_0.1 mean ARI = 0.6828.
fixed_beta_0.05 mean ARI = 0.6748.
nm_scmae_nomix mean ARI = 0.6664.
```

Caveat:

```text
worm_neuron_cell remains a failure case.
random_beta_uniform_0.1: 0.4121 ARI vs noMix 0.5167.
truncated_normal_beta_mean0.05_std0.02: 0.4064 ARI vs noMix 0.5167.
```

## Reproduction

Run each stage:

```bash
python scripts/run_beta_mechanism.py --stage stage1 --gpus 1,2,3,4,5,6
python scripts/run_beta_mechanism.py --stage stage2 --gpus 1,2,3,4,5,6
python scripts/run_beta_mechanism.py --stage stage3 --gpus 1,2,3,4,5,6
python scripts/run_beta_mechanism.py --stage stage4 --gpus 1,2,3,4,5,6
python scripts/run_beta_mechanism.py --stage full --gpus 1,2,3,4,5,6
```

Regenerate summaries:

```bash
python scripts/summarize_beta_mechanism.py \
  --root results/experiments/neighbormix_beta_mechanism_20260617
```
