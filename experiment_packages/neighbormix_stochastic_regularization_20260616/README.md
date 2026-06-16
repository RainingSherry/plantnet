# NeighborMix Stochastic Regularization Package

Package date: 2026-06-16

This package archives the NeighborMix stochastic regularization and alternative-neighborhood experiment in a GitHub-friendly form. It intentionally stores lightweight evidence and code snapshots, not large model or embedding artifacts.

## Scope

Experiment id:

```text
neighbormix_stochastic_regularization_20260616
```

Full run location on the local machine:

```text
results/experiments/neighbormix_stochastic_regularization_20260616
```

Completed matrix:

```text
8 datasets x 3 seeds x 9 methods = 216 runs
```

Datasets:

```text
Core negative-transfer: Tosches, Macosko, worm_neuron_cell
Core positive-gain: Melanoma_5K, Shekhar
Core neutral/stable: Guo
Optional validation: Wang, Pollen
```

Seeds:

```text
0, 1, 2
```

GPU policy:

```text
GPU 0 and GPU 7 were not used.
Runs were scheduled on GPU 1-6.
```

## Contents

```text
code/
  methods/DeepLearning/NeighborMix_scMAE/run_stochastic_ablation.py
  scripts/run_neighbormix_stochastic_ablation.py
  scripts/summarize_neighbormix_stochastic_ablation.py

results/
  summaries/
    main_results.csv
    group_summary.csv
    neighbor_diagnostics.csv
    interpretation.md
  run_artifacts/
    <dataset>/<method>/seed<seed>/
      args.json
      eval_fixed.csv
      eval_metrics.json
      metrics.json
      summary.json
      neighbor_diagnostics*.json/csv
      training_history.json
      per_class_neighbor_purity.csv

notes/
  experiment_readme.md

MANIFEST.files
MANIFEST.sha256
```

Excluded by design:

```text
embedding_final.npy
embeddings_base.npy
labels.npy
gene_names.npy
embedding.h5
adata*.h5ad
model.pt
stdout.log
stderr.log
```

Those files remain in the full local results directory and are too large or too operational for ordinary Git storage.

## Main Results

The required analysis outputs are in:

```text
results/summaries/main_results.csv
results/summaries/group_summary.csv
results/summaries/neighbor_diagnostics.csv
results/summaries/interpretation.md
```

The final interpretation is also copied from:

```text
notes/experiment_readme.md
```

Short conclusion:

```text
random_beta_uniform_0.1 was the best stochastic variant and beat fixed NeighborMix on 5/6 core datasets, but it failed the continuation screen because worst-case delta ARI became worse than fixed NeighborMix.

SNN was the best alternative-neighbor variant, but it did not beat fixed NeighborMix on negative-group mean delta ARI. Mutual/SNN/consensus did not provide enough evidence that replacing vanilla PCA-cosine KNN solved negative transfer.

Global random neighbors were worse than local stochastic variants, so useful regularization still depends on local neighborhood structure.
```

## Reproduction

Run the full matrix:

```bash
python scripts/run_neighbormix_stochastic_ablation.py \
  --seeds 0,1,2 \
  --datasets Tosches,Macosko,worm_neuron_cell,Melanoma_5K,Shekhar,Guo,Wang,Pollen \
  --gpus 1,2,3,4,5,6 \
  --epochs 80 \
  --batch_size 256
```

Regenerate summaries:

```bash
python scripts/summarize_neighbormix_stochastic_ablation.py \
  --root results/experiments/neighbormix_stochastic_regularization_20260616
```

## Verification

Final checks performed before packaging:

```text
main_results.csv rows: 216
neighbor_diagnostics.csv rows: 216
failed run markers: 0
running run markers: 0
Python compile check: passed for all three new scripts
Package size: about 15 MB
```
