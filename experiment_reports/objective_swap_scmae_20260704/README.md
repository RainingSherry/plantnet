# Objective-swap scMAE experiment

This directory contains external research experiments. It intentionally lives
under `experiment_reports/`, not under `methods/`, because these runners test
post-scMAE hypotheses rather than formal benchmark methods.

Question:

```text
Is the reconstruction-first scMAE objective the wrong upper-level assumption
for fine-grained scRNA clustering?
```

First prototype:

```text
input -> MLP encoder -> latent z -> prototype assignment
```

There is no decoder and no masked gene reconstruction. Training uses only
unlabeled objectives:

- two weak gene-expression views
- cross-view prototype assignment consistency
- sharpened or Sinkhorn targets
- latent std-floor anti-collapse
- mild entropy/confidence regularization

Labels are not used during training, early stopping, or model selection. They
are used only after training for fixed-K evaluation.

Run one Macosko diagnostic:

```bash
python experiment_reports/objective_swap_scmae_20260704/run_cluster_first.py \
  --data_path methods/DeepLearning/scMAEs/benchmark_data/Macosko.h5ad \
  --save_dir experiment_reports/objective_swap_scmae_20260704/runs/macosko_sharpen_seed42 \
  --dataset_name Macosko \
  --label_key resolved_label \
  --n_clusters 12 \
  --assignment_mode sharpen \
  --seed 42 \
  --gpu 1
```

Summarize completed runs:

```bash
python experiment_reports/objective_swap_scmae_20260704/summarize.py
```

Current diagnostic result:

| Dataset / model | ARI mean +/- sd | Interpretation |
|---|---:|---|
| Macosko PCA(128)+KMeans known-K | 0.8806 +/- 0.0053 | very strong linear baseline |
| Quake PCA(128)+KMeans known-K | 0.8630 +/- 0.0248 | very strong linear baseline |
| Melanoma PCA(128)+KMeans known-K | 0.6166 +/- 0.0759 | competitive but seed-sensitive baseline |
| Macosko cluster-first sharpen, KMeans on latent | 0.1755 +/- 0.0395 | latent geometry is not KMeans-friendly |
| Macosko cluster-first sharpen, direct prototype | 0.4620 +/- 0.4029 | unstable; one seed collapses to one prototype |
| Macosko cluster-first Sinkhorn, direct prototype | 0.1396 | balanced assignment is a bad assumption here |

First conclusion:

- The first post-scMAE prototype is not yet a viable method.
- Latent std-floor succeeds geometrically (`dims_std>1 = 128/128`) but does not
  make the learned geometry cluster-aligned.
- The unexpectedly strong PCA baseline changes the diagnosis: at least on
  Macosko and Quake, fine-grained structure is already visible in simple linear
  HVG space. Deep reconstruction/prototype objectives can destroy rather than
  improve that structure.
- The next priority is not adding another deep module; it is establishing a
  leakage-safe linear baseline panel and then asking which deep objectives
  preserve or damage that baseline geometry.
