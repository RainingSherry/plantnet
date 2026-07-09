# Safe-RDG-PCA V3 Branch/Stability Progress

## Purpose

V3 changes the paper route from a label-calibrated main method to a dual-track
design:

- `Safe-RDG-PCA-U`: strict unsupervised branch/stability gate for the main
  known-K benchmark.
- `Safe-RDG-PCA-C`: label-calibrated meta-selector, kept only as an appendix or
  calibrated selector analysis.

The main claim is not that Safe-RDG-PCA is a universal clustering SOTA. The
claim is:

```text
Safe-RDG-PCA studies negative-transfer-aware selective graph augmentation over a strong PCA baseline.
```

## V3-U Gate

The V3-U variant is named:

```text
safe_rdg_pca_u
```

It uses only label-free evidence from the current dataset:

- KMeans assignment stability.
- Neighbor stability.
- Graph or embedding perturbation stability.
- Cluster-size sanity.
- Hubness penalty.
- Eigengap score.
- Negative-control contrast against random, degree-shuffled, and shuffled
  gene-cell graphs.

No ARI, labels, or LODO-learned thresholds are used by V3-U.

## Output Additions

Every variant now writes:

- `stability_diagnostics.json`
- `selector_scores.json`

For non-selector variants these files can be empty. For `safe_rdg_pca_u`, they
record all candidate branch scores, negative-control scores, selected branch,
fallback reason, and thresholds.

## Development Protocol

Use the existing benchmark-compatible v2 data source and environment:

```text
result/scmae_all_methods_20260705_full/converted_data
/data/luolie/conda/envs/scclubench-main/bin/python
```

First run smoke on:

```text
Bach,Pollen,Wang,worm_neuron_cell
```

Then run 16 datasets x 3 seeds development if smoke passes. Do not use 23
formal until V3-U thresholds and weights are frozen.

## Strict Boundary

`safe_rdg_calibrated_threshold` remains useful but is not a strict unsupervised
method. It should not enter the main unsupervised benchmark table. It belongs in
a separate cross-dataset calibrated selector analysis.
