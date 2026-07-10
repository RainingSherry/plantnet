# scVICAR paper and reproducibility package

This directory is the paper-facing implementation of **scVICAR** (single-cell
Vicinal Corruption with Anchor Recovery). It deliberately leaves the retired
NeighborMix lineage as its implementation starting point and exposes the RG
runner through one matched-backbone family. Stress-test hooks are isolated by
zero-valued defaults, and frozen-source tests protect the clean execution path.

The current execution/acceptance state is tracked in `PROGRESS.md`; missing
formal results are represented as blocked generated assets rather than interim
numbers copied into the manuscript.

## Evidence boundary

- The historical 16-dataset table is development evidence, not an independent test.
- `hrvatin_geo` is excluded from formal evidence because missing annotations were
  encoded as an additional class and it is not independent of `hrvatin`.
- Confirmatory claims use six datasets that were not used to choose NeighborMix or
  RG-NeighborMix parameters.
- `scVICAR-T` is an adaptive extension, not a claim of universal superiority over
  `scVICAR-F`.

## Local and remote layout

Large datasets and run artifacts live in the authorized remote stores:

```text
<SCVICAR_DATA_ROOT>
<SCVICAR_PROJECT_ROOT>
```

The private-key path is supplied only through `SCVICAR_SSH_KEY`; it is never
written to a run record.

## Reproduction

Run commands from the repository root.

```bash
export SCVICAR_SSH_KEY=/path/to/authorized/private_key
export MPLCONFIGDIR=/tmp/scvicar-matplotlib
export NUMBA_CACHE_DIR=/tmp/scvicar-numba

# Validate/create remote layout and build canonical datasets.
python -m papers.scVICAR.code.remote_store ensure-layout
python -m papers.scVICAR.code.prepare_datasets --upload \
  --output-dir papers/scVICAR/manifests/dataset_upload

# Inspect the frozen matrix without running it.
python -m papers.scVICAR.code.orchestrate --dry-run

# Six-variant two-epoch smoke test on one dataset/seed.
python -m papers.scVICAR.code.orchestrate \
  --datasets Mouse_Pancreas_1 --seeds 42 --epochs 2 --fail-fast

# Formal 108-run confirmatory matrix.
python -m papers.scVICAR.code.orchestrate --fail-fast

# Frozen external baselines (108 runs) and graph-stress matrix (126 runs).
python -m papers.scVICAR.code.freeze_baselines
python -m papers.scVICAR.code.baseline_orchestrate
python -m papers.scVICAR.code.stress_orchestrate

# After the primary matrix is complete, run every downstream split and aggregate.
python -m papers.scVICAR.code.downstream_orchestrate
python -m papers.scVICAR.code.downstream_aggregate

# Fixed-resolution, non-oracle Leiden secondary evaluation and aggregation.
python -m papers.scVICAR.code.secondary_evaluation
python -m papers.scVICAR.code.secondary_aggregate

# Separate unfiltered-label sensitivity layer (54 runs; never pooled with primary).
python -m papers.scVICAR.code.prepare_full_label_sensitivity --upload
python -m papers.scVICAR.code.freeze_full_label_sensitivity
python -m papers.scVICAR.code.full_label_sensitivity_orchestrate --dry-run
python -m papers.scVICAR.code.full_label_sensitivity_orchestrate --fail-fast
python -m papers.scVICAR.code.full_label_sensitivity_aggregate

# Final strict paper build: refuses any incomplete evidence layer.
python -m papers.scVICAR.code.generate_manuscript_assets --require-all
```

The orchestrator downloads one canonical dataset at a time, schedules runs only
on GPUs 1--6, atomically uploads every verified result, and removes local run
artifacts only after the remote checksum passes.

## Project map

```text
code/          model, storage, downstream, aggregation, statistics and figures
configs/       immutable protocol files
experiments/   planned matrix and small status/summary files
figures/       figure contracts, source data and final vector/raster exports
manifests/     canonical dataset and remote artifact indices
manuscript/    TCBB-style LaTeX manuscript and supplement
planning/      terminology, claim-evidence and reviewer-risk ledgers
tables/        generated paper tables
tests/         unit and smoke-test assertions
```

## Safety rules

- Never use rsync `--delete`.
- Never overwrite an immutable remote dataset or completed run.
- Stop scheduling when local free space is below 5 GiB.
- Never select the best random seed.
- Labels may be used for the predefined class filter, known-K evaluation and
  post-hoc biological analysis, but not for training, graph construction or
  hyperparameter selection.
- Per-cell labels are withheld from all optimization loops. External clustering
  baselines whose published architecture requires the number of clusters
  (scDCC, scDeepCluster and scRCL) are explicitly marked as known-K training
  baselines; they are not described as fully label-free model-selection systems.
- Every external baseline must emit exact cell IDs. A run fails unless those IDs
  match the canonical H5AD row-for-row before any metric or `COMPLETED` marker is
  written.
