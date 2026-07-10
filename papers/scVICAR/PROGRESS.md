# scVICAR implementation status

Updated: 2026-07-10 20:16 CST.

## Complete

- Six canonical confirmatory H5AD files filtered, uploaded, and SHA-256 verified.
- Active primary freeze: 3ddccaae58a3a468.
- Active external-baseline freeze: b98f47e8da895988.
- Active stress-test freeze: bac5214258202a09.
- Active fixed-Leiden freeze: 9f6c9e67b4c221e0.
- Active unfiltered full-label sensitivity freeze: 0a13b1aad8fa97f5.
- Primary confirmatory matrix complete: 108/108 active-freeze runs, zero failures.
- Primary aggregation, dataset-level contrasts, and Figures 2--3 complete.
- Fixed-resolution non-oracle Leiden complete: 108/108 runs independently
  recomputed and checksum-verified; 36 dataset-variant summaries and 12 planned
  contrasts are integrated into the compiled manuscript.
- PCA baseline complete: 18/18 active-freeze runs; 19 superseded runs rejected
  by the active-freeze collector.
- External baseline matrix planned: 108 runs.
- Graph-contamination/estimator matrix planned: 126 runs.
- Recoverable atomic remote upload, one-dataset cache, and 5-GiB disk floor.
- Downstream marker recovery, marker-overlap annotation, and frozen low-label probe.
- Dataset-level statistics, Holm correction, and no-best-seed enforcement.
- Figures 1--3 and the preregistered Human_Pancreas_3 Figure 6 complete in SVG,
  PDF, PNG, and TIFF; Figures 4--5 remain fail-closed on complete aggregates.
- TCBB-style manuscript and supplement compile without missing references,
  undefined cross-references, or overfull boxes.
- A non-speculative Data and Code Availability section now records public data
  provenance, the no-H5AD-redistribution boundary, checksum traceability, and
  the planned archival code release without inventing a repository URL.
- Local protocol suite: 27 tests passed, including independent Leiden metric
  recomputation and dataset-level hierarchical summary checks.
- PCA baseline end-to-end contract and scVI one-epoch artifact contract verified.

## Running/pending evidence

- The 126-run graph-contamination/estimator stress matrix is running with
  checksum-verified resume and one worker per GPU (109/126 remote completion
  markers at 20:16 CST).
- The 108-run downstream matrix is running from the frozen primary embeddings
  (74/108 remote completion markers at 20:16 CST).
- Six unfiltered full-label H5AD files are uploaded and SHA-256 verified.  The
  separate 54-run NoMix/F/T sensitivity matrix is frozen but not yet launched;
  it is never pooled with the primary matrix.
- Remaining GPU baselines are deferred until the stress matrix releases the
  devices; this avoids cross-protocol resource contention.
- Figures 4--5 and their quantitative prose remain deliberately blocked until
  the corresponding complete matrices pass checksum and run-ID validation;
  Figure 6 is complete and visually audited.
- The Human_Pancreas_3 marker panel has been restricted to markers directly
  verified against the cited pancreas reference and remains interpretation-only.

## Formal launch command

    SCVICAR_SSH_KEY=<PATH_TO_PRIVATE_SSH_KEY> \
    MPLCONFIGDIR=/tmp/scvicar-matplotlib \
    NUMBA_CACHE_DIR=/tmp/scvicar-numba \
    OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
    PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True \
    PYTHONUNBUFFERED=1 \
    python -m papers.scVICAR.code.orchestrate \
      --datasets Blood_BoneMarrow Human_Pancreas_1 Human_Pancreas_3 \
                 Mouse_Pancreas_1 PRJNA895163 TabulaSapiens_Pancreas \
      --variants nomix random_mix fixed topology_edge_only \
                 topology_gate_only topology_full \
      --seeds 42 2024 3407 --gpus 2 3 4 5 6

The command must run outside the filesystem/GPU sandbox. Completed active-hash
runs are checksum-verified and skipped; any failure preserves staging and pauses
before the next dataset.
