# scVICAR implementation status

Updated: 2026-07-11 14:50 CST.

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
- External baseline matrix complete: 108/108 active identities accepted (18
  runs for each of PCA+KMeans, original scMAE, scVI, scDCC, repaired
  scDeepCluster, and scRCL); 19 superseded runs rejected by the collector.
- Graph-contamination/estimator matrix complete: 126/126 runs; the strict
  runtime audit accepted all 126 with zero rejected or missing runs.
- Recoverable atomic remote upload, one-dataset cache, and 5-GiB disk floor.
- Downstream marker recovery, marker-overlap annotation, and frozen low-label
  probe complete: 108/108 run identifiers, 540 marker splits, and 1,080 probe
  splits passed checksum and run-ID validation.
- Dataset-level statistics, Holm correction, and no-best-seed enforcement.
- Figures 1--6 complete in SVG, PDF, PNG, and TIFF.  Figure 5 includes all
  dataset units, complete source-data CSVs, and a verified SHA-256 manifest.
- TCBB-style manuscript and supplement compile without missing references,
  undefined cross-references, or overfull boxes.
- The manuscript argument has been rebuilt around a five-level evidence chain:
  development breadth, frozen external competitiveness, matched attribution,
  sensitivity analysis, and controlled graph contamination. The title now
  foregrounds bounded graph-vicinal anchor recovery, while topology adaptation
  is framed as conditional risk control.
- The complete-label 54-run sensitivity analysis now appears in the main
  Results section, and downstream predictive utility is explicitly separated
  from marker coherence.
- A section-matched Chinese reading draft is maintained under
  `manuscript_zh/` for author review. It preserves formulas, quantitative
  results, statistical qualifiers, and claim boundaries from the English
  submission source.
- The English main text, supplement, generated prose, and Chinese reading draft
  passed a direct-writing audit. Template contrasts (`not X, but Y`), repeated
  self-defense, and AI-like meta-narrative were removed. The author-defined
  `cell-wise perturbation budget` was renamed to the concrete
  `cell-specific mixing coefficient`; terminology provenance is archived under
  `planning/terminology_provenance_audit.md`.
- A non-speculative Data and Code Availability section now records public data
  provenance, the no-H5AD-redistribution boundary, checksum traceability, and
  the planned archival code release without inventing a repository URL.
- Local protocol suite: 30 tests passed, including independent Leiden metric
  recomputation, dataset-level hierarchical summary checks, and runtime
  resolution of public path placeholders.
- PCA baseline end-to-end contract and scVI one-epoch artifact contract verified.
- About 34 GiB of historical run products and regenerable staging files were
  moved to the matching `/data` tree and exposed through verified symbolic
  links; `/home` free space increased from about 52 to 86 GiB.
- The development benchmark is now a main-paper result: 17 complete method
  configurations, 15 valid datasets, three seeds, and 765 successful
  non-fallback runs. scVICAR-T and scVICAR-F rank first and second on mean ACC,
  NMI, ARI, and macro-F1. Their ARI gains over scMAE are +5.08% and +3.93%.
  scVICAR-T exceeds scMAE on 10/15 datasets and ranks in the top three on 11.
- Dataset-level benchmark values, paired effects, and the wider 1,055-record,
  22-model attempt registry are archived in the supplement. Eight incomplete
  or fallback implementations remain visible with status counts. The invalid
  `hrvatin_geo` label variant is excluded from formal ranking, and the final
  F/T configurations remain identified as development-selected.
- The benchmark methods now carry source citations for Louvain, Leiden, SC3,
  DEC, scNAME, scziDesk, scCDCG, PhytoCluster, and the previously cited deep
  baselines. The English paper, Chinese reading draft, evidence-boundary note,
  and generated tables use the same benchmark counts and claims.

## Final evidence status

- The stress aggregate and Figure 4 are complete.  Topology affinity detects
  injected cross-class edges (mean AUROC 0.918) and the gate tracks weighted
  purity (mean Spearman 0.520), while the relative ARI response remains
  heterogeneous (T safer in 2/3 datasets at 100% contamination).
- The separate unfiltered full-label NoMix/F/T sensitivity matrix is complete:
  54/54 runs accepted, zero rejected, and all clustering metrics independently
  recomputed.  scVICAR-T changed mean ARI versus NoMix by +0.030 across six
  datasets (5/0/1 wins/ties/losses; 95% CI -0.009 to +0.081; Holm p=0.5077).
  This sensitivity evidence is never pooled with the primary matrix.
- scDeepCluster's original adapter exposed negative Z-scored values to a ZINB
  count target and allowed a failed legacy HVG selection to retain 61,497
  genes, producing infinite loss.  The separately frozen repair_v2 aligns
  non-negative raw counts, recomputes raw-library size factors, and enforces a
  deterministic 2,000-HVG bound.  A smoke test and all 18 formal repaired runs
  have finite metrics and match the repair identity; old infinite-loss
  directories are excluded by the collector.
- The strict release audit passes every local evidence gate: primary, external
  baselines, stress, stress CUDA audit, fixed Leiden, downstream, full-label
  sensitivity, Figures 1--6, marker-panel verification, and data-license audit.
- The remote release audit passes all 24 local/remote checks. The polished
  content-addressed paper snapshot `2de2d93228f5bc30` (88 selected source
  artifacts, including the Chinese reading draft and terminology audits, plus
  metadata and checksum/completion files) is atomically archived under remote
  `paper_snapshots/`. Earlier complete snapshots remain immutable.
- The benchmark-integrated local paper snapshot is `ed5c4c58ef1664c8` with 94
  selected artifacts. It was staged under `/tmp` because the configured staging
  symlink points to a read-only `/data` target in the current session. Remote
  upload is pending an explicitly supplied `SCVICAR_SSH_KEY`; the earlier
  remote snapshot remains complete and immutable.
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
