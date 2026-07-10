# Full-label sensitivity protocol

Status: prespecified before execution; separate from the 108-run confirmatory
matrix.

## Purpose

The primary protocol removes normalized labels `unknown`, `nan`, and
`unassigned`, then removes classes with fewer than ten cells. This sensitivity
analysis asks whether the conclusions for the three principal matched models
depend on that evaluation-cohort filter.

## Data construction

- Start from the same six immutable source H5AD files and source label fields.
- Retain every source cell; do not remove missing/unknown/unassigned or small
  classes.
- Convert the source label to a string-valued `resolved_label` without merging
  categories. A missing value therefore remains the explicit string `nan`.
- Record source path, dimensions, full label distribution, SHA-256, and creation
  time in a separate `sensitivity_full_labels_v1` manifest.
- Store immutable data under
  `<SCVICAR_DATA_ROOT>/datasets/sensitivity_full_labels_v1/`.

## Frozen comparison

- Variants: NoMix, scVICAR-F, and scVICAR-T.
- Model seeds: 42, 2024, and 3407.
- Total: 6 datasets × 3 variants × 3 seeds = 54 runs.
- The backbone, preprocessing, optimizer, epoch budget, and variant parameters
  are identical to protocol_v1.
- The full observed category count supplies K only to post-hoc KMeans; labels do
  not enter optimization, clean graph construction, checkpoint selection, or
  scheduling.

## Analysis and reporting

- Report ARI, NMI, ACC, and macro-F1 after averaging model seeds within dataset.
- Use the same three paired comparisons among NoMix, F, and T, with dataset as
  the independent unit, paired hierarchical bootstrap, paired sign-flip test,
  and within-endpoint Holm correction.
- Present the results in the supplement as sensitivity evidence only.
- Do not combine these runs with the primary 108-run matrix, do not select a
  favorable filtered/unfiltered result, and retain adverse or null outcomes.
