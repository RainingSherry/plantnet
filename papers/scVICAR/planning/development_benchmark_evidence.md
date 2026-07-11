# Development benchmark evidence boundary

## What the development table establishes

The archived `result/paper_table_overall.csv` includes the invalid
`hrvatin_geo` label variant and is therefore not a formal paper table.  The
paper-facing table at `tables/development_valid15/paper_table_overall.csv`
excludes that variant and contains only methods with 45/45 successful,
non-fallback runs (15 valid datasets by three seeds).  Its descriptive means
are:

| Method | ARI | NMI | ACC | Macro-F1 |
|---|---:|---:|---:|---:|
| scVICAR-T (RG-NeighborMix-scMAE) | 0.7955 | 0.8374 | 0.8401 | 0.7541 |
| scVICAR-F (NeighborMix-scMAE) | 0.7868 | 0.8346 | 0.8353 | 0.7537 |
| scMAE | 0.7570 | 0.8224 | 0.8134 | 0.7437 |

Relative to scMAE, scVICAR-T has a 0.0385 absolute ARI gain
(5.08% relative), and scVICAR-F has a 0.0298 absolute gain (3.93%
relative).  These values support the claim that the developed configurations
were competitive across a broad benchmark.

The complete comparison contains 17 method configurations and 765 successful,
non-fallback runs. scVICAR-T ranks first and scVICAR-F second on mean ACC, NMI,
ARI, and macro-F1. At the dataset level, scVICAR-T exceeds scMAE on 10 of 15
datasets and ranks in the top three on 11; scVICAR-F exceeds scMAE on 9.

The wider lightweight archive contains 1,055 records across 22 implemented
external or backbone models: 749 successful, 168 failed, 97 fallback, and 41
pending. Eight implementations with incomplete valid-dataset coverage or
fallback outputs are listed in the supplement and excluded from the ranked
complete-case table.

## Why this is development rather than confirmation

- `objective_definition.json` labels the search `benchmark_tuned: true` and
  maximizes mean ARI over selected datasets and seeds subject to penalties.
- `selection_report_full16.md` selects the final parameter signature using the
  16-dataset archive from which the valid 15-dataset table is derived.
- Consequently, the reported full-16 score is a resubstitution/development
  estimate.  It must not be described as an untouched external test estimate.
- The methods-benchmark rows use their protocol configurations, but the current
  evidence does not establish that every external method received a search
  space and compute budget matched to the NeighborMix/RG search.  The manuscript
  must not claim that all methods attained their unknown global optimum.

## Permitted manuscript use

- Describe the table as a broad development-stage benchmark with three seeds
  and complete non-fallback coverage.
- Report absolute and relative effects, while avoiding formal generalization
  claims or p-values that treat the 48 runs as independent biological repeats.
- Put the complete 17-configuration aggregate table in the main paper and the
  dataset-level values, search history, and incomplete-attempt registry in the
  supplement.
- Use the six frozen confirmatory datasets for the primary external-validity
  claim and the separately frozen baseline matrix for the fair protocol-level
  external comparison.

## Prohibited wording

- “Independent validation on 16 datasets.”
- “All competitors were evaluated at their globally optimal parameters.”
- “scVICAR universally improves scMAE.”
- Any significance test treating datasets, seeds, or cells as interchangeable
  independent replicates.
