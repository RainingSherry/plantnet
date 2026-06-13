# RG Gate Attention Experiment Archive

Archive date: 2026-06-13

This folder packages the code, text outputs, diagnostics, and conclusions for the Reliability-Gated NeighborMix-scMAE experiment described in `BDD/260612_0025门控注意力.md`.

## Scope

The archived experiment covers:

- Phase 1 all-seeds comparison on 8 datasets, 6 RG modes, 3 seeds.
- Phase 2 sensitivity sweep on the 5 BDD representative datasets:
  `SRP182008`, `Melanoma_5K`, `Macosko`, `Tosches`, `Wang`.
- RG method implementation and runner scripts used for the experiment.
- Final BDD conclusion on whether to proceed to Phase 3 and Attention.

Large generated binaries are intentionally excluded from this GitHub archive:

- model checkpoints: `*.pt`
- embeddings and arrays: `*.npy`, `*.npz`
- HDF5 outputs and datasets: `*.h5`, `*.h5ad`

Those files remain in the original local `results/formal/rg_phase2_sensitivity_e80/` run directories.

## Folder Layout

```text
code/
  BDD/
    260612_0025门控注意力.md
  methods/DeepLearning/RG_NeighborMix_scMAE/
    RG implementation, configs, diagnostics, mixing, graph builder, runner
  scripts/
    run_ra_rg_phase.py
    run_rg_phase2_sensitivity.py
    summarize_rg_phase2_sensitivity.py

outputs/
  phase1_allseeds/
    rg_phase1_allseeds_raw.csv
    rg_phase1_allseeds_summary.csv
    rg_phase1_allseeds_group_summary.csv
    rg_phase1_allseeds_reliability_diagnostics.csv
  phase2_sensitivity/
    rg_phase2_all_sweeps_raw.csv
    rg_phase2_all_sweeps_aggregate.csv
    per-run CSV, JSON, and run.log files for the 45 Phase 2 jobs
```

The archive contains 926 files and is about 4.8 MB.

## Phase 1 Key Results

All-seeds global mean ARI:

```text
rg_none        0.620593
rg_fixed       0.619449
rg_reliability 0.625533
rg_random      0.633667
rg_far         0.597840
```

Positive-group mean ARI:

```text
rg_none        0.575297
rg_fixed       0.562365
rg_reliability 0.594522
rg_random      0.592661
```

Negative-group mean ARI:

```text
rg_none        0.483600
rg_fixed       0.471937
rg_reliability 0.455319
rg_random      0.474595
```

Reliability diagnostics indicate that the implementation did not collapse:

```text
mean effective_neighbor_count ~= 9.316
min effective_neighbor_count ~= 8.926
max max_edge_weight_p95 ~= 0.362
max fraction_effective_neighbors_lt_2 ~= 0.0033
mean_node_gate range ~= 0.0811 to 0.0964
max fraction_gate_gt_90pct_max = 0
```

## Phase 2 Key Results

The completed 5-dataset sensitivity sweep contains 45 raw rows:

```text
gate_max:      5 datasets x 4 values x seed42 = 20
neighbor_k:    5 datasets x 3 values x seed42 = 15
pseudo_weight: 5 datasets x 2 values x seed42 = 10
```

Best mean ARI by sweep:

```text
gate_max=0.10       mean ARI = 0.626381
neighbor_k=10       mean ARI = 0.627022
pseudo_weight=0.30  mean ARI = 0.627022
neighbor_k=15       mean ARI = 0.626640
```

## Conclusion

Do not proceed to Phase 3 full 18-dataset runs.

Do not proceed to the Attention version yet.

The BDD pause conditions are triggered:

1. `rg_reliability` does not reduce negative-group damage; it is worse than `rg_fixed` on the negative group.
2. `rg_random` has higher global mean ARI than `rg_reliability`, so the reliability-neighborhood mechanism is not proven.
3. `rg_far` decreases globally, so the far-neighbor negative control behaves as expected, but this is insufficient to validate reliability gating.
4. Gate and edge diagnostics show stable implementation behavior, but stability is not enough to prove mechanism efficacy.

The next step should be mechanism redesign: improve the reliability/gate formulation or add a stricter random-matched control before increasing experiment scale.

## Rebuild Summary Tables

From the repository root:

```bash
python scripts/summarize_rg_phase2_sensitivity.py
```

This regenerates:

```text
results/formal/rg_phase2_sensitivity_e80/rg_phase2_all_sweeps_raw.csv
results/formal/rg_phase2_sensitivity_e80/rg_phase2_all_sweeps_aggregate.csv
```
