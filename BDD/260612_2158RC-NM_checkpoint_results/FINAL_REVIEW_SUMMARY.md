# RC-NM v4.1 Feasibility Checkpoint Summary

Execution root:

`results/rc_nm_checkpoint_v4_1`

Formal run count:

- Stage 0: 2 runs
- Stage 1A: 16 runs
- Stage 1B: 50 runs
- Smoke run: 1 run, excluded from formal decision

Frozen Main-5 datasets:

- Wang
- worm_neuron_cell
- Macosko
- Tosches
- Pollen

Stage 1A pseudo objective decision:

- `analytic_RC_rec_only` is better than `analytic_RC_full_pseudo` on both sentinel datasets.
- Stage 1B therefore used `rec_only`.

Stage 1B decision:

- Grade: `B_continue_with_narrowed_claim`
- Mean delta ARI vs `fixed_control`: `+0.001561`
- Mean delta ARI vs `random_delta_matched`: `+0.007101`
- `analytic_RC` beats `random_delta_matched` on 2/5 dataset means.
- `analytic_RC` is positive vs `fixed_control` on 5/10 dataset-seed pairs.

Strict interpretation:

RC-NM remains a feasibility checkpoint. The evidence does not support claiming a mature method, SOTA behavior, top-tier readiness, or validated backbone agnosticism.

Main weaknesses:

- The effect over `fixed_control` is too small for a strong method claim.
- `random_delta_matched` is close to `analytic_RC`; this weakens the reliability-mechanism interpretation.
- Random perturbation matching passes BDD thresholds in only 4/10 Stage 1B random-control runs.
- Same-pipeline Main-5 `none` baselines were not run by the frozen matrix, so negative transfer vs noMix and worst-case drop vs noMix are not directly evaluable in this checkpoint.

Allowed claim:

Reliability-controlled local shrinkage is plausible as a narrow perturbation-controlled regularizer in the scMAE NeighborMix setting, but the reliability mechanism remains weak and needs stronger random matching and cleaner mechanism separation before paper-level claims.

Key tables:

- `rc_nm_checkpoint_all_runs.csv`
- `stage1A_pseudo_objective_sentinel.csv`
- `stage1B_dataset_method_summary.csv`
- `stage1B_pairwise_vs_fixed_control.csv`
- `stage1B_analytic_RC_decision_by_dataset.csv`
- `stage1B_ABC_decision.csv`
