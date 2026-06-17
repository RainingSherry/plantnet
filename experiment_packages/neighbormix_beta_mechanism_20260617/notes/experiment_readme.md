# Experiment Readme

Objective:

```text
Disentangle whether random_beta_uniform_0.1 works because mean beta is lower or because beta is stochastic.
```

Result:

```text
Stochastic beta has real average-performance support.
Controlled-variance stochastic beta is stronger than the original uniform beta.
The robust-improvement claim is not yet safe because worm_neuron_cell remains a failure case.
```

Recommended paper framing:

```text
Stochastic perturbation-strength regularization for NeighborMix-scMAE.
```

Avoid this stronger claim:

```text
Stochastic beta NeighborMix is robustly better on all datasets.
```

Most important tables:

```text
results/summaries/stage1_beta_mean_vs_randomness.csv
results/summaries/stage2_beta_variance.csv
results/summaries/stage3_local_mix_mechanism.csv
results/summaries/stage4_bad_edge_robustness.csv
results/summaries/full_benchmark_summary.csv
results/summaries/interpretation.md
```
