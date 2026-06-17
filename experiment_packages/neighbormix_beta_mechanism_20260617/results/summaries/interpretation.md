Finding 1:
Stage 1 supports a stochastic beta contribution beyond lower mean beta.

Evidence:
The key same-mean comparison is fixed_beta_0.05 vs random_beta_uniform_0.1. random_beta_uniform_0.1 improved mean ARI by 0.0638 and mean macro-F1 by 0.0257 over fixed_beta_0.05, with 12 wins and 6 losses across 18 dataset-seed pairs. At mean beta 0.025, random_beta_uniform_0.05 also beat fixed_beta_0.025 by 0.0589 ARI. At mean beta 0.1, random_beta_uniform_0.2 was worse than fixed_beta_0.1 by -0.0222 ARI.

Interpretation:
The original random_beta_uniform_0.1 result is not explained only by halving mean beta from 0.1 to 0.05. Stochasticity matters at weak to moderate perturbation strength, but too much beta variance/upper range hurts.

Next action:
Use stochastic beta as a mechanism candidate, but avoid claiming monotonic benefit from randomness.

Finding 2:
Stage 2 favors continuous low-variance stochastic beta over fixed beta and Bernoulli switching.

Evidence:
At mean beta 0.05, truncated_normal_beta_mean0.05_std0.02 had the best Stage 2 mean ARI (0.6606), mean macro-F1 (0.7205), and mean rank (2.28). uniform_beta_0.1 was second (mean ARI 0.6562). fixed_beta_0.05 was much lower (0.5925), and bernoulli_beta_0_or_0.1_p0.5 was not competitive (0.5854).

Interpretation:
The useful signal is not simply "sometimes mix, sometimes do not mix." A continuous stochastic perturbation strength is more stable.

Next action:
If a paper method is selected from this package, the stronger candidate is stochastic beta with controlled variance, not necessarily the original uniform distribution.

Finding 3:
Stage 3 supports local structure and anchor recovery, but ARI alone does not fully rule out generic matched noise.

Evidence:
anchor_target_local_mix had mean ARI 0.6606, mean macro-F1 0.7205, and mean rank 2.00. mixed_target_local_mix was lower in mean ARI (0.6482). global_random_mix_anchor_target was clearly worse (0.6170), supporting local neighborhood structure. gaussian_noise_matched_anchor_target had almost identical mean ARI (0.6606) but worse macro-F1 (0.6941), worse mean rank (2.83), and a much worse worst-case delta ARI (-0.2778 vs -0.1584).

Interpretation:
Local convex mixing and anchor recovery remain useful. However, because Gaussian matched noise is close in average ARI, the claim "benefit is not generic noise" should be framed with macro-F1, rank, and worst-case evidence rather than ARI alone.

Next action:
For a paper, include Gaussian matched noise as an important caveat/control instead of hiding it.

Finding 4:
Stage 4 gives partial evidence that random_beta_uniform_0.1 is more tolerant than fixed_beta_0.1 under severe bad-edge injection.

Evidence:
At 40% injected bad edges on the three failure datasets, random_beta_uniform_0.1 had mean ARI 0.4872, while fixed_beta_0.1 had 0.4360 and fixed_beta_0.05 had 0.4668. The ARI change from bad0 to bad0.4 was +0.0066 for random_beta_uniform_0.1, -0.0435 for fixed_beta_0.1, and +0.0097 for fixed_beta_0.05. At bad0.2, random_beta_uniform_0.1 was not best, so the robustness evidence is not perfectly monotonic.

Interpretation:
Random beta likely reduces damage from very noisy neighborhoods, but this is not a clean monotonic robustness curve.

Next action:
Use bad-edge injection as supporting mechanism evidence, not as the sole basis for a robust-improvement claim.

Finding 5:
The 8-dataset full benchmark supports stochastic beta as an average-performing main line, but not an unconditional robust-improvement claim.

Evidence:
Full benchmark completed 8 datasets x 3 seeds x 8 variants = 192 runs. truncated_normal_beta_mean0.05_std0.02 ranked first by mean ARI (0.7270), mean macro-F1 (0.7571), and mean rank (2.85). random_beta_uniform_0.1 ranked second by mean ARI (0.7234), mean macro-F1 (0.7562), and mean rank (3.44). Both beat fixed_beta_0.1, fixed_beta_0.05, and nm_scmae_nomix in mean ARI and macro-F1. random_beta_uniform_0.1 beat fixed_beta_0.1 by +0.0406 ARI on average and won 16/24 run-level comparisons.

Interpretation:
The first paper can reasonably use stochastic beta NeighborMix as the main mechanism direction, with the stronger empirical variant being controlled-variance stochastic beta. But the claim must be "better average performance and better rank," not "robust improvement on every dataset." worm_neuron_cell remains a clear failure case: random_beta_uniform_0.1 averaged 0.4121 ARI vs noMix 0.5167, and truncated_normal averaged 0.4064.

Next action:
Write the paper claim as stochastic perturbation-strength regularization with explicit failure-case analysis. Do not claim worst-case robustness unless an additional cell-level or dataset-level safeguard fixes worm_neuron_cell.
