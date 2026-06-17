# Paper Claim Tracker

| Claim | Evidence | Status |
| --- | --- | --- |
| random_beta_uniform_0.1 is not just lower mean beta | Stage 1: +0.0638 ARI vs fixed_beta_0.05 at same mean beta, 12/18 run wins | Supported |
| Low-variance stochastic beta is stronger than uniform beta | Stage 2: truncated_normal mean ARI 0.6606 vs uniform 0.6562 and fixed 0.5925 | Supported |
| Bernoulli mix/no-mix switching is the mechanism | Stage 2: Bernoulli mean ARI 0.5854, below fixed_beta_0.05 | Not supported |
| Local neighborhood structure is necessary | Stage 3: local anchor mix mean ARI 0.6606 vs global random 0.6170; full benchmark global random was last | Supported |
| Anchor-recovery is better than mixed-target expansion | Stage 3: anchor local mean ARI 0.6606 vs mixed target 0.6482 | Supported |
| Benefit is not generic matched noise | Stage 3: local anchor had better macro-F1/rank/worst-case than Gaussian, but Gaussian mean ARI was nearly equal | Partially supported |
| random beta reduces wrong-edge damage | Stage 4: at 40% bad edges, random_beta_uniform_0.1 mean ARI 0.4872 vs fixed_beta_0.1 0.4360 | Partially supported |
| stochastic beta generalizes on full benchmark | Full benchmark: truncated-normal and uniform ranked first/second by mean ARI and macro-F1 | Supported on average |
| stochastic beta is robustly better in worst case | Full benchmark: worm_neuron_cell worsens vs noMix; worst-case delta is worse than fixed_beta_0.1 | Not supported |

Decision:

```text
Proceed with stochastic beta NeighborMix only as an average-performance and mechanism-driven main line.
Do not claim robust improvement without an additional safeguard for worm_neuron_cell or boundary/high-risk cells.
```
