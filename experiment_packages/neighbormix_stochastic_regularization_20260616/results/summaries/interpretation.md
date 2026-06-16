Experiment context:
Negative-transfer group: Tosches, Macosko, worm_neuron_cell. These datasets test whether a method reduces known NeighborMix downside risk.
Positive-gain group: Melanoma_5K, Shekhar. These datasets test whether a method preserves cases where NeighborMix can help.
Neutral/stable group: Guo. This dataset tests whether a method avoids unnecessary degradation when the baseline is already stable.
Optional validation group: Wang, Pollen. These were run after the six core datasets finished, but the continuation criteria are judged on the six core datasets.

Mechanism framing:
The pseudo branch reconstructs anchor cells from mixed pseudo-cell inputs. The random variants therefore test stochastic neighborhood regularization, not reliable-cell discovery.
The alternative-neighbor variants test whether the vanilla PCA-cosine KNN graph is the bottleneck behind negative transfer.
The global-random control tests whether any gains survive after removing local neighborhood structure.

Finding 1:
Best random variant: random_beta_uniform_0.1. It beats fixed NeighborMix on 5/6 core datasets.

Evidence:
Random-variant screen:
- random_pseudo_gate_p0.5: beats fixed on 2/6 core datasets; negative mean delta ARI -0.0041; worst delta ARI -0.1731; mean macro-F1 0.7311.
- random_edge_dropout_keep0.5: beats fixed on 0/6 core datasets; negative mean delta ARI 0.0030; worst delta ARI -0.1495; mean macro-F1 0.7360.
- random_beta_uniform_0.1: beats fixed on 5/6 core datasets; negative mean delta ARI 0.0140; worst delta ARI -0.1183; mean macro-F1 0.7562.
Negative-group mean delta ARI: best random 0.0140, fixed NeighborMix 0.0129.
Worst-case delta ARI: best random -0.1183, fixed NeighborMix -0.0735.
Mean macro-F1: best random 0.7562, fixed NeighborMix 0.7476.

Interpretation:
The stochastic-neighborhood direction does not pass the continuation screen. Any isolated gains are not sufficient to claim a robust regularization benefit.

Next action:
Continue only if the detailed per-dataset table shows the gains are not concentrated in one dataset; otherwise treat this as a negative result.

Finding 2:
Best alternative-neighbor variant: snn_neighbormix.

Evidence:
Alternative-neighbor screen:
- mutual_knn_neighbormix: negative mean delta ARI -0.0044; worst delta ARI -0.1757; mean macro-F1 0.7503.
- snn_neighbormix: negative mean delta ARI 0.0125; worst delta ARI -0.1654; mean macro-F1 0.7416.
- consensus_neighbormix_threshold0.4: negative mean delta ARI -0.0283; worst delta ARI -0.1507; mean macro-F1 0.7291.
Negative-group mean delta ARI: best alternative 0.0125, fixed NeighborMix 0.0129.
Weighted same-label edge ratio: best alternative 0.9240, fixed NeighborMix 0.9226.
Perturbation norm mean: best alternative 2.6798, fixed NeighborMix 2.5364.

Interpretation:
The alternative-neighbor direction does not pass the continuation screen. Either the graph did not improve edge purity enough, or the apparent gains came from weakening the perturbation.

Next action:
Proceed to edge-level gates or MoE only if the accepted neighbor rule also improves macro-F1 or minority diagnostics; otherwise improve graph construction first.

Finding 3:

Evidence:
Global random mean delta ARI: -0.0084; fixed NeighborMix mean delta ARI: 0.0163; best random negative-group delta ARI: 0.0140.

Interpretation:
Global random neighbors are worse than local stochastic variants, so any useful regularization still depends on local neighborhood structure.

Next action:
If global random is competitive, do not build complex gates yet; first test whether the gain is generic noise regularization.

Overall next action:
Given the current screen, do not escalate directly to complex attention or MoE. The most defensible next step is a small edge-level gate only after improving graph diagnostics, because the tested mutual/SNN/consensus rules did not clearly reduce negative transfer.
