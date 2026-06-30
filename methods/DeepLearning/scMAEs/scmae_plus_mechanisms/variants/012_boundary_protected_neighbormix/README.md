# 012 Boundary-Protected NeighborMix

This variant keeps the original scMAE backbone and starts from the best stable
anchor found so far:

```text
004_neighbormix with neighbor_start_epoch = 30
```

The added mechanism is a lightweight boundary and rare-cell protection rule for
NeighborMix edges. At each graph refresh, the current scMAE embedding is
clustered with KMeans. Cells in the lowest confidence quantile are treated as
boundary cells, and cells assigned to the smallest pseudo-clusters are treated
as rare-cluster cells. A protected edge is kept only when it connects cells with
the same pseudo label and has sufficient graph reliability.

Default:

```text
neighbor_start_epoch = 30
neighbor_boundary_confidence_quantile = 0.20
neighbor_boundary_rare_quantile = 0.25
neighbor_boundary_score_threshold = 0.84
neighbor_mix_mode = first
mix_alpha = 0.90
```

This does not add a prototype loss, SwAV loss, adaptive mask, or count-aware
branch. The ablation is focused on whether protecting boundary and rare cells
can keep the stability of delayed NeighborMix while recovering the small
remaining ARI gap to the original scMAE baseline.
