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

## Melanoma_5K Screen Result

Completed on 2026-06-30 with `seeds = 42, 2024, 3407`, `epochs = 80`,
`n_top_genes = 1000`, and `batch_size = 128`.

| setting | ARI mean | ARI std | ACC mean |
|---|---:|---:|---:|
| default boundary + rare protection | 0.661512 | 0.001316 | 0.736613 |
| boundary only, rare disabled | 0.662414 | 0.002045 | 0.736391 |
| boundary only, q10 | 0.660343 | 0.003630 | 0.735653 |

The protected graph is very clean: reliable edges are almost all within the
same pseudo cluster, and variance stays low. However, the rule is too
conservative and removes useful NeighborMix signal. This variant is therefore
not promoted. The next boundary-aware attempt should use soft weighting or an
auxiliary consistency loss for boundary cells instead of hard edge deletion.
