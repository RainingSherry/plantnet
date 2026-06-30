# 007 Soft Reliable NeighborMix

This variant keeps the original scMAE backbone and focuses only on stabilizing
the useful NeighborMix signal observed in `004_neighbormix`.

Default mechanism:

```text
scMAE + warmup global mutual-KNN + reliability-weighted first-neighbor mix
```

The key difference from `004_neighbormix` is that the first reliable neighbor is
not mixed at a fixed strength. Its contribution is scaled by the edge
reliability score:

```text
beta_ij = (1 - mix_alpha) * reliability_ij ** neighbor_soft_power
x'_i = (1 - beta_ij) x_i + beta_ij x_j
```

This is intended to be a middle ground between fixed first-neighbor mixing,
which had the best mean ARI but high seed variance, and mean/weighted-mean
mixing over all mutual neighbors, which was stable but over-smoothed the
embedding. No adaptive mask, prototype, SwAV, count-aware, or fuzzy-boundary
loss is enabled by default.

## Melanoma_5K Screen Result

Completed on 2026-06-30 with `seeds = 42, 2024, 3407`, `epochs = 80`,
`n_top_genes = 1000`, and `batch_size = 128`.

| setting | ARI mean | ARI std | ACC mean |
|---|---:|---:|---:|
| default `mix_alpha=0.75` | 0.663174 | 0.004279 | 0.736982 |
| `mix_alpha=0.88` | 0.663165 | 0.002333 | 0.737056 |

Both settings are below the original scMAE Melanoma_5K ARI reference of about
`0.668`, so this variant is not promoted. The negative result is still useful:
soft edge weighting stabilizes NeighborMix, but it also removes the high-ARI
positive seed seen in fixed-strength first-neighbor mixing. Future refinements
need better edge selection, not only softer edge strength.
