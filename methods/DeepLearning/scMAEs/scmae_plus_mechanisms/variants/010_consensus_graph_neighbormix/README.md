# 010 Consensus Graph NeighborMix

This variant keeps the original scMAE backbone and tests the next NeighborMix
hypothesis after the EMA graph result: one smoothed trajectory is not enough,
but edges that persist across several recent graphs may preserve useful dynamic
refinement while rejecting short-lived noisy neighbors.

Default mechanism:

```text
scMAE + warmup global mutual-KNN + multi-update consensus edge filter
```

At each graph refresh, a normal mutual-KNN NeighborMix graph is built from the
current embedding. A candidate edge is eligible for mixing only if it appeared
as an eligible edge in at least `neighbor_consensus_min_hits` of the last
`neighbor_consensus_window` refreshed graphs.

Default:

```text
neighbor_consensus_window = 3
neighbor_consensus_min_hits = 2
```

No adaptive mask, prototype, SwAV, count-aware, or fuzzy-boundary loss is
enabled by default.

## Melanoma_5K Screen Result

Completed on 2026-06-30 with `seeds = 42, 2024, 3407`, `epochs = 80`,
`n_top_genes = 1000`, and `batch_size = 128`.

| setting | ARI mean | ARI std | ACC mean |
|---|---:|---:|---:|
| `window=3, min_hits=2` | 0.696944 | 0.056877 | 0.771918 |
| `window=3, min_hits=3` | 0.666194 | 0.002646 | 0.738090 |
| `window=3, min_hits=3, mix_alpha=0.85` | 0.662422 | 0.001301 | 0.736539 |

The loose 2-of-3 consensus preserves the high mean seen in `004_neighbormix`,
but it remains seed-unstable. The strict 3-of-3 consensus is stable and close to
the original scMAE Melanoma_5K ARI reference of about `0.668`, but it does not
cross the quick-screen threshold. Increasing mix strength hurts the stable
setting. This variant is therefore not promoted yet, but it identifies the next
search direction: confidence-adaptive consensus, with stricter rules for
boundary cells and more permissive rules for high-confidence core cells.
