# 011 Confidence-Adaptive NeighborMix

This variant keeps the original scMAE backbone and targets the failure mode
seen in `004_neighbormix` and `010_consensus_graph_neighbormix`: loose dynamic
NeighborMix can produce a high mean ARI but has large seed variance, while
strict consensus is stable but slightly below the scMAE quick-screen threshold.

Default mechanism:

```text
scMAE + warmup global mutual-KNN + confidence-adaptive consensus NeighborMix
```

At each graph refresh, the latest mutual-KNN graph is compared with recent
graphs. A candidate edge is kept when either:

```text
hit_count >= neighbor_adaptive_strict_hits
```

or:

```text
hit_count >= neighbor_adaptive_loose_hits
and average_reliability >= neighbor_adaptive_score_threshold
```

The default is therefore strict for boundary or low-confidence edges, but still
allows high-confidence core-cell edges with 2-of-3 recent graph agreement:

```text
neighbor_consensus_window = 3
neighbor_adaptive_loose_hits = 2
neighbor_adaptive_strict_hits = 3
neighbor_adaptive_score_threshold = 0.84
neighbor_mix_mode = first
mix_alpha = 0.90
```

The confidence rule is applied to edge selection. The default then uses the
same first-neighbor mixing strength as `004_neighbormix`:

```text
x'_i = mix_alpha * x_i + (1 - mix_alpha) * x_j
```

No adaptive mask, prototype, SwAV, count-aware, or fuzzy-boundary loss is
enabled by default. This keeps the ablation focused on reliable NeighborMix edge
selection only. A `soft_first` ablation is still supported through CLI flags,
but the screen showed that it over-smoothed the embedding.

## Melanoma_5K Screen Result

Completed on 2026-06-30 with `seeds = 42, 2024, 3407`, `epochs = 80`,
`n_top_genes = 1000`, and `batch_size = 128`.

| setting | ARI mean | ARI std | ACC mean |
|---|---:|---:|---:|
| default, `score_threshold=0.84`, `first`, `mix_alpha=0.90` | 0.665941 | 0.001270 | 0.737351 |
| `score_threshold=0.80`, `first`, `mix_alpha=0.90` | 0.665406 | 0.002917 | 0.736982 |
| earlier `score_threshold=0.84`, `first`, `mix_alpha=0.90` run | 0.666288 | 0.001824 | 0.737942 |

The adaptive confidence rule stabilizes NeighborMix but does not cross the
Melanoma_5K quick-screen threshold of `ARI_mean > 0.668`. It behaves like the
other conservative NeighborMix stabilizers: seed variance falls, but the large
positive dynamic feedback from `004_neighbormix` is also removed.
