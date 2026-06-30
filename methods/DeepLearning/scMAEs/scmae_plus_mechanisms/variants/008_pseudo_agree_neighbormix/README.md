# 008 Pseudo-Agree NeighborMix

This variant keeps the original scMAE backbone and targets the failure mode
seen after `004` and `007`: NeighborMix can help, but unreliable cross-cluster
edges create high seed variance, while simply softening all edges over-smooths.

Default mechanism:

```text
scMAE + warmup global mutual-KNN + KMeans pseudo-cluster agreement filter
```

At each NeighborMix graph refresh, the current scMAE embedding is clustered
with KMeans using the known `n_clusters`. A mutual KNN edge is eligible only if
both cells currently share the same pseudo cluster. The mix itself stays the
same as `004_neighbormix`:

```text
x'_i = mix_alpha * x_i + (1 - mix_alpha) * x_j
```

No adaptive mask, prototype, SwAV, count-aware, or fuzzy-boundary loss is
enabled by default. The mechanism is deliberately narrow: improve neighbor
selection before adding any extra objective.

## Melanoma_5K Screen Result

Completed on 2026-06-30 with `seeds = 42, 2024, 3407`, `epochs = 80`,
`n_top_genes = 1000`, and `batch_size = 128`.

| setting | ARI mean | ARI std | ACC mean |
|---|---:|---:|---:|
| same pseudo cluster | 0.655011 | 0.018839 | 0.726420 |
| same confident pseudo cluster, q25 | 0.662751 | 0.002924 | 0.736539 |

Neither setting reaches the original scMAE Melanoma_5K ARI reference of about
`0.668`, so this variant is not promoted. The stricter confidence filter does
stabilize the run, but it also removes the high-gain trajectory that made
`004_neighbormix` interesting. This suggests that pseudo labels alone are not a
good enough edge selector; the next refinement should use a teacher/EMA or
consensus graph rather than only the current embedding's KMeans labels.
