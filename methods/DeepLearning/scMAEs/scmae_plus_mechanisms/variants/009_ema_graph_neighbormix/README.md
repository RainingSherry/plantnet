# 009 EMA Graph NeighborMix

This variant keeps the original scMAE backbone and targets the most recent
NeighborMix diagnosis: fixed graphs are stable but lose the gain, while fully
dynamic graphs can amplify noisy trajectories.

Default mechanism:

```text
scMAE + warmup global mutual-KNN + EMA-smoothed embedding graph
```

The training loss is the same as `004_neighbormix`. The only change is how the
global KNN graph is built at refresh epochs:

```text
z_ema(t) = decay * z_ema(t-1) + (1 - decay) * z_current(t)
graph = KNN(z_ema(t))
```

This keeps beneficial graph refinement over training while making each graph
update less sensitive to the current epoch's embedding noise. No adaptive mask,
prototype, SwAV, count-aware, or fuzzy-boundary loss is enabled by default.

## Melanoma_5K Screen Result

Completed on 2026-06-30 with `seeds = 42, 2024, 3407`, `epochs = 80`,
`n_top_genes = 1000`, and `batch_size = 128`.

| setting | ARI mean | ARI std | ACC mean |
|---|---:|---:|---:|
| `neighbor_embedding_ema_decay=0.8` | 0.663759 | 0.003946 | 0.737056 |
| `neighbor_embedding_ema_decay=0.2` | 0.661184 | 0.009852 | 0.730999 |

Neither setting reaches the original scMAE Melanoma_5K ARI reference of about
`0.668`, so this variant is not promoted. The result shows that smoothing a
single dynamic graph trajectory stabilizes NeighborMix but still removes the
large positive seed seen in `004_neighbormix`. A stronger next step is a
multi-update consensus graph or separate teacher graph rather than simple EMA
of the current student embedding.
