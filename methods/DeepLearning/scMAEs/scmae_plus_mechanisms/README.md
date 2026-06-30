# scMAE Plus Mechanisms

This series keeps the original `methods/DeepLearning/scMAE` MLP encoder, mask
predictor, decoder, cross-cell replacement corruption, weighted reconstruction,
mask-prediction loss, `feature(x)` embedding, and KMeans-known-k evaluation.

The goal is no longer to port whole outside architectures. Each variant adds one
local mechanism that should make the final embedding more KMeans-friendly.

## First-stage variants

- `001_adaptive_mask`: scMAE plus gene-statistic adaptive masking.
- `002_prototype_dec`: scMAE plus delayed DEC-style prototype alignment.
- `003_swav_assignment`: scMAE plus two-view SwAV-style swapped assignment.
- `004_neighbormix`: scMAE plus warmup global mutual-KNN NeighborMix.
- `005_neighbormix_prototype`: NeighborMix plus delayed DEC prototypes.

All non-scMAE losses start at zero and warm up. The scMAE loss always remains the
main objective.

## Benchmark rule

Use `benchmark.py`. The screen default is 80 epochs and three seeds:

```bash
python methods/DeepLearning/scMAEs/scmae_plus_mechanisms/benchmark.py \
  --stage screen \
  --variants 001_adaptive_mask 002_prototype_dec \
  --gpu 1
```

Screen success requires Melanoma_5K mean ARI above the original scMAE reference
ARI, currently about `0.668`; strong success is about `0.710`.

