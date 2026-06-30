# 001 Adaptive Mask

This variant keeps the original scMAE encoder, mask predictor, decoder,
cross-cell replacement corruption, and KMeans-known-k evaluation. The only
default change is replacing uniform Bernoulli gene masking with a gene-statistic
adaptive probability vector.

Default mechanism:

```text
scMAE + variance_adaptive mask
```

All prototype, SwAV, NeighborMix, count-aware, and fuzzy-boundary mechanisms are
off by default.

