# 004 NeighborMix

This variant warms up the original scMAE embedding, builds a global mutual-KNN
graph, and adds a light NeighborMix reconstruction regularizer. It does not
replace the MLP encoder with a graph model.

Default mechanism:

```text
scMAE + warmup global KNN NeighborMix
```

The mix loss ramps from zero and uses only reliable mutual-nearest neighbors.

