# PCA+KMeans known-K

Traditional fixed-K baseline for the formal benchmark.

The runner applies the shared `scMAE_family` preprocessing, projects the data
with PCA, and evaluates KMeans with the benchmark-supplied number of clusters.
Labels are used only to score the final clustering result.

This is a known-K baseline, not a K-free clustering method.

