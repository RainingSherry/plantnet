# Frozen statistical analysis plan (protocol_v1)

## Units and order of aggregation

The independent analysis unit is the dataset.  Cells, train/evaluation splits,
and model seeds are repeated measurements and are never counted as independent
biological replicates.  Downstream metrics are averaged over the five frozen
split seeds within each model run, then over the three model seeds within each
dataset.  Clustering metrics are averaged over model seeds within each dataset.

## Confirmatory comparisons

The three planned paired comparisons are:

1. scVICAR-F versus matched-backbone NoMix (scMAE objective);
2. scVICAR-T versus matched-backbone NoMix;
3. scVICAR-T versus scVICAR-F.

ARI is the primary clustering endpoint. NMI, ACC, macro-F1, silhouette,
runtime, and peak GPU memory are secondary. Marker recovery, marker-overlap
annotation, and low-label linear-probe endpoints are reported separately and
are not pooled into a composite score.

For each endpoint we report the paired dataset-level mean difference,
win/tie/loss counts, a hierarchical bootstrap 95% confidence interval, and a
paired sign-flip permutation p-value.  Holm correction is applied within the
family of the three planned comparisons for each endpoint.  Ties use an
absolute difference tolerance of `1e-12`.

## Reporting constraints

- All three model seeds are retained; there is no best-seed selection.
- Known-K KMeans uses the preregistered K and `n_init=20`.
- Fixed-resolution Leiden is secondary; no label-guided resolution sweep is
  permitted.
- `hrvatin_geo` is excluded from formal summaries because its label semantics
  are not an independent valid dataset.
- Development experiments on the original 16 datasets are disclosed as model
  development evidence in the supplement and are not combined with the six
  confirmatory datasets.
- If downstream evidence is null or adverse, the complete result is retained
  and the manuscript claim is narrowed.
