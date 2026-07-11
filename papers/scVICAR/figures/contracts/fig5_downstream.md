# Figure 5 contract: downstream biological utility

- **Core conclusion:** scVICAR embeddings preserve modest but reproducible low-label annotation utility, while marker recovery and marker-overlap annotation remain dataset dependent.
- **Archetype:** quantitative grid with a dominant low-label validation panel.
- **Target/output:** TCBB double-column figure, approximately 183 mm wide and at most 70 mm high; Python/matplotlib only; editable SVG and PDF plus 600-dpi TIFF and PNG preview.
- **Panel a:** dataset-level Recovery@100 for scMAE, scVICAR-F, and scVICAR-T.
- **Panel b:** dataset-level marker-overlap annotation macro-F1 for the same three models.
- **Panel c (hero):** frozen linear-probe macro-F1 at 10% and 30% labeled cells, with dataset-level paired trajectories or estimates and uncertainty that do not treat seeds or splits as biological replicates.
- **Evidence hierarchy:** panel c is the principal downstream evidence; panels a and b test whether the utility extends to marker-based interpretation and annotation.
- **Statistics:** aggregate split seeds within model seed and model seeds within dataset first; use datasets as independent units; report paired effects, hierarchical-bootstrap 95% CIs, paired permutation tests, and Holm correction in source tables/manuscript rather than decorating the plot with seed-level significance stars.
- **Source data:** frozen 108-row dataset/variant/task aggregate, downstream contrast table, configuration hashes, and SHA-256 manifest.
- **Image integrity:** no image panels; no smoothing, clipping, selective seed display, or best-seed selection.
- **Reviewer risks:** small absolute effects can look larger as relative percentages; ceiling effects may compress probe differences; marker endpoints may be heterogeneous or negative. Axes and prose must show absolute values, with relative percentages only as secondary interpretation.

## Final QA

- Complete aggregate: 108 unique dataset--variant--task/fraction rows (36 marker and 72 probe rows).
- All applicable marker and probe values are finite; invalid-marker cluster fraction is 0 and annotation/probe coverage is 1 for every aggregate row.
- Dataset-level points are displayed in every panel; panel c overlays cross-dataset means without treating seeds or splits as independent units.
- Exports: editable SVG/PDF, 600-dpi TIFF, and PNG preview; source-data CSVs pass their SHA-256 manifest.
- Visual inspection at final width found no clipped labels, overlap, redundant panel legends, or misleading significance decoration.
