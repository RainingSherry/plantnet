# rank14_cicl_cluster_aware_contrast_full

Independent-full scMAE candidate based on CICL.

Read sources:
- `00_scMAE改良方法整理总报告.md`, rank 14 CICL section.
- `02_整理索引.csv`, rank 14 row.
- `014_高_CICL_scRNA-seq_Data_Clustering_by_Cluster-aware_Iterative_Contrastive_Learning.pdf`.
- GitHub checks: index URL `https://github.com/Alunethy/CIRCLE` returned HTTP 200; report URL `https://github.com/WHY-17/Circle` returned HTTP 404.

Adaptation:
- Retains scMAE mask prediction and masked expression reconstruction.
- Adds two masked views, a projection head, and trainable Student-t cluster centers.
- After warmup, KMeans assigns pseudo labels and initializes cluster centers.
- Instance contrastive loss uses false-negative filtering so same high-confidence pseudo-cluster cells are not treated as negatives.
- Cluster-aware contrastive loss pulls together high-confidence same-pseudo-cluster cells.

NeighborMix relation:
- Independent and potentially complementary.
- This candidate does not mix cells; `mixed_cell_fraction=0.0`.

Smoke/screen outputs are candidate evidence only and must not be appended to `全benchmark结果.csv`.

