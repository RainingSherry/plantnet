# rank13_masked_sc_cluster_target_full

Independent-full scMAE candidate based on Masked Modeling for Single-cell Clustering.

Read sources:
- `00_scMAE改良方法整理总报告.md`, rank 13 section.
- `02_整理索引.csv`, rank 13 row.
- `013_高_Masked_Modeling_for_Single-cell_Clustering_of_scRNA-seq_Data.pdf`.

No GitHub URL is listed in the index, so this implementation is reconstructed from the paper and report.

Adaptation:
- Retains scMAE mask prediction and masked expression reconstruction.
- Warmup trains only the scMAE body.
- After warmup, KMeans initializes trainable cluster centers on the learned embedding.
- A DEC-style target distribution supplies high-confidence cluster-aware KL supervision.
- Low-confidence/boundary cells are gated out of the hard cluster term.

NeighborMix relation:
- Independent and potentially complementary.
- This candidate does not mix cells; `mixed_cell_fraction=0.0`.

Smoke/screen results are candidate evidence only and must not be appended to `全benchmark结果.csv`.

