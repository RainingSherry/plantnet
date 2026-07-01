# rank15_scagc_adaptive_graph_full

Independent-full scMAE candidate based on scAGC.

Read sources:
- `00_scMAE改良方法整理总报告.md`, rank 15 section.
- `02_整理索引.csv`, rank 15 row.
- `015_高_scAGC_Learning_Adaptive_Cell_Graphs_with_Contrastive_Guidance_for_Single-Cell_Clustering.pdf`.

No GitHub URL is listed in the index, so this implementation is reconstructed from the paper and report.

Adaptation:
- Retains scMAE mask prediction and masked expression reconstruction.
- Builds a local KNN graph from encoder input and estimates edge confidence by distance.
- Uses a shallow residual graph adapter, not a deep GNN stack.
- Adds edge reconstruction, neighbor consistency, and DropEdge robustness losses.
- Optional local mix is confidence-gated and vetoed by boundary/rare risk; it is not forced on all cells.

NeighborMix relation:
- Complementary.
- Mix strength goes to 0 when edge confidence is low or rare-risk is high.
- `mixed_cell_fraction` is reported in `diagnostics.json`.

Smoke/screen outputs are candidate evidence only and must not be appended to `全benchmark结果.csv`.

