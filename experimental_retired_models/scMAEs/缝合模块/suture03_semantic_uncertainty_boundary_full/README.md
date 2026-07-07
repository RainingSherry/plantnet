# suture03_semantic_uncertainty_boundary_full

This candidate adapts the Semantic Information Decoupling idea into a cell-level core / boundary / rare-risk latent gate for scMAE.

It keeps scMAE mask prediction and masked expression reconstruction. The new mechanism is a lightweight latent adapter: core cells may be refined, rare-risk cells receive a separate branch, and boundary-like cells are guarded by a conservative delta penalty.

The reference SID module is image-based; this implementation does not use 2D image convolution or gene-vector reshaping. Boundary and rare-risk soft targets are computed from expression-space neighbor distance, dropout, and expression concentration.

NeighborMix relation: independent. No cells are mixed, so `mixed_cell_fraction=0.0`.
