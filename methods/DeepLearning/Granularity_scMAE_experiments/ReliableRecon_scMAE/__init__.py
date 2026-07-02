"""ReliableRecon-scMAE — backbone-level reliability-weighted reconstruction (NEGATIVE RESULT).

Reconstructed 2026-07-01 for the record. This line was a clean negative result
(see ../EXPERIMENT_LOG.md §6). Kept so the mechanism and its failure reason are
not lost.

Idea: attack the RECONSTRUCTION task (upstream of clustering, since Macosko is
backbone-level). Per-cell-per-gene LOCAL reliability r_ig from a DECOUPLED
raw-data PCA-KNN graph (variance of gene g among cell i's neighbors, relative to
per-gene median, floor 0.2). Use r to precision-weight the recon loss.
LOCAL not global (a marker is globally high-variance but locally stable -> stays
reliable; global variance would kill markers). lambda blends to vanilla.

RESULT: lambda=0 vs lambda=1, ARI: Melanoma 0.564->0.492, Quake 0.538->0.521,
Macosko 0.264->0.173. Monotonic degradation on all three.
WHY: reconstruction fidelity != clustering quality. Down-weighting locally-noisy
genes discards exactly the boundary/transition/rare discriminative signal.
"""
