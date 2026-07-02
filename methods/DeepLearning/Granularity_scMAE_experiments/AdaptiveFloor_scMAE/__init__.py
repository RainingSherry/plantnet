"""AdaptiveFloor-scMAE — DEC clustering + PRINCIPLED adaptive per-dimension variance floor.

Context: a plain VICReg std-hinge (fixed std>=1) already rescues DEC on
fine-grained scRNA (see AdaptiveSwitch results: Macosko 0.343->0.695), but the
std-floor itself is NOT novel (it is VICReg's variance term). The novel
contribution here is a THEORY-DRIVEN ADAPTIVE floor.

Theory (why DEC collapses, why a per-dim std-floor is the exact fix):
  DEC minimizes KL(P||Q), P = sharpen(Q). The gradient pulls points toward
  centroids and centroids toward members -> a CONTRACTION on the latent. The
  encoder has free scale (it can shrink ||z|| while keeping softmax assignments
  ~unchanged). On fine-grained data (large k, subtypes close in expr space),
  most centroids crowd a small region -> per-dimension marginal variance Var_d
  contracts toward 0 on many dims -> effective dim drops (128->~60) -> KMeans
  cannot separate fine subtypes. Along each dim, the DEC update is a contraction
  map whose STABLE fixed point is Var_d=0 (degenerate). A per-dim std-floor adds
  a repulsion that MOVES the stable fixed point to Var_d = target_d^2.
  Covariance-decorrelation only penalizes OFF-diagonal cov -> blind to isotropic
  scale collapse. KoLeo acts on L2-normalized space -> does not control raw
  per-dim variance KMeans uses. (Both empirically fail: 0.33/0.36 vs 0.695.)

Adaptive floor (the theory says the floor should restore the data's NATURAL
scale, not a magic 1.0):
  floor_mode="ref"  : target_d = pre-DEC (post-warmup, reconstruction-only)
                      per-dim std. "Do not let DEC compress any dim below the
                      reconstruction geometry's own spread." Parameter-free,
                      per-dimension, self-adaptive to each run's warmup scale
                      -> also expected to reduce Macosko seed variance.
  floor_mode="fixed": legacy VICReg std>=1 (baseline for the ablation).
"""
