from __future__ import annotations

import numpy as np
from sklearn.cluster import KMeans


def compute_gene_modules(expr: np.ndarray, n_modules: int, seed: int = 42,
                         max_cells: int = 8000) -> np.ndarray:
    """Cluster genes into co-expression modules (solves 'genes have no order').

    Gene features = each gene's expression profile across (subsampled) cells,
    standardized. KMeans on gene-space groups co-expressed genes together.
    Returns module_of: int array (G,) with module id in [0, n_modules).

    Decoupled + precomputed once from training data (like the reliability graph),
    so it does not depend on the live embedding.
    """
    n, g = expr.shape
    rng = np.random.default_rng(seed)
    if n > max_cells:
        idx = rng.choice(n, size=max_cells, replace=False)
        sub = expr[idx]
    else:
        sub = expr
    # gene vectors: (g, cells); standardize each gene
    gv = sub.T.astype(np.float64)
    gv = (gv - gv.mean(axis=1, keepdims=True)) / (gv.std(axis=1, keepdims=True) + 1e-8)
    m = min(int(n_modules), g)
    km = KMeans(n_clusters=m, n_init=10, random_state=seed).fit(gv)
    module_of = km.labels_.astype(np.int64)
    # guard: ensure contiguous ids in [0, m)
    _, remap = np.unique(module_of, return_inverse=True)
    return remap.astype(np.int64)
