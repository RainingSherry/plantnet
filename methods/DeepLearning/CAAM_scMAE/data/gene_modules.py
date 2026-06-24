from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import scipy.sparse as sp
from sklearn.cluster import KMeans
from sklearn.decomposition import TruncatedSVD
from sklearn.preprocessing import StandardScaler


def build_gene_modules(
    x: np.ndarray,
    n_gene_modules: int,
    module_svd_dim: int,
    module_seed: int,
    save_dir: Path | None = None,
) -> tuple[np.ndarray, sp.csr_matrix]:
    n_genes = int(x.shape[1])
    n_modules = max(1, min(int(n_gene_modules), n_genes))
    x_copy = np.asarray(x, dtype=np.float32)
    x_scaled = StandardScaler(with_mean=True, with_std=True).fit_transform(x_copy)
    dim = max(1, min(int(module_svd_dim), min(x_scaled.shape) - 1)) if min(x_scaled.shape) > 1 else 1
    if n_genes <= n_modules:
        ids = np.arange(n_genes, dtype=np.int64)
        n_modules = n_genes
    else:
        svd = TruncatedSVD(n_components=dim, random_state=module_seed)
        gene_embedding = svd.fit_transform(x_scaled.T)
        ids = KMeans(n_clusters=n_modules, n_init=20, random_state=module_seed).fit_predict(gene_embedding)
        ids = ids.astype(np.int64)
    rows = np.arange(n_genes, dtype=np.int64)
    data = np.ones(n_genes, dtype=np.float32)
    assignment = sp.csr_matrix((data, (rows, ids)), shape=(n_genes, n_modules), dtype=np.float32)
    if save_dir is not None:
        save_dir.mkdir(parents=True, exist_ok=True)
        np.save(save_dir / "gene_module_ids.npy", ids)
        sp.save_npz(save_dir / "gene_module_assignment.npz", assignment)
        with open(save_dir / "gene_module_builder.json", "w", encoding="utf-8") as handle:
            json.dump(
                {
                    "n_genes": n_genes,
                    "n_gene_modules": int(n_modules),
                    "module_svd_dim": int(dim),
                    "module_seed": int(module_seed),
                    "label_free": True,
                },
                handle,
                indent=2,
            )
    return ids, assignment


def normalized_assignment_dense(assignment: sp.csr_matrix) -> np.ndarray:
    arr = assignment.toarray().astype(np.float32)
    denom = arr.sum(axis=0, keepdims=True)
    denom[denom == 0.0] = 1.0
    return arr / denom

