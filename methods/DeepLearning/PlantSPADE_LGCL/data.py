from dataclasses import dataclass
from typing import Optional, Tuple

import numpy as np

if not hasattr(np, "string_"):
    np.string_ = np.bytes_

import scanpy as sc
import scipy.sparse as sp
from sklearn.decomposition import TruncatedSVD
from sklearn.preprocessing import LabelEncoder, StandardScaler, normalize

from .utils import LABEL_CANDIDATES


@dataclass
class LGCLDatasetBundle:
    adata: sc.AnnData
    support: sp.csr_matrix
    amplitude: sp.csr_matrix
    global_embedding: np.ndarray
    labels: Optional[np.ndarray]
    label_names: Optional[np.ndarray]
    label_key: Optional[str]
    gene_names: np.ndarray
    input_mode: str
    support_density: float


def _ensure_csr(matrix) -> sp.csr_matrix:
    if sp.issparse(matrix):
        out = matrix.tocsr().astype(np.float32)
    else:
        out = sp.csr_matrix(np.asarray(matrix, dtype=np.float32))
    out.data = np.nan_to_num(out.data, nan=0.0, posinf=0.0, neginf=0.0)
    out.data[out.data < 0] = 0.0
    out.eliminate_zeros()
    out.sort_indices()
    return out


def _sample_values(matrix, max_rows: int = 256) -> np.ndarray:
    sample = matrix[: min(max_rows, matrix.shape[0])]
    if sp.issparse(sample):
        if sample.nnz == 0:
            return np.array([], dtype=np.float32)
        return sample.data.astype(np.float32, copy=False)
    return np.asarray(sample, dtype=np.float32).ravel()


def _looks_like_raw_counts(matrix) -> bool:
    values = _sample_values(matrix)
    values = values[np.isfinite(values)]
    if values.size == 0:
        return False
    values = values[: min(values.size, 100000)]
    return bool(np.all(values >= 0.0) and np.allclose(values, np.round(values), atol=1e-4))


def _infer_labels(adata: sc.AnnData) -> Tuple[Optional[np.ndarray], Optional[np.ndarray], Optional[str]]:
    for candidate in LABEL_CANDIDATES:
        if candidate in adata.obs.columns:
            raw = adata.obs[candidate].astype(str).to_numpy()
            encoder = LabelEncoder()
            labels = encoder.fit_transform(raw).astype(np.int64)
            return labels, encoder.classes_, candidate
    return None, None, None


def _select_source(adata: sc.AnnData, input_mode: str):
    if input_mode not in {"auto", "raw", "log1p"}:
        raise ValueError(f"Unknown input_mode: {input_mode}")

    if input_mode in {"auto", "raw"} and adata.raw is not None:
        return adata.raw.X, np.asarray(adata.raw.var_names), adata.raw.var.copy(), "raw"

    inferred = "raw" if _looks_like_raw_counts(adata.X) else "log1p"
    if input_mode == "raw" and not _looks_like_raw_counts(adata.X):
        raise ValueError("--input_mode raw was requested, but adata.X does not look like raw counts and adata.raw is absent.")
    if input_mode == "log1p":
        inferred = "log1p"
    return adata.X, np.asarray(adata.var_names), adata.var.copy(), inferred


def _tfidf_svd(
    amplitude: sp.csr_matrix,
    n_components: int,
    seed: int,
    n_iter: int = 7,
) -> np.ndarray:
    n_cells, n_genes = amplitude.shape
    max_components = max(1, min(n_components, n_cells - 1, n_genes - 1))
    x = amplitude.astype(np.float32).tocsr(copy=True)
    row_sum = np.asarray(x.sum(axis=1)).ravel().astype(np.float32)
    inv_row = np.divide(1.0, row_sum, out=np.zeros_like(row_sum), where=row_sum > 0)
    x = x.multiply(inv_row[:, None]).tocsr()

    df = np.diff(x.tocsc().indptr).astype(np.float32)
    idf = np.log1p(float(n_cells) / (1.0 + df)).astype(np.float32)
    x = x.multiply(idf).tocsr()

    svd = TruncatedSVD(n_components=max_components, n_iter=n_iter, random_state=seed)
    emb = svd.fit_transform(x).astype(np.float32)
    emb = np.nan_to_num(emb, nan=0.0, posinf=0.0, neginf=0.0)
    if emb.shape[1] > 1:
        emb = StandardScaler().fit_transform(emb).astype(np.float32)
    emb = normalize(emb, norm="l2", axis=1, copy=False).astype(np.float32)
    return emb


def load_lgcl_dataset(
    file_path: str,
    input_mode: str = "auto",
    n_top_genes: int = 2000,
    target_sum: float = 1e4,
    svd_dim: int = 32,
    svd_iter: int = 7,
    seed: int = 42,
) -> LGCLDatasetBundle:
    adata = sc.read_h5ad(file_path)
    source_x, gene_names, var, inferred_mode = _select_source(adata, input_mode)
    counts = _ensure_csr(source_x)

    work = sc.AnnData(X=counts.copy())
    work.obs_names = adata.obs_names.copy()
    work.var_names = gene_names.copy()
    work.obs = adata.obs.copy()
    work.var = var.copy()

    if inferred_mode == "raw":
        sc.pp.normalize_total(work, target_sum=target_sum)
        sc.pp.log1p(work)
    else:
        probe = _sample_values(work.X)
        if probe.size and float(np.nanmax(probe)) > 30.0:
            sc.pp.normalize_total(work, target_sum=target_sum)
            sc.pp.log1p(work)

    selected_idx = np.arange(counts.shape[1], dtype=np.int64)
    if n_top_genes and n_top_genes > 0 and work.n_vars > n_top_genes:
        sc.pp.highly_variable_genes(work, flavor="seurat", n_top_genes=n_top_genes, subset=True)
        selected = np.asarray(work.var_names)
        gene_to_idx = {gene: idx for idx, gene in enumerate(gene_names)}
        selected_idx = np.asarray([gene_to_idx[gene] for gene in selected], dtype=np.int64)
        gene_names = selected

    amplitude = _ensure_csr(work.X)
    support = counts[:, selected_idx]
    support = support.copy()
    support.data = np.ones_like(support.data, dtype=np.float32)
    support = support.astype(np.float32).tocsr()
    support.eliminate_zeros()
    support.sort_indices()

    labels, label_names, label_key = _infer_labels(work)
    global_embedding = _tfidf_svd(amplitude, n_components=svd_dim, seed=seed, n_iter=svd_iter)
    support_density = float(support.nnz / max(1, support.shape[0] * support.shape[1]))

    if label_key is not None:
        print(f"Detected labels from '{label_key}' with {len(label_names)} classes")
    else:
        print("Warning: no known label column found; evaluation will be skipped")

    return LGCLDatasetBundle(
        adata=work,
        support=support,
        amplitude=amplitude,
        global_embedding=global_embedding,
        labels=labels,
        label_names=label_names,
        label_key=label_key,
        gene_names=np.asarray(gene_names),
        input_mode=inferred_mode,
        support_density=support_density,
    )
