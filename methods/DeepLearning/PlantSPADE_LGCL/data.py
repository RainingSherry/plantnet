from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Optional, Tuple

import numpy as np

if not hasattr(np, "string_"):
    np.string_ = np.bytes_

import scanpy as sc
import scipy.sparse as sp
from sklearn.decomposition import TruncatedSVD
from sklearn.preprocessing import LabelEncoder, StandardScaler, normalize

from .utils import LABEL_CANDIDATES, ensure_dir, save_json


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
    counts_source: str
    support_density: float
    profile: Dict
    preprocess_config: Dict


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


def _infer_labels(
    adata: sc.AnnData,
    label_key: Optional[str] = None,
) -> Tuple[Optional[np.ndarray], Optional[np.ndarray], Optional[str]]:
    if label_key and label_key != "auto":
        if label_key not in adata.obs.columns:
            raise KeyError(f"Configured label_key={label_key!r} is absent. Available obs columns: {list(adata.obs.columns)}")
        raw = adata.obs[label_key].astype(str).to_numpy()
        encoder = LabelEncoder()
        labels = encoder.fit_transform(raw).astype(np.int64)
        return labels, encoder.classes_, label_key
    for candidate in LABEL_CANDIDATES:
        if candidate in adata.obs.columns:
            raw = adata.obs[candidate].astype(str).to_numpy()
            encoder = LabelEncoder()
            labels = encoder.fit_transform(raw).astype(np.int64)
            return labels, encoder.classes_, candidate
    return None, None, None


def _subsample_adata(
    adata: sc.AnnData,
    label_key: Optional[str],
    per_class_max: Optional[int],
    fallback_max: Optional[int],
    seed: int,
) -> Tuple[sc.AnnData, Dict]:
    per_class_max = int(per_class_max or 0)
    fallback_max = int(fallback_max or 0)
    if per_class_max <= 0 and (fallback_max <= 0 or adata.n_obs <= fallback_max):
        return adata, {"enabled": False}

    rng = np.random.default_rng(seed)
    labels, _, resolved_label_key = _infer_labels(adata, label_key)
    selected = []
    mode = "none"
    if per_class_max > 0 and labels is not None:
        mode = "per_class"
        for value in np.unique(labels):
            idx = np.flatnonzero(labels == value)
            if idx.size > per_class_max:
                idx = rng.choice(idx, size=per_class_max, replace=False)
            selected.append(idx)
        selected_idx = np.sort(np.concatenate(selected)) if selected else np.arange(adata.n_obs)
    elif fallback_max > 0 and adata.n_obs > fallback_max:
        mode = "fallback"
        selected_idx = np.sort(rng.choice(adata.n_obs, size=fallback_max, replace=False))
    else:
        return adata, {"enabled": False}

    out = adata[selected_idx].copy()
    return out, {
        "enabled": True,
        "mode": mode,
        "label_key": resolved_label_key,
        "original_n_cells": int(adata.n_obs),
        "subsampled_n_cells": int(out.n_obs),
        "per_class_max": int(per_class_max),
        "fallback_max": int(fallback_max),
        "seed": int(seed),
    }


def _var_gene_names(var) -> Tuple[np.ndarray, object]:
    var = var.copy()
    gene_names = np.asarray(var.index).astype(str)
    for name_col in ["features", "gene_name", "gene_symbols", "symbol", "_index"]:
        if name_col in var.columns:
            candidate = var[name_col].astype(str).to_numpy()
            if len(candidate) == len(gene_names) and len(np.unique(candidate)) == len(candidate):
                gene_names = candidate
                break
    var.index = gene_names
    return gene_names, var


def _select_count_source(adata: sc.AnnData, input_mode: str):
    if input_mode not in {"auto", "raw", "log1p"}:
        raise ValueError(f"Unknown input_mode: {input_mode}")

    if input_mode in {"auto", "raw"} and "counts" in adata.layers:
        gene_names, var = _var_gene_names(adata.var)
        return adata.layers["counts"], gene_names, var, "layers[counts]", "raw"

    if input_mode in {"auto", "raw"} and adata.raw is not None:
        gene_names, var = _var_gene_names(adata.raw.var)
        return adata.raw.X, gene_names, var, "adata.raw.X", "raw"

    inferred = "raw" if _looks_like_raw_counts(adata.X) else "log1p"
    if input_mode == "raw" and not _looks_like_raw_counts(adata.X):
        raise ValueError("--input_mode raw was requested, but adata.X does not look like raw counts and adata.raw is absent.")
    if input_mode == "log1p":
        inferred = "log1p"
    gene_names, var = _var_gene_names(adata.var)
    source = "adata.X" if inferred == "raw" else "adata.X_log1p_fallback"
    return adata.X, gene_names, var, source, inferred


def _quantiles(values: np.ndarray) -> Dict[str, float]:
    values = np.asarray(values, dtype=np.float64)
    values = values[np.isfinite(values)]
    if values.size == 0:
        return {key: 0.0 for key in ["min", "q05", "q25", "median", "q75", "q95", "max", "mean"]}
    qs = np.quantile(values, [0.05, 0.25, 0.5, 0.75, 0.95])
    return {
        "min": float(np.min(values)),
        "q05": float(qs[0]),
        "q25": float(qs[1]),
        "median": float(qs[2]),
        "q75": float(qs[3]),
        "q95": float(qs[4]),
        "max": float(np.max(values)),
        "mean": float(np.mean(values)),
    }


def profile_anndata(
    adata: sc.AnnData,
    dataset_name: str,
    label_key: Optional[str] = None,
    input_mode: str = "auto",
) -> Dict:
    selected_x, _, _, counts_source, inferred_mode = _select_count_source(adata, input_mode)
    counts = _ensure_csr(selected_x)
    labels, label_names, resolved_label_key = _infer_labels(adata, label_key)

    per_cell_nnz = np.diff(counts.indptr)
    library_size = np.asarray(counts.sum(axis=1)).ravel()
    label_counts = {}
    if resolved_label_key is not None:
        label_counts = {
            str(key): int(value)
            for key, value in adata.obs[resolved_label_key].astype(str).value_counts(dropna=False).sort_index().items()
        }

    obs_fields = set(map(str, adata.obs.columns))
    field_candidates = {
        "batch": ["batch", "Batch", "sample", "sample_id", "donor", "orig.ident"],
        "tissue": ["tissue", "Tissue", "organ", "Organ"],
        "genotype": ["genotype", "Genotype", "ecotype", "cultivar"],
        "condition": ["condition", "Condition", "treatment", "Treatment", "stage", "Stage"],
    }
    present_fields = {
        key: sorted([field for field in candidates if field in obs_fields])
        for key, candidates in field_candidates.items()
    }

    return {
        "dataset_name": dataset_name,
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "label_key": resolved_label_key,
        "n_cell_types": int(len(label_names)) if label_names is not None else 0,
        "cell_type_counts": label_counts,
        "x_looks_like_raw_counts": bool(_looks_like_raw_counts(adata.X)),
        "counts_source": counts_source,
        "counts_source_looks_like_raw_counts": bool(_looks_like_raw_counts(counts)),
        "input_mode": inferred_mode,
        "has_raw": bool(adata.raw is not None),
        "has_layers_counts": bool("counts" in adata.layers),
        "has_layers_log1p_norm": bool("log1p_norm" in adata.layers),
        "layers": sorted(map(str, adata.layers.keys())),
        "sparsity": float(1.0 - counts.nnz / max(1, counts.shape[0] * counts.shape[1])),
        "per_cell_nonzero_genes": _quantiles(per_cell_nnz),
        "library_size": _quantiles(library_size),
        "obs_fields_present": present_fields,
    }


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
    label_key: Optional[str] = None,
    subsample_per_class_max: Optional[int] = None,
    subsample_fallback_max: Optional[int] = None,
) -> LGCLDatasetBundle:
    adata = sc.read_h5ad(file_path)
    adata, subsample_config = _subsample_adata(
        adata,
        label_key=label_key,
        per_class_max=subsample_per_class_max,
        fallback_max=subsample_fallback_max,
        seed=seed,
    )
    dataset_name = Path(file_path).stem
    profile = profile_anndata(adata, dataset_name=dataset_name, label_key=label_key, input_mode=input_mode)
    source_x, gene_names, var, counts_source, inferred_mode = _select_count_source(adata, input_mode)
    counts = _ensure_csr(source_x)

    work = sc.AnnData(X=counts.copy())
    work.obs_names = adata.obs_names.copy()
    work.obs = adata.obs.copy()
    work.var = var.copy()
    work.var_names = gene_names.copy()

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

    labels, label_names, resolved_label_key = _infer_labels(work, label_key)
    global_embedding = _tfidf_svd(amplitude, n_components=svd_dim, seed=seed, n_iter=svd_iter)
    support_density = float(support.nnz / max(1, support.shape[0] * support.shape[1]))

    preprocess_config = {
        "file_path": str(file_path),
        "counts_source": counts_source,
        "input_mode": inferred_mode,
        "support_definition": "M_cg = 1[count_cg > 0]",
        "amplitude_definition": f"normalize_total(target_sum={target_sum}) + log1p on selected genes",
        "n_top_genes": int(n_top_genes) if n_top_genes is not None else 0,
        "selected_n_genes": int(len(gene_names)),
        "target_sum": float(target_sum),
        "svd_dim": int(global_embedding.shape[1]),
        "svd_iter": int(svd_iter),
        "seed": int(seed),
        "label_key": resolved_label_key,
        "subsample": subsample_config,
    }

    if resolved_label_key is not None:
        print(f"Detected labels from '{resolved_label_key}' with {len(label_names)} classes")
    else:
        print("Warning: no known label column found; evaluation will be skipped")

    return LGCLDatasetBundle(
        adata=work,
        support=support,
        amplitude=amplitude,
        global_embedding=global_embedding,
        labels=labels,
        label_names=label_names,
        label_key=resolved_label_key,
        gene_names=np.asarray(gene_names),
        input_mode=inferred_mode,
        counts_source=counts_source,
        support_density=support_density,
        profile=profile,
        preprocess_config=preprocess_config,
    )


def write_dataset_artifacts(bundle: LGCLDatasetBundle, output_dir: str) -> None:
    output = Path(ensure_dir(output_dir))
    save_json(bundle.profile, str(output / "dataset_profile.json"))
    save_json(bundle.preprocess_config, str(output / "preprocess_config.json"))
    with open(output / "selected_genes.txt", "w", encoding="utf-8") as handle:
        for gene in bundle.gene_names:
            handle.write(f"{gene}\n")
