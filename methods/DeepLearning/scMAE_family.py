from __future__ import annotations

import argparse
import os
import random
from dataclasses import dataclass
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import torch
from scipy.optimize import linear_sum_assignment
from sklearn.cluster import KMeans
from sklearn.metrics import (
    adjusted_rand_score,
    completeness_score,
    f1_score,
    fowlkes_mallows_score,
    homogeneity_score,
    normalized_mutual_info_score,
    v_measure_score,
)
from sklearn.preprocessing import LabelEncoder
from torch.utils.data import DataLoader, Dataset

from methods.DeepLearning.PlantSPADE_LGCL.utils import save_json


LABEL_CANDIDATES = [
    "maintype",
    "cell_type",
    "Celltype",
    "celltype",
    "label",
    "labels",
    "cell_label",
    "Cluster",
    "cluster",
    "clusters",
    "Seurat_clusters",
    "true_label",
]


@dataclass
class DataBundle:
    adata: sc.AnnData
    data: np.ndarray
    labels: np.ndarray
    label_names: np.ndarray
    label_key: str
    gene_names: np.ndarray
    profile: dict
    preprocess_config: dict


class IndexedExpressionDataset(Dataset):
    def __init__(self, data: np.ndarray, labels: np.ndarray):
        self.data = torch.as_tensor(data, dtype=torch.float32)
        self.labels = torch.as_tensor(labels, dtype=torch.long)

    def __len__(self) -> int:
        return int(self.data.shape[0])

    def __getitem__(self, idx: int):
        return int(idx), self.data[idx], self.labels[idx]


def str2bool(value):
    if isinstance(value, bool):
        return value
    value = str(value).strip().lower()
    if value in {"1", "true", "t", "yes", "y"}:
        return True
    if value in {"0", "false", "f", "no", "n"}:
        return False
    raise argparse.ArgumentTypeError(f"Expected true/false, got {value!r}")


def set_seed(seed: int) -> None:
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed(seed)
        torch.cuda.manual_seed_all(seed)
        torch.backends.cudnn.deterministic = True
        torch.backends.cudnn.benchmark = False


def get_device(gpu: int, no_cuda: bool) -> torch.device:
    if no_cuda or not torch.cuda.is_available():
        return torch.device("cpu")
    visible = os.environ.get("CUDA_VISIBLE_DEVICES")
    if visible:
        visible_ids = [item.strip() for item in visible.split(",") if item.strip()]
        if set(visible_ids).intersection({"0", "7"}):
            raise ValueError("CUDA_VISIBLE_DEVICES includes forbidden physical GPU 0 or 7.")
        if len(visible_ids) == 1:
            return torch.device("cuda:0")
        if str(gpu) in visible_ids:
            return torch.device(f"cuda:{visible_ids.index(str(gpu))}")
        if 0 <= gpu < len(visible_ids):
            return torch.device(f"cuda:{gpu}")
        raise ValueError(f"--gpu {gpu} is outside isolated CUDA_VISIBLE_DEVICES={visible!r}.")
    if gpu in {0, 7}:
        raise ValueError("Physical GPU 0 and GPU 7 are forbidden. Use 1,2,3,4,5,6 or --no_cuda.")
    return torch.device(f"cuda:{gpu}")


def ensure_csr(matrix) -> sp.csr_matrix:
    if sp.issparse(matrix):
        out = matrix.tocsr().astype(np.float32)
    else:
        out = sp.csr_matrix(np.asarray(matrix, dtype=np.float32))
    out.data = np.nan_to_num(out.data, nan=0.0, posinf=0.0, neginf=0.0)
    out.data[out.data < 0.0] = 0.0
    out.eliminate_zeros()
    out.sort_indices()
    return out


def sample_values(matrix, max_rows: int = 256) -> np.ndarray:
    sample = matrix[: min(max_rows, matrix.shape[0])]
    if sp.issparse(sample):
        return sample.data.astype(np.float32, copy=False) if sample.nnz else np.array([], dtype=np.float32)
    return np.asarray(sample, dtype=np.float32).ravel()


def looks_like_raw_counts(matrix) -> bool:
    values = sample_values(matrix)
    values = values[np.isfinite(values)]
    if values.size == 0:
        return False
    values = values[: min(values.size, 100000)]
    return bool(np.all(values >= 0.0) and np.allclose(values, np.round(values), atol=1e-4))


def var_gene_names(var) -> tuple[np.ndarray, object]:
    var = var.copy()
    gene_names = np.asarray(var.index).astype(str)
    for name_col in ["gene_name", "features", "gene_symbols", "symbol", "_index"]:
        if name_col in var.columns:
            candidate = var[name_col].astype(str).to_numpy()
            if len(candidate) == len(gene_names) and len(np.unique(candidate)) == len(candidate):
                gene_names = candidate
                break
    var.index = gene_names
    return gene_names, var


def select_count_source(adata: sc.AnnData, input_mode: str):
    if input_mode in {"auto", "raw"} and "counts" in adata.layers:
        gene_names, var = var_gene_names(adata.var)
        return adata.layers["counts"], gene_names, var, "layers[counts]", "raw"
    if input_mode in {"auto", "raw"} and adata.raw is not None:
        gene_names, var = var_gene_names(adata.raw.var)
        return adata.raw.X, gene_names, var, "adata.raw.X", "raw"
    inferred = "raw" if looks_like_raw_counts(adata.X) else "log1p"
    if input_mode == "raw" and inferred != "raw":
        raise ValueError("--input_mode raw was requested, but no raw-looking X/raw/layers[counts] source is available.")
    if input_mode == "log1p":
        inferred = "log1p"
    gene_names, var = var_gene_names(adata.var)
    source = "adata.X" if inferred == "raw" else "adata.X_log1p_fallback"
    return adata.X, gene_names, var, source, inferred


def resolve_labels(adata: sc.AnnData, label_key: str):
    if label_key and label_key != "auto":
        if label_key not in adata.obs.columns:
            raise KeyError(f"Configured label_key={label_key!r} is absent. Available obs columns: {list(adata.obs.columns)}")
        raw = adata.obs[label_key].astype(str).to_numpy()
        encoder = LabelEncoder()
        return encoder.fit_transform(raw).astype(np.int64), encoder.classes_, label_key
    for candidate in LABEL_CANDIDATES:
        if candidate in adata.obs.columns:
            raw = adata.obs[candidate].astype(str).to_numpy()
            encoder = LabelEncoder()
            return encoder.fit_transform(raw).astype(np.int64), encoder.classes_, candidate
    raise KeyError(f"No label column found. Available obs columns: {list(adata.obs.columns)}")


def dense_float32(matrix) -> np.ndarray:
    if sp.issparse(matrix):
        arr = matrix.toarray()
    else:
        arr = np.asarray(matrix)
    arr = np.asarray(arr, dtype=np.float32)
    return np.nan_to_num(arr, nan=0.0, posinf=0.0, neginf=0.0)


def normalize_total_log1p(matrix, target_sum: float) -> sp.csr_matrix:
    x = ensure_csr(matrix).astype(np.float32)
    row_sum = np.asarray(x.sum(axis=1)).ravel().astype(np.float32)
    scale = np.divide(
        float(target_sum),
        row_sum,
        out=np.zeros_like(row_sum, dtype=np.float32),
        where=row_sum > 0.0,
    )
    x = x.multiply(scale[:, None]).tocsr()
    x.data = np.log1p(x.data).astype(np.float32, copy=False)
    x.eliminate_zeros()
    return x


def load_scmae_dataset(
    file_path: str,
    input_mode: str,
    n_top_genes: int,
    target_sum: float,
    scale_input: bool,
    label_key: str,
    seed: int,
) -> DataBundle:
    adata = sc.read_h5ad(file_path)
    source_x, gene_names, var, counts_source, inferred_mode = select_count_source(adata, input_mode)
    counts = ensure_csr(source_x)
    work = sc.AnnData(X=counts.copy(), obs=adata.obs.copy(), var=var.copy())
    work.obs_names = adata.obs_names.copy()
    work.var_names = gene_names.copy()

    if inferred_mode == "raw":
        work.X = normalize_total_log1p(work.X, target_sum=target_sum)
    elif sample_values(work.X).size and float(np.nanmax(sample_values(work.X))) > 30.0:
        work.X = normalize_total_log1p(work.X, target_sum=target_sum)

    if n_top_genes and n_top_genes > 0 and work.n_vars > n_top_genes:
        sc.pp.highly_variable_genes(work, flavor="seurat", n_top_genes=n_top_genes, subset=True)

    if scale_input:
        sc.pp.scale(work)

    data = dense_float32(work.X)
    labels, label_names, resolved_label_key = resolve_labels(work, label_key)
    gene_names = np.asarray(work.var_names).astype(str)

    label_counts = {
        str(key): int(value)
        for key, value in work.obs[resolved_label_key].astype(str).value_counts(dropna=False).sort_index().items()
    }
    profile = {
        "dataset_name": Path(file_path).stem,
        "n_cells": int(work.n_obs),
        "n_genes_original": int(adata.n_vars),
        "n_genes": int(work.n_vars),
        "label_key": resolved_label_key,
        "n_cell_types": int(len(label_names)),
        "cell_type_counts": label_counts,
        "counts_source": counts_source,
        "input_mode": inferred_mode,
        "has_raw": bool(adata.raw is not None),
        "has_layers_counts": bool("counts" in adata.layers),
        "scale_input": bool(scale_input),
    }
    preprocess_config = {
        "file_path": str(file_path),
        "counts_source": counts_source,
        "input_mode": inferred_mode,
        "normalization": f"normalize_total(target_sum={target_sum}) + log1p when raw",
        "hvg": {"n_top_genes": int(n_top_genes), "flavor": "seurat"},
        "scale_input": bool(scale_input),
        "selected_n_genes": int(work.n_vars),
        "seed": int(seed),
        "label_key": resolved_label_key,
    }
    return DataBundle(
        adata=work,
        data=data,
        labels=labels,
        label_names=label_names,
        label_key=resolved_label_key,
        gene_names=gene_names,
        profile=profile,
        preprocess_config=preprocess_config,
    )


def apply_scmae_noise(x: torch.Tensor, mask_ratio: float) -> tuple[torch.Tensor, torch.Tensor]:
    should_swap = torch.bernoulli(float(mask_ratio) * torch.ones_like(x))
    if x.shape[0] <= 1:
        replacement = x
    else:
        replacement = x[torch.randperm(x.shape[0], device=x.device)]
    corrupted = torch.where(should_swap.bool(), replacement, x)
    mask = (corrupted != x).float()
    return corrupted, mask


@torch.no_grad()
def extract_embedding(model, loader: DataLoader, device: torch.device) -> tuple[np.ndarray, np.ndarray]:
    model.eval()
    embeddings = []
    labels = []
    for _, x, y in loader:
        z = model.feature(x.to(device))
        embeddings.append(z.detach().cpu().numpy())
        labels.append(y.numpy())
    emb = np.concatenate(embeddings, axis=0).astype(np.float32)
    labels_np = np.concatenate(labels, axis=0).astype(np.int64)
    emb = np.nan_to_num(emb, nan=0.0, posinf=0.0, neginf=0.0)
    return emb, labels_np


def save_embedding_h5(path: Path, embedding: np.ndarray, labels: np.ndarray) -> None:
    with h5py.File(path, "w") as handle:
        handle.create_dataset("X", data=embedding.astype(np.float32))
        handle.create_dataset("labels", data=labels.astype(np.int64))


def best_map(y_true: np.ndarray, y_pred: np.ndarray) -> np.ndarray:
    true_values = np.unique(y_true)
    pred_values = np.unique(y_pred)
    n = max(len(true_values), len(pred_values))
    counts = np.zeros((n, n), dtype=np.int64)
    for i, true_label in enumerate(true_values):
        for j, pred_label in enumerate(pred_values):
            counts[i, j] = int(np.sum((y_true == true_label) & (y_pred == pred_label)))
    rows, cols = linear_sum_assignment(-counts)
    mapped = np.zeros_like(y_pred, dtype=np.int64)
    for row, col in zip(rows, cols):
        if row < len(true_values) and col < len(pred_values):
            mapped[y_pred == pred_values[col]] = true_values[row]
    return mapped


def compute_kmeans_metrics(labels: np.ndarray, pred: np.ndarray) -> tuple[dict, np.ndarray]:
    mapped = best_map(labels, pred)
    metrics = {
        "acc": float(np.mean(mapped == labels)),
        "nmi": float(normalized_mutual_info_score(labels, pred)),
        "ari": float(adjusted_rand_score(labels, pred)),
        "f1_macro": float(f1_score(labels, mapped, average="macro", zero_division=0)),
        "fmi": float(fowlkes_mallows_score(labels, pred)),
        "v_measure": float(v_measure_score(labels, pred)),
        "homogeneity": float(homogeneity_score(labels, pred)),
        "completeness": float(completeness_score(labels, pred)),
        "n_pred_clusters": int(len(np.unique(pred))),
        "silhouette": float("nan"),
        "protocol": "fixed",
        "cluster_method": "kmeans_known_k",
        "uses_known_k": True,
    }
    return metrics, mapped.astype(np.int64)


def write_kmeans_known_k_outputs(
    output_dir: Path,
    dataset: str,
    method: str,
    seed: int,
    embedding: np.ndarray,
    labels: np.ndarray,
    n_clusters: int,
    extra: dict,
) -> dict:
    pred = KMeans(n_clusters=n_clusters, n_init=20, random_state=seed).fit_predict(embedding)
    metrics, mapped = compute_kmeans_metrics(labels, pred.astype(np.int64))
    fixed = {"kmeans_known_k": metrics}
    row = {
        "dataset": dataset,
        "method": method,
        "seed": int(seed),
        **extra,
        **metrics,
    }
    pd.DataFrame([row]).to_csv(output_dir / "eval_fixed.csv", index=False)
    payload = {
        "dataset": dataset,
        "method": method,
        "seed": int(seed),
        "n_clusters": int(n_clusters),
        "fixed": fixed,
        "oracle": {},
        "sweep": [],
    }
    save_json(payload, str(output_dir / "eval_metrics.json"))
    np.save(output_dir / "eval_kmeans_known_k.npy", pred.astype(np.int64))
    np.save(output_dir / "eval_kmeans_known_k_mapped.npy", mapped)
    return {"fixed": fixed, "preds": {"kmeans_known_k": pred.astype(np.int64)}}
