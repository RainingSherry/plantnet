from dataclasses import dataclass
from typing import Optional

import numpy as np

if not hasattr(np, "string_"):
    np.string_ = np.bytes_

import scanpy as sc
import scipy.sparse as sp
import torch
from sklearn.preprocessing import LabelEncoder
from torch.utils.data import DataLoader, Dataset

from .utils import LABEL_CANDIDATES


@dataclass
class SCDatasetBundle:
    adata: sc.AnnData
    values: np.ndarray
    support: np.ndarray
    labels: Optional[np.ndarray]
    label_names: Optional[np.ndarray]
    gene_names: np.ndarray
    x_max: float
    input_mode: str


class CellExpressionDataset(Dataset):
    def __init__(self, values: np.ndarray, support: np.ndarray):
        self.values = torch.tensor(values, dtype=torch.float32)
        self.support = torch.tensor(support, dtype=torch.float32)

    def __len__(self) -> int:
        return self.values.shape[0]

    def __getitem__(self, idx: int):
        return self.values[idx], self.support[idx]


class FullDataset(Dataset):
    def __init__(self, values: np.ndarray, support: np.ndarray, labels: Optional[np.ndarray]):
        self.values = torch.tensor(values, dtype=torch.float32)
        self.support = torch.tensor(support, dtype=torch.float32)
        self.labels = None if labels is None else torch.tensor(labels, dtype=torch.long)

    def __len__(self) -> int:
        return self.values.shape[0]

    def __getitem__(self, idx: int):
        label = -1 if self.labels is None else self.labels[idx]
        return self.values[idx], self.support[idx], label


def _to_dense(matrix) -> np.ndarray:
    if sp.issparse(matrix):
        return matrix.toarray()
    return np.asarray(matrix)


def _infer_labels(adata: sc.AnnData):
    label_col = None
    for candidate in LABEL_CANDIDATES:
        if candidate in adata.obs.columns:
            label_col = candidate
            break
    if label_col is None:
        return None, None
    labels = adata.obs[label_col].astype(str).to_numpy()
    encoder = LabelEncoder()
    encoded = encoder.fit_transform(labels)
    return encoded.astype(np.int64), encoder.classes_


def _prepare_matrix(adata: sc.AnnData, input_mode: str):
    if input_mode == "raw":
        if adata.raw is not None:
            X = _to_dense(adata.raw.X).astype(np.float32)
            gene_names = np.asarray(adata.raw.var_names)
        else:
            X = _to_dense(adata.X).astype(np.float32)
            gene_names = np.asarray(adata.var_names)
    else:
        X = _to_dense(adata.X).astype(np.float32)
        gene_names = np.asarray(adata.var_names)
    return X, gene_names


def load_sc_dataset(
    file_path: str,
    input_mode: str = "auto",
    n_top_genes: int = 2000,
    normalize_total: bool = True,
    target_sum: float = 1e4,
):
    adata = sc.read_h5ad(file_path)
    inferred_mode = input_mode
    if input_mode == "auto":
        X_sample = _to_dense(adata.X[: min(256, adata.n_obs)]).astype(np.float32)
        inferred_mode = "raw" if np.allclose(X_sample, np.round(X_sample)) else "log1p"

    X_input, gene_names = _prepare_matrix(adata, inferred_mode)
    support = (X_input > 0).astype(np.float32)

    work = sc.AnnData(X=X_input.copy())
    work.obs_names = adata.obs_names.copy()
    work.var_names = gene_names.copy()
    work.obs = adata.obs.copy()
    work.var = adata.raw.var.copy() if inferred_mode == "raw" and adata.raw is not None else adata.var.copy()

    if inferred_mode == "raw":
        if normalize_total:
            sc.pp.normalize_total(work, target_sum=target_sum)
        sc.pp.log1p(work)
    elif normalize_total:
        X_cur = _to_dense(work.X)
        if X_cur.max() > 30:
            sc.pp.normalize_total(work, target_sum=target_sum)
            sc.pp.log1p(work)

    if n_top_genes is not None and n_top_genes > 0 and work.n_vars > n_top_genes:
        sc.pp.highly_variable_genes(work, flavor="seurat", n_top_genes=n_top_genes, subset=True)
        selected = np.asarray(work.var_names)
        if inferred_mode == "raw" and adata.raw is not None:
            gene_to_idx = {g: i for i, g in enumerate(gene_names)}
        else:
            gene_to_idx = {g: i for i, g in enumerate(gene_names)}
        selected_idx = np.array([gene_to_idx[g] for g in selected], dtype=np.int64)
        support = support[:, selected_idx]
        gene_names = selected

    values = _to_dense(work.X).astype(np.float32)
    values = np.where(np.isfinite(values), values, 0.0)
    x_max = float(values.max()) if values.size > 0 else 1.0
    x_max = max(x_max, 1e-6)
    values = np.clip(values / x_max, 0.0, 1.0)

    labels, label_names = _infer_labels(work)

    return SCDatasetBundle(
        adata=work,
        values=values,
        support=support.astype(np.float32),
        labels=labels,
        label_names=label_names,
        gene_names=np.asarray(gene_names),
        x_max=x_max,
        input_mode=inferred_mode,
    )


def make_dataloader(values, support, labels=None, batch_size: int = 256, shuffle: bool = True):
    dataset = FullDataset(values, support, labels)
    return DataLoader(dataset, batch_size=batch_size, shuffle=shuffle, drop_last=False)
