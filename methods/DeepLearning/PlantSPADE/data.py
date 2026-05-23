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


# ---------- 数据结构 ----------

@dataclass
class SCDatasetBundle:
    """单细胞数据集封装，包含表达值、观测支持掩码、标签等."""
    adata: sc.AnnData
    values: np.ndarray          # 预处理后的表达矩阵（已归一化、log化、[0,1]缩放）
    support: np.ndarray         # 二值掩码：哪些(cell,gene)位置有观测值（非零）
    labels: Optional[np.ndarray]
    label_names: Optional[np.ndarray]
    gene_names: np.ndarray
    x_max: float               # 原始表达最大值，用于反归一化
    input_mode: str


# ---------- PyTorch Dataset ----------

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
    """将稀疏矩阵转为密集numpy数组."""
    if sp.issparse(matrix):
        return matrix.toarray()
    return np.asarray(matrix)


def _infer_labels(adata: sc.AnnData):
    """从AnnData的obs列中查找细胞类型标签（支持多种标准列名）."""
    for candidate in LABEL_CANDIDATES:
        if candidate in adata.obs.columns:
            labels = adata.obs[candidate].astype(str).to_numpy()
            encoder = LabelEncoder()
            encoded = encoder.fit_transform(labels)
            return encoded.astype(np.int64), encoder.classes_, candidate
    return None, None, None


def _looks_like_raw_counts(x: np.ndarray) -> bool:
    """判断是否为原始计数数据：值为非负整数→原始计数；浮点数→已log化."""
    finite = x[np.isfinite(x)]
    if finite.size == 0:
        return False
    finite = finite[: min(finite.size, 100000)]
    return bool(np.all(finite >= 0) and np.allclose(finite, np.round(finite), atol=1e-4))


def _prepare_input_matrix(adata: sc.AnnData, input_mode: str):
    """从AnnData提取表达矩阵：原始计数优先使用raw层."""
    if input_mode == "raw" and adata.raw is not None:
        return _to_dense(adata.raw.X).astype(np.float32), np.asarray(adata.raw.var_names), adata.raw.var.copy()
    return _to_dense(adata.X).astype(np.float32), np.asarray(adata.var_names), adata.var.copy()


# ---------- 主数据加载函数 ----------


def load_sc_dataset(
    file_path: str,
    input_mode: str = "auto",
    n_top_genes: int = 2000,
    normalize_total: bool = True,
    target_sum: float = 1e4,
) -> SCDatasetBundle:
    """加载单细胞h5ad数据，分离表达值与观测支持掩码，并做标准化预处理."""
    adata = sc.read_h5ad(file_path)

    # 步骤1：自动检测数据类型——整数非负→原始计数；浮点数→已log化
    inferred_mode = input_mode
    if inferred_mode == "auto":
        sample = _to_dense(adata.X[: min(256, adata.n_obs)]).astype(np.float32)
        inferred_mode = "raw" if _looks_like_raw_counts(sample) else "log1p"

    # 步骤2：提取表达矩阵（优先使用raw层用于原始计数数据）
    x_input, gene_names, var = _prepare_input_matrix(adata, inferred_mode)
    x_input = np.nan_to_num(x_input, nan=0.0, posinf=0.0, neginf=0.0)
    x_input = np.clip(x_input, 0.0, None)

    # 步骤3：建立二进制支持掩码——表达量>0的位置标记为1（已观测）
    support = (x_input > 0).astype(np.float32)

    # 步骤4：为标准化/HVG筛选构建工作AnnData
    work = sc.AnnData(X=x_input.copy())
    work.obs_names = adata.obs_names.copy()
    work.var_names = gene_names.copy()
    work.obs = adata.obs.copy()
    work.var = var

    # 步骤5：Library size标准化 + log1p变换（原始计数）或条件标准化（已log大数据）
    if inferred_mode == "raw":
        if normalize_total:
            sc.pp.normalize_total(work, target_sum=target_sum)
        sc.pp.log1p(work)
    elif normalize_total:
        values_probe = _to_dense(work.X)
        if values_probe.max(initial=0.0) > 30:
            sc.pp.normalize_total(work, target_sum=target_sum)
            sc.pp.log1p(work)

    # 步骤6：筛选高变异基因（HVG），降维至n_top_genes个基因
    if n_top_genes and n_top_genes > 0 and work.n_vars > n_top_genes:
        sc.pp.highly_variable_genes(work, flavor="seurat", n_top_genes=n_top_genes, subset=True)
        selected = np.asarray(work.var_names)
        gene_to_idx = {gene: idx for idx, gene in enumerate(gene_names)}
        selected_idx = np.array([gene_to_idx[gene] for gene in selected], dtype=np.int64)
        support = support[:, selected_idx]
        gene_names = selected

    # 步骤7：[0,1]归一化：用全局最大值缩放，便于网络训练
    values = _to_dense(work.X).astype(np.float32)
    values = np.nan_to_num(values, nan=0.0, posinf=0.0, neginf=0.0)
    values = np.clip(values, 0.0, None)
    x_max = float(values.max(initial=0.0))
    x_max = max(x_max, 1e-6)
    values = np.clip(values / x_max, 0.0, 1.0).astype(np.float32)

    # 步骤8：从标准obs列提取细胞类型标签
    labels, label_names, label_col = _infer_labels(work)
    if label_col is not None:
        print(f"Detected labels from '{label_col}' with {len(label_names)} classes")
    else:
        print("Warning: no known label column found; evaluation will be skipped")

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


def make_dataloader(
    values: np.ndarray,
    support: np.ndarray,
    labels: Optional[np.ndarray] = None,
    batch_size: int = 256,
    shuffle: bool = True,
    drop_last: bool = False,
) -> DataLoader:
    dataset = FullDataset(values, support, labels)
    return DataLoader(dataset, batch_size=batch_size, shuffle=shuffle, drop_last=drop_last)

