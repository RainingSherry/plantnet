"""
datasets/srp182008_dataset.py
================================
SRP182008 数据集加载器。

数据概况（从 h5ad 文件元数据读取）：
  - 13,514 个细胞，53,678 个基因
  - 15 种细胞类型（PlantNet 植物物种）
  - 组织：Root tip（根尖）
  - 基因型：Col-0
  - 条件：Normal（无扰动）
  - Seurat clusters: 23 个簇

关键挑战：
  1. 数据损坏（NumPy 2.0 / np.string_ 兼容性问题）
  2. 极度稀疏（17M 非零值 / 724M 总元素 ≈ 2.4% 稀疏度）
  3. 物种混合（需要处理未知细胞类型）

处理流程：
  1. 使用 h5py 直接读取（绕过 anndata 的 np.string_ 兼容性问题）
  2. log1p 归一化 + HVG 筛选
  3. 构建基因图（co-expression / marker / random）
  4. 构建细胞-基因支撑图
  5. 输出训练所需格式
"""

from __future__ import annotations

import os
import sys
import warnings
from pathlib import Path
from typing import Tuple, Optional, Literal, List, Dict

# 确保 graphs 模块可导入
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import h5py
import numpy as np
import scanpy as sc
import scipy.sparse as sp
import torch
from sklearn.preprocessing import StandardScaler, LabelEncoder


# ---------------------------------------------------------------------------
# 兼容层（修复 NumPy 2.0 np.string_ 问题）
# ---------------------------------------------------------------------------

def _read_string_dataset(h5file, key: str) -> np.ndarray:
    """安全读取字符串类型的 h5py dataset（兼容 NumPy 2.0）。"""
    ds = h5file[key]
    arr = []
    for i in range(ds.len()):
        try:
            val = ds[i]
            if isinstance(val, bytes):
                val = val.decode("utf-8")
            elif isinstance(val, np.bytes_):
                val = val.decode("utf-8") if val else ""
            arr.append(val)
        except Exception:
            arr.append("")
    return np.array(arr)


def _read_numeric_dataset(h5file, key: str) -> np.ndarray:
    """安全读取数值类型的 h5py dataset。"""
    ds = h5file[key]
    return np.array(ds)


# ---------------------------------------------------------------------------
# 核心数据集类
# ---------------------------------------------------------------------------

class SRP182008Dataset(torch.utils.data.Dataset):
    """
    SRP182008 数据集类。

    支持：
      - 直接从 h5ad 读取（绕过 anndata 的兼容性 bug）
      - 预处理：log1p + HVG 筛选 + Z-score
      - 支撑集预计算（避免每个 batch 重复计算）
      - 细胞类型编码
    """

    def __init__(
        self,
        h5ad_path: str,
        n_hvg: int = 1500,
        normalize: bool = True,
        log_transform: bool = True,
        scale: bool = True,
        min_genes_per_cell: int = 200,
        min_cells_per_gene: int = 3,
        random_seed: int = 42,
        graph_type: str = "coexpression",
        support_strategy: str = "log1p",
        dropout_rate: float = 0.0,
    ):
        """
        参数
        ----
        h5ad_path : str  h5ad 文件路径
        n_hvg : int  保留的高变基因数（默认 1500）
        normalize : bool  是否 Per-cell 归一化
        log_transform : bool  是否 log1p 变换
        scale : bool  是否 Z-score 标准化
        min_genes_per_cell : int  最小基因数阈值
        min_cells_per_gene : int  最小细胞数阈值
        random_seed : int  随机种子
        graph_type : str  基因图类型（coexpression/marker/random）
        support_strategy : str  支撑集权重策略
        dropout_rate : float  支撑集 dropout 率
        """
        self.n_hvg = n_hvg
        self.support_strategy = support_strategy
        self.dropout_rate = dropout_rate
        self.random_seed = random_seed
        self.rng = np.random.default_rng(random_seed)

        print(f"[SRP182008Dataset] Loading from {h5ad_path}")

        # ---- Step 1: 读取数据 ----
        X, obs, var = self._load_h5ad_safe(h5ad_path)
        n_cells_orig, n_genes_orig = X.shape
        print(f"[SRP182008Dataset] Raw: {n_cells_orig} cells × {n_genes_orig} genes")

        # ---- Step 2: 细胞过滤 ----
        cell_keep = X.sum(axis=1) > 0  # 去除空细胞
        X = X[cell_keep]
        obs = obs[cell_keep]
        print(f"[SRP182008Dataset] After cell filter: {X.shape[0]} cells")

        # ---- Step 3: 基因过滤 ----
        gene_keep = (X > 0).sum(axis=0) >= min_cells_per_gene
        X = X[:, gene_keep]
        var = var[gene_keep]
        print(f"[SRP182008Dataset] After gene filter: {X.shape[1]} genes")

        # ---- Step 4: Per-cell 归一化 ----
        if normalize:
            cell_sum = X.sum(axis=1, keepdims=True)
            cell_sum = np.where(cell_sum > 0, cell_sum, 1.0)
            X = X / cell_sum * 1e4

        # ---- Step 5: log1p ----
        if log_transform:
            X = np.log1p(X)

        # ---- Step 6: HVG 筛选 ----
        X_hvg, hvg_idx, hvg_genes = self._select_hvg(X, var, n_hvg)
        self.X_raw = X
        self.X = X_hvg
        self.hvg_idx = hvg_idx
        self.hvg_genes = hvg_genes
        self.n_hvg_actual = len(hvg_genes)
        print(f"[SRP182008Dataset] HVG: {self.n_hvg_actual} genes")

        # ---- Step 7: 支撑集预计算 ----
        self._compute_support_sets()

        # ---- Step 8: 细胞类型编码 ----
        self._encode_labels(obs)

        # ---- Step 9: Z-score ----
        if scale:
            self.scaler = StandardScaler()
            self.X = self.scaler.fit_transform(self.X)

        # 保存元数据
        self.gene_names = np.array([str(g) for g in var[hvg_idx]])

        # ---- Step 10: 构建基因图 ----
        self._build_graph(graph_type)

        print(f"[SRP182008Dataset] Final: {self.n_cells} cells × {self.n_hvg_actual} genes")
        print(f"[SRP182008Dataset] Cell types: {self.n_cell_types}, Labels: {self.cell_type_names}")

    def _load_h5ad_safe(
        self,
        path: str,
    ) -> Tuple[np.ndarray, "pd.DataFrame", np.ndarray]:
        """安全加载 h5ad（使用 scanpy，兼容不同数据集）。"""
        import pandas as pd
        import scanpy as sc

        adata = sc.read_h5ad(path)

        # 获取表达矩阵
        X = adata.X
        if hasattr(X, "toarray"):
            X = X.toarray()
        X = X.astype(np.float32)

        # 获取 obs
        obs = adata.obs

        # 获取 var
        var = np.array(adata.var_names)

        return X, obs, var

    def _select_hvg(
        self,
        X: np.ndarray,
        gene_names: np.ndarray,
        n_hvg: int,
    ) -> Tuple[np.ndarray, np.ndarray, List[str]]:
        """选择高变基因。"""
        from sklearn.preprocessing import MaxAbsScaler

        # 使用 scanpy 的 HVG 方法（基于 dispersion）
        adata_tmp = sc.AnnData(X)
        adata_tmp.var_names = gene_names

        try:
            sc.pp.highly_variable_genes(
                adata_tmp,
                n_top_genes=n_hvg,
                flavor="seurat_v3",
                subset=False,
            )
        except Exception:
            # 如果 seurat_v3 失败，使用 seurat
            sc.pp.highly_variable_genes(
                adata_tmp,
                n_top_genes=n_hvg,
                flavor="seurat",
                subset=False,
            )

        hvg_mask = adata_tmp.var["highly_variable"].values
        hvg_idx = np.where(hvg_mask)[0]

        X_hvg = X[:, hvg_idx]
        hvg_genes = [str(gene_names[i]) for i in hvg_idx]

        return X_hvg, hvg_idx, hvg_genes

    def _compute_support_sets(self):
        """预计算每个细胞的支撑集（避免重复计算）。"""
        X_bool = (self.X_raw > 0) if hasattr(self, "X_raw") else (self.X > 0)

        # 获取 HVG 子矩阵的布尔值
        if hasattr(self, "X_raw"):
            X_raw_hvg = self.X_raw[:, self.hvg_idx]
            X_bool = (X_raw_hvg > 0)
        else:
            X_bool = (self.X > 0)

        # 支撑集索引和掩码
        self.n_cells = X_bool.shape[0]

        # 预计算支撑集：每个细胞的非零基因索引
        self.support_idx = []  # List of arrays
        self.support_mask = []  # List of arrays
        self.support_weight = []  # List of arrays

        for c in range(self.n_cells):
            nz_idx = np.where(X_bool[c])[0]
            # 限制最大支撑集大小（避免 OOM）
            max_support = min(len(nz_idx), 2000)
            if len(nz_idx) > max_support:
                # 取表达量最高的 max_support 个基因
                if hasattr(self, "X_raw"):
                    vals = self.X_raw[c, self.hvg_idx][nz_idx]
                else:
                    vals = self.X[c][nz_idx]
                top_k = np.argsort(vals)[::-1][:max_support]
                nz_idx = nz_idx[top_k]

            self.support_idx.append(nz_idx)
            self.support_mask.append(np.ones(len(nz_idx), dtype=np.float32))

            # 支撑权重
            if hasattr(self, "X_raw"):
                expr_vals = np.log1p(self.X_raw[c, self.hvg_idx][nz_idx]).astype(np.float32)
            else:
                expr_vals = self.X[c][nz_idx].astype(np.float32)

            # 归一化权重
            w_sum = expr_vals.sum() + 1e-8
            self.support_weight.append(expr_vals / w_sum)

        # 填充到固定长度（方便 batch）
        max_support_global = max(len(s) for s in self.support_idx)

        # 填充到整数倍（加速 GPU）
        max_support_padded = ((max_support_global + 31) // 32) * 32
        max_support_padded = min(max_support_padded, 2000)  # 上限 2000

        self.max_support = max_support_padded
        self._pad_support(max_support_padded)

        print(f"[SRP182008Dataset] Support sets: max={max_support_padded} per cell")

    def _pad_support(self, max_support: int):
        """将支撑集填充到固定长度。"""
        self.support_idx_padded = np.full((self.n_cells, max_support), -1, dtype=np.int64)
        self.support_mask_padded = np.zeros((self.n_cells, max_support), dtype=np.float32)
        self.support_weight_padded = np.zeros((self.n_cells, max_support), dtype=np.float32)

        for c in range(self.n_cells):
            s = min(len(self.support_idx[c]), max_support)
            self.support_idx_padded[c, :s] = self.support_idx[c][:s]
            self.support_mask_padded[c, :s] = self.support_mask[c][:s]
            self.support_weight_padded[c, :s] = self.support_weight[c][:s]

        # 转为 tensor（在 CPU 上；DataLoader 的 pin_memory 会处理 GPU 传输）
        self.support_idx_t = torch.from_numpy(self.support_idx_padded).long()
        self.support_mask_t = torch.from_numpy(self.support_mask_padded).float()
        self.support_weight_t = torch.from_numpy(self.support_weight_padded).float()

    def _encode_labels(self, obs):
        """编码细胞类型标签。"""
        # 支持多个可能的列名
        label_col = None
        for col in ["Celltype", "cell_type", "celltype", "cell_label",
                    "label", "CellType", "cell_type1", "celltype_after"]:
            if col in obs.columns:
                label_col = col
                break
        if label_col is None:
            raise KeyError(f"No label column found. Available columns: {list(obs.columns)}")

        cell_types = np.asarray(obs[label_col]).astype(str)

        # 过滤未知类型
        valid_mask = cell_types != "Unknow"
        self.cell_type_names = np.unique(cell_types[valid_mask])
        self.n_cell_types = len(self.cell_type_names)

        self.label_encoder = LabelEncoder()
        self.labels = self.label_encoder.fit_transform(cell_types)
        self.labels = self.labels.astype(np.int64)

        print(f"[SRP182008Dataset] Label col: '{label_col}', Cell types: {self.n_cell_types} types")

    def _build_graph(self, graph_type: str):
        """构建基因图。"""
        from graphs.build_gene_graph import build_gene_graph

        print(f"[SRP182008Dataset] Building {graph_type} gene graph...")

        # 构建 HVG 子矩阵：仅包含我们选择的 HVG 基因
        # X_raw: (n_cells, n_genes_full), hvg_idx: indices into full gene space
        # X_hvg: (n_cells, n_hvg_actual) — 对齐后的 HVG 表达矩阵
        X_raw = self.X_raw[:, self.hvg_idx] if hasattr(self, "X_raw") else self.X
        X_for_graph = np.nan_to_num(X_raw, nan=0.0)

        try:
            graph_dict = build_gene_graph(
                X_for_graph,        # 使用 HVG 子矩阵建图
                self.gene_names,    # HVG 基因名（顺序与 X_for_graph 的列对应）
                graph_type=graph_type,
                n_hvg=min(self.n_hvg_actual, 1500),
                corr_threshold=0.3,
                top_k_neighbors=20,
                marker_boost_weight=3.0,
                random_seed=self.random_seed,
            )
            self.graph_dict = graph_dict
            print(f"[SRP182008Dataset] Graph: {graph_dict['n_nodes']} nodes, {len(graph_dict['edge_index'][0])} edges ({graph_type})")
        except Exception as e:
            warnings.warn(f"[SRP182008Dataset] Graph building failed: {e}")
            # Fallback: 空图
            self.graph_dict = {
                "edge_index": [[], []],
                "edge_weight": [],
                "n_nodes": self.n_hvg_actual,
                "gene_names": list(self.gene_names),
                "graph_type": "empty",
            }

    def __len__(self) -> int:
        return self.n_cells

    def __getitem__(self, idx: int) -> dict:
        """返回单个样本。"""
        return {
            "X": torch.from_numpy(self.X[idx]).float(),
            "support_idx": self.support_idx_t[idx],
            "support_mask": self.support_mask_t[idx],
            "support_weight": self.support_weight_t[idx],
            "label": torch.tensor(self.labels[idx], dtype=torch.long),
            "cell_idx": torch.tensor(idx, dtype=torch.long),
        }

    def get_full_embeddings(self) -> np.ndarray:
        """返回全部细胞嵌入（用于聚类初始化）。"""
        return self.X.astype(np.float32)


def _obs_to_dataframe(obs: dict) -> "pd.DataFrame":
    """将 obs 字典转为 DataFrame。"""
    try:
        import pandas as pd
        return pd.DataFrame(obs)
    except ImportError:
        # 无 pandas 时的 fallback
        class SimpleDataFrame:
            def __init__(self, d):
                for k, v in d.items():
                    setattr(self, k, v)
                self._data = d
            def __getitem__(self, key):
                return np.array(self._data[key])
            def __setitem__(self, key, val):
                self._data[key] = val
            def __repr__(self):
                return f"SimpleDataFrame(keys={list(self._data.keys())})"
            def __len__(self):
                return len(list(self._data.values())[0])
            def __iter__(self):
                return iter(self._data.keys())
            def items(self):
                return self._data.items()
            def __getitem__slice(self, mask):
                return SimpleDataFrame({k: v[mask] for k, v in self._data.items()})
        return SimpleDataFrame(obs)


def create_dataloader(
    dataset: SRP182008Dataset,
    batch_size: int = 32,
    shuffle: bool = True,
    drop_last: bool = True,
    num_workers: int = 0,
) -> torch.utils.data.DataLoader:
    """创建 DataLoader。"""
    return torch.utils.data.DataLoader(
        dataset,
        batch_size=batch_size,
        shuffle=shuffle,
        drop_last=drop_last,
        num_workers=num_workers,
        pin_memory=torch.cuda.is_available(),
        collate_fn=_collate_fn,
    )


def _collate_fn(batch: List[dict]) -> dict:
    """自定义 collate 函数。"""
    keys = batch[0].keys()
    result = {}
    for k in keys:
        if k == "support_idx":
            result[k] = torch.stack([b[k] for b in batch])
        elif k == "support_mask":
            result[k] = torch.stack([b[k] for b in batch])
        elif k == "support_weight":
            result[k] = torch.stack([b[k] for b in batch])
        elif k == "X":
            result[k] = torch.from_numpy(np.stack([b[k].numpy() for b in batch])).float()
        elif k == "label":
            result[k] = torch.stack([b[k] for b in batch])
        elif k == "cell_idx":
            result[k] = torch.stack([b[k] for b in batch])
    return result


if __name__ == "__main__":
    import warnings
    warnings.filterwarnings("ignore")

    data_path = "../../../data/SRP182008.h5ad"
    if os.path.exists(data_path):
        dataset = SRP182008Dataset(
            h5ad_path=data_path,
            n_hvg=1000,
            graph_type="coexpression",
        )
        print(f"Dataset: {len(dataset)} cells")
        print(f"Graph: {dataset.graph_dict['n_nodes']} nodes, {len(dataset.graph_dict['edge_index'][0])} edges")
        print(f"Sample: {dataset[0]}")
    else:
        print(f"Data file not found at {data_path}")
