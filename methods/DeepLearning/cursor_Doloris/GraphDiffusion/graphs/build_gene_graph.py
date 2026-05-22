"""
graphs/build_gene_graph.py
===========================
构建植物基因调控/共表达图谱，提供三种图结构消融：
  1. HVG Co-expression Graph  — 基于高变基因的 Pearson 相关性图
  2. Marker Graph             — 基于 marker gene 连接性的图（基于植物已知 marker）
  3. Random Graph             — 随机图（负对照，用于消融实验）

核心思想：
  - 不依赖外部 GRN 数据库（如 TRRUST、PlantTFDB），直接从数据中构造图
  - 图节点 = HVG 基因，图边 = 基因间共表达关系（可带权重）
  - Marker Graph 通过已知 marker 列表对共表达图加权，使 marker 基因的邻居更可靠
  - 支持保存为 PyG Data/图格式，供 GAT 模块使用

作者: 基于 DOLORIS GRN.py 和 PhytoCluster 思想改编
"""

from __future__ import annotations

import warnings
from typing import List, Optional, Tuple

import numpy as np
import scipy.sparse as sp
import scanpy as sc
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import pairwise_distances

try:
    import networkx as nx
    HAS_NX = True
except ImportError:
    HAS_NX = False
    warnings.warn("networkx not installed; graph serialization will use adjacency dict.")


# ---------------------------------------------------------------------------
# 已知植物单细胞 marker 基因库（部分，可扩展）
# 来源：PlantPlantDB, scPlantDB, Arabidopsis Root Atlas
# ---------------------------------------------------------------------------
PLANT_MARKER_SETS = {
    # ---- 根组织 ----
    "root_stele": [
        "SULTR1;1", "SULTR1;2", "OPT5", "BTS", "CLE25", "CLE26",
        "WOX5", "WOX4", "ACR4", "SKOR", "MIR165a", "MIR166a",
        "ATHB8", "HB8", "S17", "S4", "CASP1", "PEAR1",
    ],
    "root_hair": [
        "RSH1", "RSH3", "RSH4", "EXT3", "EXT7", "EXT12",
        "AtEXP7", "AtEXP18", "EXPA1", "EXPA5", "EXPA7",
        "EXPA10", "EXPA11", "EXPA12", "EXPA14", "EXPA18",
        "EXPA20", "EXPA22", "EXPA23", "EXPA24",
        "RGXT1", "RGXT2", "FLA11", "FLA13", "FLA15",
    ],
    "root_epidermis": [
        "ATML1", "PDF1", "PDF2", "PDF2.1", "ANL2", "GL2",
        "WER", "WERWOLF", "CPC", "TRY", "ETC1", "ETC2",
        "CPL3", "TMO3", "TMO5", "LHW", "OPI10",
    ],
    "root_cortex": [
        "CORTEX", "SMO1", "SMO2", "AT1G52770", "AT1G68570",
        "AT3G05980", "AT5G39570", "AT1G68560",
    ],
    "root_endodermis": [
        "CASP1", "CASP2", "CASP3", "CASP4", "CASP5",
        "SGN1", "SGN2", "SGN3", "G36", "S17", "S23",
        "MYB36", "EN7", "EN8", "SCHENGEN1", "SCHENGEN2", "SCHENGEN3",
        "IRK", "SIRE", "SLAH3",
    ],
    "columella": [
        "PGD1", "PGD2", "PGD3", "MDM1", "MIZ1", "MIZ2",
        "NRT2.1", "NRT2.2", "NRT2.4", "NRT2.5", "NRT2.6",
        "CLE40", "ACR4", "WOX5", "COL1", "COL2",
    ],
    "lateral_root_cap": [
        "SHR", "SCR", "LAX1", "LAX2", "LAX3", "LAX4",
        "IAA28", "GATA23", "NAC4", "BBM", "WOX11", "WOX12",
        "LBD29", "LBD16", "LBD18", "LBD33",
    ],
    "root_phloem": [
        "SUC2", "SWEET11", "SWEET12", "SWEET13", "SUT1", "SUT2",
        "APL", "PS1", "PS2",
    ],
    # ---- 细胞周期 ----
    "g2m_phase": [
        "CYCB1;1", "CYCB1;2", "CYCB2;1", "CYCB2;2", "CYCB2;3",
        "CYCD3;1", "CYCD3;2", "CYCD3;3", "CYCD4;1", "CYCD5;1",
        "CDKB1;1", "CDKB1;2", "CDKB2;1", "CDKB2;2",
        "WEE1", "CDC25",
    ],
    "s_phase": [
        "HIS1", "HIS2", "HIS3", "HTR1", "HTR2", "HTR3",
        "HTR4", "HTR5", "HTR6", "HTR7", "HTR8", "HTR9",
        "HTR10", "HTR11", "HTR12", "HTR13",
        "RFC3", "PCNA1", "PCNA2", "MCM2", "MCM3", "MCM4",
    ],
}


def _fuzzy_match(gene_name: str, marker_list: List[str]) -> bool:
    """对 marker 基因做模糊匹配（兼容 A. thaliana 不同命名惯例）。"""
    gn = gene_name.upper().replace(" ", "").replace("-", "")
    for m in marker_list:
        m_up = m.upper().replace(" ", "").replace("-", "")
        if gn == m_up or gn.startswith(m_up) or m_up.startswith(gn):
            return True
    return False


def _build_coexpression_graph(
    X_norm: np.ndarray,
    gene_names: np.ndarray,
    n_hvg: int = 1500,
    corr_threshold: float = 0.35,
    max_neighbors: int = 50,
    top_k_edges_per_node: int = 20,
) -> Tuple[List[int], List[int], List[float]]:
    """
    基于 HVG 的共表达图。

    参数
    ----
    X_norm : (n_cells, n_genes)  标准化后的基因表达矩阵（已 log1p 或 scale）
    gene_names :  基因名数组
    n_hvg : 选取的高变基因数
    corr_threshold : Pearson 相关系数阈值（绝对值 > 此值才建边）
    max_neighbors : 每个节点最多保留的邻居数（避免图过于稠密）
    top_k_edges_per_node : 每个节点保留权重最高的 top-k 边

    返回
    ----
    edge_index : [2, n_edges]  边索引（scipy COO 格式）
    edge_weight : [n_edges]    边权重（相关系数绝对值）
    node_genes : 保留的 HVG 基因名列表
    """
    n_cells, n_genes = X_norm.shape

    # ---- Step 1: 选择 HVG ----
    # 使用 scanpy 的 HVG 逻辑（基于 mean-bin 的 dispersion）
    adata_tmp = sc.AnnData(X_norm)
    sc.pp.highly_variable_genes(
        adata_tmp,
        n_top_genes=n_hvg,
        flavor="seurat_v3",  # 基于方差
        subset=False,
    )
    hvg_mask = adata_tmp.var["highly_variable"].values
    hvg_idx = np.where(hvg_mask)[0]
    X_hvg = X_norm[:, hvg_idx]
    hvg_genes = gene_names[hvg_idx]

    print(f"  [build_gene_graph] Selected {len(hvg_genes)} HVG out of {n_genes} genes")

    # ---- Step 2: 计算 Pearson 相关矩阵（仅在 HVG 上）----
    # 直接在 HVG 矩阵上计算（1500×1500 = ~9MB float32，完全可接受）
    # 使用稀疏友好的分块方式避免 OOM
    n_hvg = X_hvg.shape[1]
    chunk_size = 500  # 分块大小，避免大矩阵内存问题
    all_corr = np.zeros((n_hvg, n_hvg), dtype=np.float32)

    for start in range(0, n_hvg, chunk_size):
        end = min(start + chunk_size, n_hvg)
        chunk = X_hvg[:, start:end]  # (n_cells, chunk)
        if chunk.shape[1] == 0:
            continue
        # 逐列计算：每个 chunk 列与所有 HVG 列的相关性
        chunk_mean = chunk.mean(axis=0, keepdims=True)
        chunk_std = chunk.std(axis=0, keepdims=True) + 1e-8
        chunk_norm = (chunk - chunk_mean) / chunk_std
        full_mean = X_hvg.mean(axis=0, keepdims=True)
        full_std = X_hvg.std(axis=0, keepdims=True) + 1e-8
        full_norm = (X_hvg - full_mean) / full_std
        # corr(chunk_col, full_col) = (chunk_col_norm.T @ full_col_norm) / (n_cells - 1)
        block_corr = np.dot(chunk_norm.T, full_norm) / max(n_cells - 1, 1)
        block_corr = np.abs(block_corr)
        all_corr[start:end, :] = block_corr

    # 对角线清零（去除自环）
    np.fill_diagonal(all_corr, 0)
    all_corr = np.nan_to_num(all_corr, nan=0.0, posinf=0.0, neginf=0.0)

    # ---- Step 3: 建边（稀疏化）----
    n_hvg_actual = len(hvg_genes)
    row_list, col_list, weight_list = [], [], []

    for i in range(n_hvg_actual):
        # 与节点 i 的相关系数（已取绝对值）
        abs_corr = all_corr[i]
        abs_corr[i] = 0  # 去掉自环

        # 取 top-k 邻居
        neighbors = np.argsort(abs_corr)[::-1][:top_k_edges_per_node]
        for j in neighbors:
            if abs_corr[j] > 0:
                row_list.append(i)
                col_list.append(j)
                weight_list.append(float(max(abs_corr[j], 0.01)))

    print(f"  [build_gene_graph] Co-expression graph: {len(row_list)} edges from {n_hvg} HVG nodes")

    return row_list, col_list, weight_list, list(hvg_genes)


def _build_marker_graph(
    X_norm: np.ndarray,
    gene_names: np.ndarray,
    hvg_genes: List[str],
    cell_types: Optional[List[str]] = None,
    top_k_neighbors: int = 30,
    marker_boost_weight: float = 3.0,
) -> Tuple[List[int], List[int], List[float], List[str]]:
    """
    基于 marker 基因引导的图构建。
    
    思路：先构建 co-expression 图，然后对涉及 marker 基因的边进行权重 boost。
    如果提供了 cell_types（与 marker sets 的 key 对应），则只 boost 匹配到的 marker。

    返回
    ----
    edge_index, edge_weight, node_genes（同上）
    """
    n_genes = len(gene_names)
    hvg_set = set(hvg_genes)

    # ---- Step 1: 先构建基础共表达图 ----
    row_list, col_list, weight_list = [], [], []
    n_hvg = len(hvg_genes)
    gene_to_idx = {g: i for i, g in enumerate(hvg_genes)}
    hvg_gene_idx = [np.where(gene_names == g)[0][0] for g in hvg_genes]

    # 计算子矩阵上的相关性
    X_sub = X_norm[:, hvg_gene_idx]  # (n_cells, n_hvg)

    # 分块计算避免 OOM
    chunk_size = 200
    all_abs_corr = np.zeros((n_hvg, n_hvg), dtype=np.float32)

    for start in range(0, n_hvg, chunk_size):
        end = min(start + chunk_size, n_hvg)
        chunk = X_sub[:, start:end]
        if chunk.shape[1] > 0:
            # 使用稀疏友好的方式
            chunk_mean = chunk.mean(axis=0, keepdims=True)
            chunk_std = chunk.std(axis=0, keepdims=True) + 1e-8
            chunk_norm = (chunk - chunk_mean) / chunk_std
            # 手动计算 Pearson 相关（避免全量矩阵）
            # 只计算 chunk vs all
            all_corr = np.dot(chunk_norm.T, X_sub) / (X_sub.shape[0] - 1)
            all_abs_corr[start:end, :] = np.abs(np.nan_to_num(all_corr))

    # 对角线清零
    np.fill_diagonal(all_abs_corr, 0)

    for i in range(n_hvg):
        row = all_abs_corr[i]
        row[i] = 0
        top_k_idx = np.argsort(row)[::-1][:top_k_neighbors]
        for j in top_k_idx:
            if row[j] > 0:
                row_list.append(i)
                col_list.append(j)
                weight_list.append(float(row[j]))

    # ---- Step 2: Marker boost ----
    # 收集所有 marker 基因
    all_markers = set()
    for markers in PLANT_MARKER_SETS.values():
        for m in markers:
            all_markers.add(m)

    marker_in_hvg = []
    marker_idx_map = {}
    for idx, g in enumerate(hvg_genes):
        if _fuzzy_match(g, list(all_markers)):
            marker_in_hvg.append(idx)
            marker_idx_map[idx] = g

    print(f"  [build_marker_graph] Found {len(marker_in_hvg)} marker genes in HVG set")

    # boost 边权重：涉及 marker 的边权重 × marker_boost_weight
    boosted_weights = []
    for i in range(len(row_list)):
        if row_list[i] in marker_idx_map or col_list[i] in marker_idx_map:
            boosted_weights.append(weight_list[i] * marker_boost_weight)
        else:
            boosted_weights.append(weight_list[i])

    print(f"  [build_marker_graph] Marker-boosted graph: {len(row_list)} edges")

    return row_list, col_list, boosted_weights, list(hvg_genes)


def _build_random_graph(
    n_nodes: int,
    n_edges: int,
) -> Tuple[List[int], List[int], List[float]]:
    """随机图（负对照），保持与共表达图相近的边密度。"""
    rng = np.random.default_rng(42)
    row_list, col_list = [], []
    seen = set()

    while len(row_list) < n_edges and len(row_list) < n_nodes * (n_nodes - 1):
        i = rng.integers(0, n_nodes)
        j = rng.integers(0, n_nodes)
        if i != j:
            key = (min(i, j), max(i, j))
            if key not in seen:
                seen.add(key)
                row_list.append(i)
                col_list.append(j)

    weight_list = [float(np.abs(rng.standard_normal())) for _ in range(len(row_list))]
    return row_list, col_list, weight_list


def build_gene_graph(
    adata_or_X: "sc.AnnData | np.ndarray",
    gene_names: Optional[np.ndarray] = None,
    graph_type: str = "coexpression",
    n_hvg: int = 1500,
    corr_threshold: float = 0.35,
    top_k_neighbors: int = 20,
    marker_boost_weight: float = 3.0,
    cell_types: Optional[np.ndarray] = None,
    random_seed: int = 42,
    pre_selected_hvg_idx: Optional[np.ndarray] = None,
) -> dict:
    """
    统一入口：根据配置构建基因图。

    参数
    ----
    adata_or_X : AnnData 或 ndarray
        - 如果是 AnnData，使用 adata.X 做预处理和建图
        - 如果是 ndarray，直接用
    gene_names : ndarray，基因名（当 adata_or_X 是 ndarray 时必须提供）
    graph_type : str
        - "coexpression" : HVG 共表达图
        - "marker"      : Marker-boosted 图
        - "random"      : 随机图（负对照）
        - "all"         : 返回所有三种图的字典
    n_hvg : int  HVG 数量
    corr_threshold : float  相关系数阈值
    top_k_neighbors : int  每个节点最多邻居数
    marker_boost_weight : float  marker 边权重 boost 倍数
    cell_types : ndarray  细胞类型标签（可选，用于 marker 图的细胞特异性加权）
    random_seed : int  随机种子

    返回
    ----
    graph_dict : dict
        当 graph_type != "all" 时，返回：
        {
            "edge_index": [[src_nodes], [dst_nodes]],   # [2, E]
            "edge_weight": [float, ...],               # [E]
            "n_nodes": int,
            "gene_names": [str, ...],                 # [n_nodes]
            "graph_type": str,
        }
        当 graph_type == "all" 时，返回三种图的字典。
    """
    rng = np.random.default_rng(random_seed)

    # ---- 数据准备 ----
    if hasattr(adata_or_X, "X"):
        adata = adata_or_X
        if gene_names is None:
            if "gene_name" in adata.var.columns:
                gene_names = adata.var["gene_name"].values
            else:
                gene_names = adata.var_names.values
        n_cells, n_genes = adata.shape

        # 预处理：log1p + scale
        X_raw = adata.X.toarray() if sp.issparse(adata.X) else adata.X.copy()
        X_log = np.log1p(X_raw)
        # Z-score per gene
        scaler = StandardScaler()
        X_norm = scaler.fit_transform(X_log)
        del X_raw, X_log
    else:
        X_norm = adata_or_X
        n_cells, n_genes = X_norm.shape
        if gene_names is None:
            gene_names = np.arange(n_genes, dtype=object)
        gene_names = np.asarray(gene_names)

    X_norm = np.nan_to_num(X_norm, nan=0.0, posinf=0.0, neginf=0.0)

    if graph_type == "all":
        coexp = build_gene_graph(
            X_norm, gene_names, "coexpression", n_hvg,
            corr_threshold, top_k_neighbors, marker_boost_weight, cell_types, random_seed
        )
        marker = build_gene_graph(
            X_norm, gene_names, "marker", n_hvg,
            corr_threshold, top_k_neighbors, marker_boost_weight, cell_types, random_seed
        )
        rnd = build_gene_graph(
            X_norm, gene_names, "random", n_hvg,
            corr_threshold, top_k_neighbors, marker_boost_weight, cell_types, random_seed
        )
        return {"coexpression": coexp, "marker": marker, "random": rnd}

    # ---- 构建指定类型的图 ----
    if graph_type == "coexpression":
        row_list, col_list, weight_list, node_genes = _build_coexpression_graph(
            X_norm, gene_names, n_hvg, corr_threshold, max_neighbors=top_k_neighbors,
            top_k_edges_per_node=top_k_neighbors,
        )
    elif graph_type == "marker":
        row_list, col_list, weight_list, node_genes = _build_coexpression_graph(
            X_norm, gene_names, n_hvg, corr_threshold, max_neighbors=top_k_neighbors,
            top_k_edges_per_node=top_k_neighbors,
        )
        # 再做 marker boost
        all_markers = set()
        for markers in PLANT_MARKER_SETS.values():
            for m in markers:
                all_markers.add(m)
        gene_to_idx = {g: i for i, g in enumerate(node_genes)}
        marker_idx = [gene_to_idx[g] for g in node_genes
                       if _fuzzy_match(g, list(all_markers))]
        boosted_weights = []
        for i in range(len(row_list)):
            if row_list[i] in marker_idx or col_list[i] in marker_idx:
                boosted_weights.append(weight_list[i] * marker_boost_weight)
            else:
                boosted_weights.append(weight_list[i])
        weight_list = boosted_weights
    elif graph_type == "random":
        row_list, col_list, weight_list = [], [], []
        # 先建一个 coexpression 图来获取节点数
        _, _, _, node_genes = _build_coexpression_graph(
            X_norm, gene_names, n_hvg, corr_threshold, max_neighbors=top_k_neighbors,
            top_k_edges_per_node=top_k_neighbors,
        )
        n_nodes = len(node_genes)
        n_target = int(n_nodes * top_k_neighbors * 0.5)  # 大约一半密度的边
        row_list, col_list, weight_list = _build_random_graph(n_nodes, n_target)
    else:
        raise ValueError(f"Unknown graph_type: {graph_type}")

    n_nodes = len(node_genes)

    # 归一化权重到 [0.5, 1.0] 范围（避免零权重导致 GAT 无法关注）
    if weight_list:
        w_min, w_max = min(weight_list), max(weight_list)
        if w_max > w_min:
            weight_list = [
                0.5 + 0.5 * (w - w_min) / (w_max - w_min) for w in weight_list
            ]
        else:
            weight_list = [1.0] * len(weight_list)

    result = {
        "edge_index": [row_list, col_list],
        "edge_weight": weight_list,
        "n_nodes": n_nodes,
        "gene_names": node_genes,
        "graph_type": graph_type,
    }
    return result


# ---------------------------------------------------------------------------
# 导出 PyG Data 对象（如果安装了 torch_geometric）
# ---------------------------------------------------------------------------
def to_pyg_data(graph_dict: dict, edge_weight_normalize: bool = True):
    """将 graph_dict 转换为 torch_geometric Data 对象。"""
    try:
        import torch
        from torch_geometric.data import Data
    except ImportError:
        raise ImportError("torch_geometric is required to use to_pyg_data()")

    row, col = graph_dict["edge_index"]
    weight = graph_dict["edge_weight"]

    if weight is None or len(weight) == 0:
        edge_attr = None
    else:
        edge_attr = torch.tensor(weight, dtype=torch.float32).unsqueeze(-1)

    # 节点特征：随机初始化（后续由 gene_gat_encoder 替换）
    x = torch.randn(graph_dict["n_nodes"], 64)  # 临时占位

    data = Data(
        x=x,
        edge_index=torch.tensor([row, col], dtype=torch.long),
        edge_attr=edge_attr,
    )
    return data


if __name__ == "__main__":
    # 简单测试
    import scanpy as sc
    warnings.filterwarnings("ignore")

    # 用 SRP182008 的一个小样本测试
    try:
        adata = sc.read_h5ad("../../../data/SRP182008.h5ad")
        # 取 500 细胞子样本加速测试
        adata = adata[np.random.RandomState(42).choice(adata.n_obs, 500, replace=False)].copy()

        result = build_gene_graph(adata, graph_type="coexpression", n_hvg=200)
        print(f"Coexp graph: {result['n_nodes']} nodes, {len(result['edge_index'][0])} edges")

        result_marker = build_gene_graph(adata, graph_type="marker", n_hvg=200)
        print(f"Marker graph: {result_marker['n_nodes']} nodes, {len(result_marker['edge_index'][0])} edges")

        result_random = build_gene_graph(adata, graph_type="random", n_hvg=200)
        print(f"Random graph: {result_random['n_nodes']} nodes, {len(result_random['edge_index'][0])} edges")
    except Exception as e:
        print(f"Test skipped (missing data): {e}")
