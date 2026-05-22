import json
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import scanpy as sc
import scipy.sparse as sp
from sklearn.decomposition import TruncatedSVD
from sklearn.feature_extraction.text import TfidfTransformer
from sklearn.metrics.pairwise import cosine_similarity


@dataclass
class GeneGraphBundle:
    gene_names: List[str]
    edge_index: np.ndarray
    edge_weight: np.ndarray
    edge_type: np.ndarray
    gene_features: np.ndarray
    marker_modules: Dict[str, List[int]]
    metadata: Dict[str, object]


EDGE_TYPES = {
    "coexpression": 0,
    "marker": 1,
    "orthology": 2,
}


def _ensure_numpy_compat() -> None:
    if not hasattr(np, "string_"):
        np.string_ = np.bytes_
    if not hasattr(np, "unicode_"):
        np.unicode_ = np.str_


def _infer_label_col(adata: sc.AnnData, label_col: Optional[str]) -> str:
    if label_col and label_col in adata.obs.columns:
        return label_col
    for candidate in ["Celltype", "celltype", "cell_type", "celltype_after", "cell_label", "label", "Seurat_clusters"]:
        if candidate in adata.obs.columns:
            return candidate
    raise KeyError(f"No label column found. Available columns: {list(adata.obs.columns)}")


def load_adata(data_path: str) -> sc.AnnData:
    _ensure_numpy_compat()
    return sc.read_h5ad(data_path)


def select_hvgs(
    adata: sc.AnnData,
    n_top_genes: int = 2000,
    target_sum: float = 1e4,
) -> sc.AnnData:
    adata = adata.copy()
    if adata.raw is not None:
        adata.X = adata.raw.X.copy()
    sc.pp.filter_genes(adata, min_counts=1)
    sc.pp.filter_cells(adata, min_counts=1)
    sc.pp.normalize_total(adata, target_sum=target_sum)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=n_top_genes, subset=True)
    return adata


def build_gene_features(
    adata: sc.AnnData,
    feature_dim: int = 64,
    seed: int = 42,
) -> np.ndarray:
    x = adata.X
    if sp.issparse(x):
        x = x.tocsr()
    transformer = TfidfTransformer(norm="l2", use_idf=True, smooth_idf=True, sublinear_tf=True)
    tfidf = transformer.fit_transform(x)
    max_components = min(feature_dim, tfidf.shape[0] - 1, tfidf.shape[1] - 1)
    if max_components < 2:
        dense = tfidf.toarray() if sp.issparse(tfidf) else np.asarray(tfidf)
        return dense.astype(np.float32)
    svd = TruncatedSVD(n_components=max_components, random_state=seed)
    gene_features = svd.fit_transform(tfidf.T)
    if gene_features.shape[1] < feature_dim:
        pad = np.zeros((gene_features.shape[0], feature_dim - gene_features.shape[1]), dtype=np.float32)
        gene_features = np.concatenate([gene_features.astype(np.float32), pad], axis=1)
    return gene_features.astype(np.float32)


def build_coexpression_edges(
    adata: sc.AnnData,
    corr_threshold: float = 0.35,
    top_k: int = 20,
) -> List[Tuple[int, int, float, int]]:
    x = adata.X
    if sp.issparse(x):
        x = x.toarray()
    x = np.asarray(x, dtype=np.float32)
    x = x - x.mean(axis=0, keepdims=True)
    std = x.std(axis=0, keepdims=True)
    std[std == 0] = 1.0
    x = x / std
    corr = np.corrcoef(x, rowvar=False)
    corr = np.nan_to_num(corr, nan=0.0, posinf=0.0, neginf=0.0)
    np.fill_diagonal(corr, 0.0)

    edges: List[Tuple[int, int, float, int]] = []
    n_genes = corr.shape[0]
    for gene_idx in range(n_genes):
        row = corr[gene_idx]
        candidate_idx = np.where(np.abs(row) >= corr_threshold)[0]
        if candidate_idx.size == 0:
            candidate_idx = np.argsort(np.abs(row))[-top_k:]
        candidate_idx = candidate_idx[candidate_idx != gene_idx]
        if candidate_idx.size > top_k:
            keep = np.argsort(np.abs(row[candidate_idx]))[-top_k:]
            candidate_idx = candidate_idx[keep]
        for nbr in candidate_idx:
            weight = float(abs(row[nbr]))
            if weight > 0:
                edges.append((gene_idx, int(nbr), weight, EDGE_TYPES["coexpression"]))
    return edges


def infer_marker_modules(
    adata: sc.AnnData,
    label_col: str,
    top_markers_per_group: int = 30,
    min_groups_per_gene: int = 2,
) -> Dict[str, List[int]]:
    tmp = adata.copy()
    sc.tl.rank_genes_groups(tmp, groupby=label_col, method="wilcoxon")
    names = tmp.uns["rank_genes_groups"]["names"]
    modules: Dict[str, List[int]] = {}
    gene_to_idx = {g: i for i, g in enumerate(tmp.var_names.tolist())}
    for group in names.dtype.names:
        selected = []
        for gene in names[group][:top_markers_per_group]:
            if gene in gene_to_idx:
                selected.append(gene_to_idx[gene])
        if selected:
            modules[str(group)] = sorted(set(selected))
    gene_membership = {}
    for module_name, gene_ids in modules.items():
        for gid in gene_ids:
            gene_membership.setdefault(gid, 0)
            gene_membership[gid] += 1
    for module_name, gene_ids in list(modules.items()):
        modules[module_name] = [gid for gid in gene_ids if gene_membership.get(gid, 0) >= min_groups_per_gene or len(gene_ids) <= 5]
    return modules


def build_marker_edges(modules: Dict[str, List[int]]) -> List[Tuple[int, int, float, int]]:
    edges: List[Tuple[int, int, float, int]] = []
    for _, gene_ids in modules.items():
        if len(gene_ids) < 2:
            continue
        weight = 1.0 / np.log2(len(gene_ids) + 1)
        for i, src in enumerate(gene_ids):
            for dst in gene_ids[i + 1:]:
                edges.append((int(src), int(dst), float(weight), EDGE_TYPES["marker"]))
                edges.append((int(dst), int(src), float(weight), EDGE_TYPES["marker"]))
    return edges


def build_orthology_like_edges(
    gene_names: Sequence[str],
    gene_features: np.ndarray,
    top_k: int = 5,
    similarity_threshold: float = 0.75,
) -> List[Tuple[int, int, float, int]]:
    prefix_groups: Dict[str, List[int]] = {}
    for idx, gene in enumerate(gene_names):
        prefix = gene[:6]
        prefix_groups.setdefault(prefix, []).append(idx)

    similarity = cosine_similarity(gene_features)
    np.fill_diagonal(similarity, 0.0)
    edges: List[Tuple[int, int, float, int]] = []
    for members in prefix_groups.values():
        if len(members) < 2:
            continue
        for src in members:
            sims = similarity[src, members]
            order = np.argsort(sims)[::-1]
            kept = 0
            for rel_idx in order:
                dst = members[rel_idx]
                if src == dst:
                    continue
                weight = float(sims[rel_idx])
                if weight < similarity_threshold:
                    continue
                edges.append((int(src), int(dst), weight, EDGE_TYPES["orthology"]))
                kept += 1
                if kept >= top_k:
                    break
    return edges


def merge_edges(edges: Sequence[Tuple[int, int, float, int]]) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    merged: Dict[Tuple[int, int], Tuple[float, int]] = {}
    for src, dst, weight, edge_type in edges:
        key = (int(src), int(dst))
        if key not in merged or weight > merged[key][0]:
            merged[key] = (float(weight), int(edge_type))
    ordered = sorted(merged.items())
    edge_index = np.array([[src, dst] for (src, dst), _ in ordered], dtype=np.int64).T
    edge_weight = np.array([value[0] for _, value in ordered], dtype=np.float32)
    edge_type = np.array([value[1] for _, value in ordered], dtype=np.int64)
    return edge_index, edge_weight, edge_type


def build_gene_graph_bundle(
    data_path: str,
    label_col: str = "Celltype",
    n_top_genes: int = 2000,
    corr_threshold: float = 0.35,
    coexpr_top_k: int = 20,
    feature_dim: int = 64,
    seed: int = 42,
) -> GeneGraphBundle:
    adata = load_adata(data_path)
    label_col = _infer_label_col(adata, label_col)
    adata = select_hvgs(adata, n_top_genes=n_top_genes)
    gene_names = adata.var_names.tolist()
    gene_features = build_gene_features(adata, feature_dim=feature_dim, seed=seed)
    coexpr_edges = build_coexpression_edges(adata, corr_threshold=corr_threshold, top_k=coexpr_top_k)
    marker_modules = infer_marker_modules(adata, label_col=label_col)
    marker_edges = build_marker_edges(marker_modules)
    orthology_edges = build_orthology_like_edges(gene_names, gene_features)
    edge_index, edge_weight, edge_type = merge_edges(coexpr_edges + marker_edges + orthology_edges)
    metadata = {
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "edge_type_counts": {
            "coexpression": len(coexpr_edges),
            "marker": len(marker_edges),
            "orthology": len(orthology_edges),
        },
        "label_col": label_col,
        "n_top_genes": n_top_genes,
    }
    return GeneGraphBundle(
        gene_names=gene_names,
        edge_index=edge_index,
        edge_weight=edge_weight,
        edge_type=edge_type,
        gene_features=gene_features,
        marker_modules=marker_modules,
        metadata=metadata,
    )


def save_gene_graph_bundle(bundle: GeneGraphBundle, output_dir: str) -> None:
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        output_path / "gene_graph.npz",
        edge_index=bundle.edge_index,
        edge_weight=bundle.edge_weight,
        edge_type=bundle.edge_type,
        gene_features=bundle.gene_features,
        gene_names=np.array(bundle.gene_names, dtype=object),
    )
    with open(output_path / "marker_modules.json", "w") as handle:
        json.dump({k: list(map(int, v)) for k, v in bundle.marker_modules.items()}, handle, indent=2)
    with open(output_path / "metadata.json", "w") as handle:
        json.dump(bundle.metadata, handle, indent=2)


def load_gene_graph_bundle(output_dir: str) -> GeneGraphBundle:
    output_path = Path(output_dir)
    graph = np.load(output_path / "gene_graph.npz", allow_pickle=True)
    with open(output_path / "marker_modules.json") as handle:
        marker_modules = json.load(handle)
    with open(output_path / "metadata.json") as handle:
        metadata = json.load(handle)
    return GeneGraphBundle(
        gene_names=list(graph["gene_names"].tolist()),
        edge_index=graph["edge_index"],
        edge_weight=graph["edge_weight"],
        edge_type=graph["edge_type"],
        gene_features=graph["gene_features"],
        marker_modules={k: list(map(int, v)) for k, v in marker_modules.items()},
        metadata=metadata,
    )
