from __future__ import annotations

from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
import scanpy as sc

from ..support_gene_attention import SparseAttentionWeights


def attention_top_genes_per_cluster(
    attention_weights: SparseAttentionWeights,
    cluster_labels: np.ndarray,
    support,
    amplitude,
    gene_names: np.ndarray,
    top_n: int = 20,
) -> pd.DataFrame:
    rows = []
    n_genes = support.shape[1]
    labels = np.asarray(cluster_labels, dtype=np.int64)
    support = support.tocsr()
    amplitude = amplitude.tocsr()
    global_support_rate = np.asarray(support.mean(axis=0)).ravel()
    global_mean_expression = np.asarray(amplitude.mean(axis=0)).ravel()

    for cluster_id in sorted(np.unique(labels)):
        cluster_cells = np.flatnonzero(labels == cluster_id)
        other_cells = np.flatnonzero(labels != cluster_id)
        if cluster_cells.size == 0:
            continue

        edge_mask = labels[attention_weights.cells] == cluster_id
        if np.any(edge_mask):
            mean_attention = np.bincount(
                attention_weights.genes[edge_mask],
                weights=attention_weights.weights[edge_mask],
                minlength=n_genes,
            ).astype(np.float64)
            mean_attention /= float(cluster_cells.size)
        else:
            mean_attention = np.zeros(n_genes, dtype=np.float64)

        top_genes = np.argsort(-mean_attention, kind="stable")[:top_n]
        support_block = support[cluster_cells][:, top_genes]
        amplitude_block = amplitude[cluster_cells][:, top_genes]
        support_rate = np.asarray(support_block.mean(axis=0)).ravel()
        mean_expression = np.asarray(amplitude_block.mean(axis=0)).ravel()

        if other_cells.size:
            other_support = np.asarray(support[other_cells][:, top_genes].mean(axis=0)).ravel()
            other_expression = np.asarray(amplitude[other_cells][:, top_genes].mean(axis=0)).ravel()
        else:
            other_support = global_support_rate[top_genes]
            other_expression = global_mean_expression[top_genes]

        cluster_specificity = support_rate - other_support
        expression_enrichment = (mean_expression + 1e-8) / (other_expression + 1e-8)

        for local_idx, gene_idx in enumerate(top_genes):
            rows.append(
                {
                    "cluster_id": int(cluster_id),
                    "gene_name": str(gene_names[gene_idx]),
                    "mean_attention": float(mean_attention[gene_idx]),
                    "support_rate": float(support_rate[local_idx]),
                    "mean_expression": float(mean_expression[local_idx]),
                    "cluster_specificity": float(cluster_specificity[local_idx]),
                    "expression_enrichment": float(expression_enrichment[local_idx]),
                }
            )
    return pd.DataFrame(rows)


def wilcoxon_deg_top_genes(
    amplitude,
    cluster_labels: np.ndarray,
    gene_names: np.ndarray,
    top_n: int = 50,
) -> pd.DataFrame:
    adata = sc.AnnData(amplitude.copy())
    adata.var_names = np.asarray(gene_names).astype(str)
    adata.obs["cluster"] = pd.Categorical(np.asarray(cluster_labels).astype(str))
    try:
        sc.tl.rank_genes_groups(adata, groupby="cluster", method="wilcoxon", n_genes=top_n)
    except Exception as exc:
        return pd.DataFrame([{"error": f"wilcoxon_failed: {exc}"}])

    result = adata.uns["rank_genes_groups"]
    rows = []
    for cluster_id in result["names"].dtype.names:
        names = result["names"][cluster_id]
        scores = result["scores"][cluster_id]
        pvals_adj = result["pvals_adj"][cluster_id]
        if "logfoldchanges" in result:
            logfoldchanges = result["logfoldchanges"][cluster_id]
        else:
            logfoldchanges = [np.nan] * len(names)
        for rank, gene in enumerate(names[:top_n], start=1):
            rows.append(
                {
                    "cluster_id": int(cluster_id) if str(cluster_id).isdigit() else str(cluster_id),
                    "rank": int(rank),
                    "gene_name": str(gene),
                    "score": float(scores[rank - 1]),
                    "pvals_adj": float(pvals_adj[rank - 1]),
                    "logfoldchange": float(logfoldchanges[rank - 1]) if len(logfoldchanges) >= rank else float("nan"),
                }
            )
    return pd.DataFrame(rows)


def attention_deg_overlap(attention_df: pd.DataFrame, deg_df: pd.DataFrame, top_n: int = 20) -> pd.DataFrame:
    if attention_df.empty or deg_df.empty or "error" in deg_df.columns:
        return pd.DataFrame()
    rows = []
    for cluster_id in sorted(set(attention_df["cluster_id"]).intersection(set(deg_df["cluster_id"]))):
        attn_genes = set(attention_df[attention_df["cluster_id"] == cluster_id].head(top_n)["gene_name"])
        deg_genes = set(deg_df[deg_df["cluster_id"] == cluster_id].head(top_n)["gene_name"])
        overlap = sorted(attn_genes.intersection(deg_genes))
        rows.append(
            {
                "cluster_id": cluster_id,
                "attention_top_n": int(len(attn_genes)),
                "deg_top_n": int(len(deg_genes)),
                "overlap_n": int(len(overlap)),
                "overlap_genes": ";".join(overlap),
            }
        )
    return pd.DataFrame(rows)


def write_marker_outputs(
    output_dir: str,
    attention_weights: SparseAttentionWeights,
    cluster_labels: np.ndarray,
    support,
    amplitude,
    gene_names: np.ndarray,
    top_n_attention: int = 20,
    top_n_deg: int = 50,
) -> None:
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    attention_df = attention_top_genes_per_cluster(
        attention_weights,
        cluster_labels,
        support,
        amplitude,
        gene_names,
        top_n=top_n_attention,
    )
    attention_df.to_csv(output / "attention_top_genes_per_cluster.csv", index=False)
    deg_df = wilcoxon_deg_top_genes(amplitude, cluster_labels, gene_names, top_n=top_n_deg)
    deg_df.to_csv(output / "wilcoxon_deg_top_genes_per_cluster.csv", index=False)
    overlap_df = attention_deg_overlap(attention_df, deg_df, top_n=top_n_attention)
    overlap_df.to_csv(output / "attention_deg_overlap.csv", index=False)
