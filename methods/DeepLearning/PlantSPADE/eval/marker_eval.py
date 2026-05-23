import numpy as np
import scanpy as sc


def marker_gene_enrichment(adata, pred_key="scspade_cluster", label_key=None, top_n=20):
    if label_key is None or label_key not in adata.obs.columns or pred_key not in adata.obs.columns:
        return {"marker_overlap": float("nan")}

    sc.tl.rank_genes_groups(adata, groupby=label_key, method="wilcoxon")
    true_markers = set()
    for group in adata.uns["rank_genes_groups"]["names"].dtype.names:
        true_markers.update(list(adata.uns["rank_genes_groups"]["names"][group][:top_n]))

    sc.tl.rank_genes_groups(adata, groupby=pred_key, method="wilcoxon")
    pred_markers = set()
    for group in adata.uns["rank_genes_groups"]["names"].dtype.names:
        pred_markers.update(list(adata.uns["rank_genes_groups"]["names"][group][:top_n]))

    overlap = len(true_markers & pred_markers) / max(len(true_markers), 1)
    return {"marker_overlap": float(overlap)}
