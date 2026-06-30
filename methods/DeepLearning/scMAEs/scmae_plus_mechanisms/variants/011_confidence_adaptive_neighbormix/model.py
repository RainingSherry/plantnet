from __future__ import annotations

from common.backbone import build_mechanism_scmae


def build_model(num_genes: int, args, n_clusters: int):
    return build_mechanism_scmae(num_genes=num_genes, args=args, n_clusters=n_clusters)
