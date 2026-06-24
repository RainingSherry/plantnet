from __future__ import annotations

import numpy as np
import torch


def summarize_mask(mask: torch.Tensor, eligibility: torch.Tensor | None = None) -> dict:
    m = mask.detach().float().cpu()
    per_gene = m.mean(dim=0).numpy()
    sorted_x = np.sort(per_gene)
    denom = sorted_x.sum() + 1.0e-8
    n = len(sorted_x)
    gini = float(np.sum((2 * np.arange(1, n + 1) - n - 1) * sorted_x) / (n * denom)) if n else 0.0
    p = per_gene / (per_gene.sum() + 1.0e-8)
    entropy = float(-(p * np.log(p + 1.0e-8)).sum() / np.log(max(2, len(p))))
    out = {
        "selected_mask_ratio": float(m.mean().item()),
        "per_gene_mask_rate_max": float(per_gene.max()) if per_gene.size else 0.0,
        "per_gene_mask_rate_min": float(per_gene.min()) if per_gene.size else 0.0,
        "normalized_mask_entropy": entropy,
        "mask_gini": max(0.0, gini),
        "top_masked_genes": np.argsort(-per_gene)[:10].astype(int).tolist(),
    }
    if eligibility is not None:
        out["eligible_ratio"] = float(eligibility.detach().float().mean().cpu())
    return out

