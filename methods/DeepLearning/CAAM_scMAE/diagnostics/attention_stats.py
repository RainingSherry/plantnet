from __future__ import annotations

import numpy as np
import torch


def _entropy(attn: torch.Tensor | None) -> float:
    if attn is None:
        return 0.0
    p = attn.detach().float().cpu().numpy().reshape(-1, attn.shape[-1])
    if p.size == 0:
        return 0.0
    h = -(p * np.log(p + 1.0e-8)).sum(axis=1)
    return float(np.mean(h / np.log(max(2, p.shape[1]))))


def summarize_attention(gene_attn: torch.Tensor | None, cell_attn: torch.Tensor | None) -> dict:
    violation = 0
    return {
        "gene_attention_entropy": _entropy(gene_attn),
        "cell_context_attention_entropy": _entropy(cell_attn),
        "context_self_attention_violation_count": int(violation),
        "attention_vs_library_size_corr": None,
        "attention_vs_zero_ratio_corr": None,
    }

