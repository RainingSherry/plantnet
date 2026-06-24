from __future__ import annotations

import torch
import torch.nn as nn

from .relaxed_topk import relaxed_topk_straight_through


class AdversarialMaskGenerator(nn.Module):
    """Generator that only scores mask positions. It never emits replacement values."""

    def __init__(self, n_genes: int, hidden_dim: int, mask_ratio: float) -> None:
        super().__init__()
        self.n_genes = int(n_genes)
        self.mask_ratio = float(mask_ratio)
        self.cell_mlp = nn.Sequential(
            nn.Linear(n_genes, hidden_dim),
            nn.LayerNorm(hidden_dim),
            nn.Mish(),
            nn.Linear(hidden_dim, hidden_dim),
        )
        self.gene_embedding = nn.Parameter(torch.empty(n_genes, hidden_dim))
        nn.init.normal_(self.gene_embedding, std=0.02)
        self.value_projection = nn.Linear(1, hidden_dim)
        self.score_mlp = nn.Sequential(nn.Linear(hidden_dim, hidden_dim), nn.Mish(), nn.Linear(hidden_dim, 1))

    def forward(self, x: torch.Tensor, eligibility: torch.Tensor, tau: float, add_gumbel: bool = True):
        b, g = x.shape
        if g != self.n_genes:
            raise ValueError(f"Expected {self.n_genes} genes, got {g}")
        cell_context = self.cell_mlp(x).unsqueeze(1)
        gene = self.gene_embedding.unsqueeze(0)
        value = self.value_projection(x.unsqueeze(-1))
        logits = self.score_mlp(torch.tanh(cell_context + gene + value)).squeeze(-1)
        logits = logits.masked_fill(~eligibility.bool(), -1.0e9)
        base_k = int(round(self.mask_ratio * g))
        eligible_count = eligibility.long().sum(dim=1)
        k_i = torch.minimum(eligible_count, torch.full((b,), base_k, dtype=torch.long, device=x.device))
        hard, soft, st = relaxed_topk_straight_through(logits, k_i, tau, eligibility, add_gumbel=add_gumbel)
        deficit = (base_k - k_i).clamp_min(0)
        return logits, hard, soft, st, {
            "mask_type": "adversarial",
            "budget_per_cell": base_k,
            "k_i": k_i.detach(),
            "budget_deficit": deficit.detach(),
            "budget_deficit_rate": float((deficit > 0).float().mean().detach().cpu()),
        }
