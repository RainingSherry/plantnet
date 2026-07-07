from __future__ import annotations

import torch
from torch import nn


class HeterogeneousDomainAdapter(nn.Module):
    """DAP-MAE style pretrain adapter with domain-specific projections."""

    def __init__(self, dim: int, n_domains: int = 3, dropout: float = 0.0):
        super().__init__()
        self.n_domains = int(n_domains)
        self.branches = nn.ModuleList(
            [
                nn.Sequential(
                    nn.Linear(dim, dim),
                    nn.LayerNorm(dim),
                    nn.GELU(),
                    nn.Dropout(dropout),
                    nn.Linear(dim, dim),
                )
                for _ in range(self.n_domains)
            ]
        )
        self.fusion_gate = nn.Sequential(nn.Linear(dim, self.n_domains), nn.Softmax(dim=-1))

    def forward(self, z: torch.Tensor, domain_id: torch.Tensor) -> torch.Tensor:
        branch_out = torch.stack([branch(z) for branch in self.branches], dim=1)
        if domain_id is None:
            weights = self.fusion_gate(z).unsqueeze(-1)
            return (weights * branch_out).sum(dim=1)
        out = torch.zeros_like(z)
        for idx in range(self.n_domains):
            mask = domain_id.long() == idx
            if mask.any():
                out[mask] = branch_out[mask, idx]
        return out


class DAPMAEScMAE(nn.Module):
    """scMAE with DAP-MAE-inspired heterogeneous domain adapter and DFG heads."""

    def __init__(self, input_dim: int, hidden_size: int = 128, decoder_hidden: int = 128, n_domains: int = 3, dropout: float = 0.0):
        super().__init__()
        self.input_dim = int(input_dim)
        self.hidden_size = int(hidden_size)
        self.n_domains = int(n_domains)
        self.encoder = nn.Sequential(
            nn.Dropout(dropout),
            nn.Linear(input_dim, 256),
            nn.LayerNorm(256),
            nn.GELU(),
            nn.Linear(256, hidden_size),
            nn.LayerNorm(hidden_size),
            nn.GELU(),
        )
        self.domain_token = nn.Embedding(n_domains, hidden_size)
        self.adapter = nn.Sequential(
            HeterogeneousDomainAdapter(hidden_size, n_domains, dropout),
            HeterogeneousDomainAdapter(hidden_size, n_domains, dropout),
        )
        self.adapter_norm = nn.LayerNorm(hidden_size)
        self.mask_predictor = nn.Linear(hidden_size, input_dim)
        self.decoder = nn.Sequential(nn.Linear(hidden_size + input_dim, decoder_hidden), nn.GELU(), nn.Linear(decoder_hidden, input_dim))
        self.domain_feature_head = nn.Sequential(nn.Linear(hidden_size, hidden_size), nn.LayerNorm(hidden_size), nn.GELU(), nn.Linear(hidden_size, hidden_size))
        self.domain_classifier = nn.Linear(hidden_size, n_domains)

    def encode(self, x: torch.Tensor, domain_id: torch.Tensor | None) -> torch.Tensor:
        z = self.encoder(x)
        if domain_id is not None:
            z = z + self.domain_token(domain_id.long().to(x.device))
        for adapter in self.adapter:
            z = adapter((z, domain_id)) if isinstance(adapter, nn.Identity) else adapter(z, domain_id)
        return self.adapter_norm(z)

    def forward(self, x: torch.Tensor, domain_id: torch.Tensor | None = None) -> dict[str, torch.Tensor]:
        z = self.encode(x, domain_id)
        mask_logits = self.mask_predictor(z)
        recon = self.decoder(torch.cat([z, mask_logits], dim=1))
        return {
            "embedding": z,
            "reconstruction": recon,
            "mask_logits": mask_logits,
            "domain_feature": self.domain_feature_head(z),
            "domain_logits": self.domain_classifier(z),
        }
