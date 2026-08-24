from __future__ import annotations

from typing import Dict, Optional, Tuple

import torch
import torch.nn as nn
import torch.nn.functional as F


class FlexibleScMAE(nn.Module):
    """A compact scMAE-compatible model with optional variant heads."""

    def __init__(
        self,
        num_genes: int,
        hidden_size: int,
        dropout: float,
        encoder_kind: str,
        token_bins: int,
        n_prototypes: int,
        use_gene_gate: bool,
    ):
        super().__init__()
        if num_genes <= 0:
            raise ValueError(f"num_genes must be positive, got {num_genes}")
        if hidden_size <= 0:
            raise ValueError(f"hidden_size must be positive, got {hidden_size}")
        self.num_genes = int(num_genes)
        self.hidden_size = int(hidden_size)
        self.token_bins = int(token_bins)
        self.encoder_kind = encoder_kind
        self.use_gene_gate = bool(use_gene_gate)

        self.gene_gate = nn.Parameter(torch.zeros(num_genes)) if use_gene_gate else None
        width = max(256, hidden_size * 2)
        self.input_norm = nn.LayerNorm(num_genes)
        self.context_proj = nn.Linear(num_genes, num_genes)

        if encoder_kind == "gated_sequence":
            self.sequence_gate = nn.Sequential(
                nn.Linear(num_genes, num_genes),
                nn.Sigmoid(),
            )
            self.encoder = nn.Sequential(
                nn.Dropout(dropout),
                nn.Linear(num_genes, width),
                nn.LayerNorm(width),
                nn.Mish(inplace=True),
                nn.Linear(width, hidden_size),
                nn.LayerNorm(hidden_size),
                nn.Mish(inplace=True),
                nn.Linear(hidden_size, hidden_size),
            )
        else:
            self.sequence_gate = None
            self.encoder = nn.Sequential(
                nn.Dropout(dropout),
                nn.Linear(num_genes, width),
                nn.LayerNorm(width),
                nn.Mish(inplace=True),
                nn.Linear(width, hidden_size),
                nn.LayerNorm(hidden_size),
                nn.Mish(inplace=True),
                nn.Linear(hidden_size, hidden_size),
            )

        self.mask_predictor = nn.Linear(hidden_size, num_genes)
        self.decoder = nn.Linear(hidden_size + num_genes, num_genes)
        self.token_head = nn.Linear(hidden_size, num_genes * token_bins)
        self.gene_feature_head = nn.Linear(hidden_size, 4)
        self.prototypes = nn.Parameter(torch.randn(max(2, n_prototypes), hidden_size) * 0.02)

    def _check_expr(self, x: torch.Tensor, name: str) -> None:
        if x.ndim != 2:
            raise ValueError(f"{name} must have shape [batch, genes], got {tuple(x.shape)}")
        if x.shape[1] != self.num_genes:
            raise ValueError(f"{name}.shape[1] must be {self.num_genes}, got {x.shape[1]}")

    def _prepare_input(self, x: torch.Tensor, context: Optional[torch.Tensor]) -> torch.Tensor:
        self._check_expr(x, "x")
        h = x
        if context is not None:
            self._check_expr(context, "context")
            h = h + 0.2 * torch.tanh(self.context_proj(context))
        if self.gene_gate is not None:
            h = h * (1.0 + torch.sigmoid(self.gene_gate))
        h = self.input_norm(h)
        if self.sequence_gate is not None:
            h = h * self.sequence_gate(h)
        return h

    def forward_mask(
        self,
        x: torch.Tensor,
        context: Optional[torch.Tensor] = None,
    ) -> Tuple[torch.Tensor, torch.Tensor, torch.Tensor, torch.Tensor, torch.Tensor]:
        h = self._prepare_input(x, context)
        latent = self.encoder(h)
        mask_logits = self.mask_predictor(latent)
        reconstruction = self.decoder(torch.cat([latent, mask_logits], dim=1))
        token_logits = self.token_head(latent).view(x.shape[0], self.num_genes, self.token_bins)
        gene_features = self.gene_feature_head(latent)
        return latent, mask_logits, reconstruction, token_logits, gene_features

    def feature(self, x: torch.Tensor) -> torch.Tensor:
        h = self._prepare_input(x, None)
        return self.encoder(h)

    def soft_assignments(self, latent: torch.Tensor) -> torch.Tensor:
        z = F.normalize(latent, dim=1)
        p = F.normalize(self.prototypes, dim=1)
        return torch.softmax(torch.matmul(z, p.t()) / 0.2, dim=1)

    def regularization_terms(self) -> Dict[str, torch.Tensor]:
        if self.gene_gate is None:
            return {}
        return {"gene_gate_l1": torch.sigmoid(self.gene_gate).mean()}

