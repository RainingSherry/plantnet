from __future__ import annotations

from typing import Dict, Tuple

import torch
import torch.nn as nn
import torch.nn.functional as F

from experimental_retired_models.NeighborMix_scMAE.model import AutoEncoder as _BaseAutoEncoder


class CutAwareAutoEncoder(_BaseAutoEncoder):
    """scMAE encoder with a lightweight clustering head for cut/OT losses."""

    def __init__(self, num_genes: int, n_clusters: int, edge_feature_dim: int = 5, **kwargs):
        dropout_rate = float(kwargs.get("dropout", 0.0))
        super().__init__(num_genes=num_genes, **kwargs)
        if n_clusters <= 1:
            raise ValueError(f"n_clusters must be > 1 for cut-aware training, got {n_clusters}.")
        self.n_clusters = int(n_clusters)
        self.cluster_head = nn.Linear(self.hidden_size, self.n_clusters)
        self.edge_feature_dim = int(edge_feature_dim)
        self.edge_gate = nn.Sequential(
            nn.Linear(self.hidden_size * 4 + self.edge_feature_dim, self.hidden_size),
            nn.GELU(),
            nn.Dropout(dropout_rate),
            nn.Linear(self.hidden_size, 1),
        )

    def cluster_logits(self, latent: torch.Tensor) -> torch.Tensor:
        return self.cluster_head(latent)

    def cluster_probs(self, latent: torch.Tensor, temperature: float = 1.0) -> torch.Tensor:
        return torch.softmax(self.cluster_logits(latent) / max(float(temperature), 1e-6), dim=1)

    def edge_gate_logits(
        self,
        src_latent: torch.Tensor,
        dst_latent: torch.Tensor,
        edge_features: torch.Tensor,
    ) -> torch.Tensor:
        if src_latent.shape != dst_latent.shape:
            raise ValueError("src_latent and dst_latent must have the same shape.")
        if edge_features.shape[0] != src_latent.shape[0]:
            raise ValueError("edge_features must align with latent edge rows.")
        features = torch.cat(
            [
                src_latent,
                dst_latent,
                torch.abs(src_latent - dst_latent),
                src_latent * dst_latent,
                edge_features.to(dtype=src_latent.dtype, device=src_latent.device),
            ],
            dim=1,
        )
        return self.edge_gate(features).squeeze(-1)

    def edge_gate_scores(
        self,
        src_latent: torch.Tensor,
        dst_latent: torch.Tensor,
        edge_features: torch.Tensor,
        temperature: float = 1.0,
    ) -> torch.Tensor:
        logits = self.edge_gate_logits(src_latent, dst_latent, edge_features)
        return torch.sigmoid(logits / max(float(temperature), 1e-6))

    def loss_mask_weighted(
        self,
        x: torch.Tensor,
        y: torch.Tensor,
        mask: torch.Tensor,
        sample_weight: torch.Tensor | None = None,
        mask_loss_scale: float = 1.0,
    ) -> Tuple[torch.Tensor, torch.Tensor, Dict[str, torch.Tensor]]:
        self._check_expression_shape(x, "x")
        self._check_expression_shape(y, "y")
        self._check_expression_shape(mask, "mask")
        if x.shape != y.shape or x.shape != mask.shape:
            raise ValueError("x, y, and mask must have identical shapes.")

        mask = mask.to(dtype=x.dtype, device=x.device)
        y = y.to(dtype=x.dtype, device=x.device)
        latent, mask_logits, reconstruction = self.forward_mask(x)
        raw_mse = F.mse_loss(reconstruction, y, reduction="none")
        weights = mask * self.masked_data_weight + (1.0 - mask) * (1.0 - self.masked_data_weight)
        weighted_mse = weights * raw_mse
        if self.normalize_reconstruction_by_weight:
            rec_per = weighted_mse.sum(dim=1) / weights.sum(dim=1).clamp_min(1e-8)
        else:
            rec_per = weighted_mse.mean(dim=1)
        rec_per = (1.0 - self.mask_loss_weight) * rec_per
        mask_per = F.binary_cross_entropy_with_logits(mask_logits, mask, reduction="none").mean(dim=1)
        mask_per = self.mask_loss_weight * mask_per
        total_per = rec_per + float(mask_loss_scale) * mask_per

        if sample_weight is None:
            loss = total_per.mean()
        else:
            w = sample_weight.to(dtype=x.dtype, device=x.device).view(-1)
            loss = (total_per * w).sum() / w.sum().clamp_min(1e-8)
        parts = {
            "reconstruction_loss": rec_per.mean().detach(),
            "mask_loss": mask_per.mean().detach(),
            "total_loss": loss.detach(),
            "mask_positive_rate": mask.mean().detach(),
            "per_sample_loss": total_per.detach(),
        }
        return latent, loss, parts
