from __future__ import annotations

import torch
from torch import nn
import torch.nn.functional as F


class GeneStructuralMaskPolicy(nn.Module):
    """Gene-level structural masking policy adapted from SMMM's saliency idea."""

    def __init__(self, n_genes: int, mask_prob: float = 0.4, temperature: float = 1.0):
        super().__init__()
        self.n_genes = int(n_genes)
        self.mask_prob = float(mask_prob)
        self.temperature = float(temperature)
        self.policy = nn.Sequential(
            nn.Linear(4, 32),
            nn.GELU(),
            nn.Linear(32, 1),
        )
        self.global_bias = nn.Parameter(torch.zeros(n_genes))

    @staticmethod
    def _normalize(v: torch.Tensor) -> torch.Tensor:
        return (v - v.mean()) / (v.std(unbiased=False) + 1e-6)

    def forward(self, gene_mean: torch.Tensor, gene_var: torch.Tensor, gene_dropout: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
        # High variance, high dropout, low mean genes are more likely to be markers or rare-cell signals.
        marker_risk = self._normalize(gene_var) + self._normalize(gene_dropout) - 0.5 * self._normalize(gene_mean)
        marker_risk = torch.sigmoid(marker_risk)
        stats = torch.stack(
            [
                self._normalize(gene_mean),
                self._normalize(gene_var),
                self._normalize(gene_dropout),
                self._normalize(marker_risk),
            ],
            dim=1,
        )
        learned = self.policy(stats).squeeze(1) + self.global_bias
        raw = learned / max(self.temperature, 1e-4)
        # Mask unprotected genes more often and protect marker-risk genes from destructive masking.
        protect = marker_risk.detach()
        prob = torch.sigmoid(raw) * (1.0 - 0.65 * protect)
        prob = prob * (self.mask_prob / (prob.mean().detach() + 1e-6))
        prob = prob.clamp(0.02, 0.90)
        return prob, marker_risk


class StructuralMaskScMAE(nn.Module):
    """Independent scMAE-style model with a structural gene mask policy."""

    def __init__(
        self,
        input_dim: int,
        hidden_dim: int = 512,
        latent_dim: int = 32,
        dropout: float = 0.1,
        mask_prob: float = 0.4,
        policy_weight: float = 1.0,
    ):
        super().__init__()
        self.input_dim = int(input_dim)
        self.policy_weight = float(policy_weight)
        self.mask_policy = GeneStructuralMaskPolicy(input_dim, mask_prob=mask_prob)
        self.encoder = nn.Sequential(
            nn.Linear(input_dim, hidden_dim),
            nn.LayerNorm(hidden_dim),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim, hidden_dim // 2),
            nn.LayerNorm(hidden_dim // 2),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim // 2, latent_dim),
        )
        self.mask_predictor = nn.Sequential(
            nn.Linear(latent_dim, hidden_dim // 2),
            nn.GELU(),
            nn.Linear(hidden_dim // 2, input_dim),
        )
        self.decoder = nn.Sequential(
            nn.Linear(latent_dim, hidden_dim // 2),
            nn.LayerNorm(hidden_dim // 2),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim // 2, hidden_dim),
            nn.GELU(),
            nn.Linear(hidden_dim, input_dim),
        )

    def corrupt(self, x: torch.Tensor, gene_mean: torch.Tensor, gene_var: torch.Tensor, gene_dropout: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor, torch.Tensor]:
        mask_prob, marker_risk = self.mask_policy(gene_mean, gene_var, gene_dropout)
        if self.policy_weight <= 0:
            mask_prob = torch.full_like(mask_prob, self.mask_policy.mask_prob)
        random = torch.rand_like(x)
        mask = (random < mask_prob.unsqueeze(0)).float()
        corrupted = x * (1.0 - mask)
        return corrupted, mask, mask_prob, marker_risk

    def forward(self, x: torch.Tensor, gene_mean: torch.Tensor, gene_var: torch.Tensor, gene_dropout: torch.Tensor) -> dict:
        corrupted, mask, mask_prob, marker_risk = self.corrupt(x, gene_mean, gene_var, gene_dropout)
        latent = self.encoder(corrupted)
        return {
            "latent": latent,
            "reconstruction": self.decoder(latent),
            "mask_logits": self.mask_predictor(latent),
            "mask": mask,
            "mask_prob": mask_prob,
            "marker_risk": marker_risk,
            "corrupted": corrupted,
        }

    def feature(self, x: torch.Tensor) -> torch.Tensor:
        return self.encoder(x)

    @staticmethod
    def entropy_from_prob(prob: torch.Tensor) -> torch.Tensor:
        p = prob.clamp(1e-5, 1.0 - 1e-5)
        return -(p * torch.log(p) + (1.0 - p) * torch.log(1.0 - p))
