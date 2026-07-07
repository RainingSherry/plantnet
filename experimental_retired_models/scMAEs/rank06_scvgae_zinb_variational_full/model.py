from __future__ import annotations

import torch
from torch import nn
import torch.nn.functional as F


class ScVGAEZINBScMAE(nn.Module):
    """
    scVGAE-inspired variational scMAE.

    The model keeps scMAE mask prediction and masked log-expression
    reconstruction, while adding a variational latent distribution and a ZINB
    count decoder. Counts are supplied by the runner from a non-scaled count
    branch; scaled expression is used only by the encoder.
    """

    def __init__(self, num_genes: int, hidden_size: int = 128, latent_dim: int = 128, dropout: float = 0.1):
        super().__init__()
        self.num_genes = int(num_genes)
        self.hidden_size = int(hidden_size)
        self.latent_dim = int(latent_dim)
        self.encoder_base = nn.Sequential(
            nn.Dropout(float(dropout)),
            nn.Linear(self.num_genes, 256),
            nn.LayerNorm(256),
            nn.Mish(),
            nn.Linear(256, self.hidden_size),
            nn.LayerNorm(self.hidden_size),
            nn.Mish(),
        )
        self.z_mu = nn.Linear(self.hidden_size, self.latent_dim)
        self.z_logvar = nn.Linear(self.hidden_size, self.latent_dim)
        self.latent_to_hidden = nn.Sequential(
            nn.Linear(self.latent_dim, self.hidden_size),
            nn.LayerNorm(self.hidden_size),
            nn.Mish(),
        )
        self.mask_predictor = nn.Linear(self.hidden_size, self.num_genes)
        self.expr_decoder = nn.Linear(self.hidden_size + self.num_genes, self.num_genes)
        self.zinb_mean = nn.Sequential(nn.Linear(self.hidden_size, self.num_genes), nn.Softplus())
        self.zinb_dropout = nn.Linear(self.hidden_size, self.num_genes)
        self.theta_raw = nn.Parameter(torch.zeros(self.num_genes))

    def encode(self, x: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
        h = self.encoder_base(x)
        return self.z_mu(h), self.z_logvar(h).clamp(min=-8.0, max=8.0)

    def reparameterize(self, mu: torch.Tensor, logvar: torch.Tensor) -> torch.Tensor:
        if self.training:
            eps = torch.randn_like(mu)
            return mu + eps * torch.exp(0.5 * logvar)
        return mu

    def forward(self, x: torch.Tensor, size_factor: torch.Tensor) -> dict[str, torch.Tensor]:
        mu_z, logvar_z = self.encode(x)
        z = self.reparameterize(mu_z, logvar_z)
        h = self.latent_to_hidden(z)
        mask_logits = self.mask_predictor(h)
        reconstruction = self.expr_decoder(torch.cat([h, mask_logits], dim=1))
        mean_scale = self.zinb_mean(h) + 1e-4
        mean = mean_scale / mean_scale.sum(dim=1, keepdim=True).clamp_min(1e-6)
        mean = mean * size_factor[:, None].clamp_min(1.0)
        theta = F.softplus(self.theta_raw)[None, :].expand_as(mean) + 1e-4
        dropout_logits = self.zinb_dropout(h)
        return {
            "latent": mu_z,
            "sampled_latent": z,
            "logvar": logvar_z,
            "hidden": h,
            "mask_logits": mask_logits,
            "reconstruction": reconstruction,
            "zinb_mean": mean,
            "zinb_theta": theta,
            "zinb_dropout_logits": dropout_logits,
        }

    @torch.no_grad()
    def feature(self, x: torch.Tensor) -> torch.Tensor:
        mu, _ = self.encode(x)
        return mu

    def corrupt(self, x: torch.Tensor, mask_prob: float) -> tuple[torch.Tensor, torch.Tensor]:
        should_swap = torch.bernoulli(float(mask_prob) * torch.ones_like(x)).bool()
        replacement = x[torch.randperm(x.shape[0], device=x.device)] if x.shape[0] > 1 else x
        corrupted = torch.where(should_swap, replacement, x)
        mask = (corrupted != x).float()
        empty = mask.sum(dim=1) == 0
        if bool(empty.any()):
            cols = torch.randint(0, x.shape[1], (int(empty.sum()),), device=x.device)
            corrupted[empty, cols] = replacement[empty, cols]
            mask[empty, cols] = (corrupted[empty, cols] != x[empty, cols]).float()
        return corrupted, mask

