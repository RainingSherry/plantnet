from __future__ import annotations

import torch
from torch import nn
import torch.nn.functional as F


def sinusoidal_time_embedding(t: torch.Tensor, dim: int) -> torch.Tensor:
    half = dim // 2
    freqs = torch.exp(torch.linspace(0, -torch.log(torch.tensor(10000.0, device=t.device)), half, device=t.device))
    args = t.float().view(-1, 1) * freqs.view(1, -1)
    emb = torch.cat([torch.sin(args), torch.cos(args)], dim=1)
    if emb.shape[1] < dim:
        emb = F.pad(emb, (0, dim - emb.shape[1]))
    return emb


class ScVAEDerScMAE(nn.Module):
    """scMAE backbone with VAE latent adapter and latent diffusion denoiser."""

    def __init__(self, input_dim: int, hidden_size: int = 128, decoder_hidden: int = 128, dropout: float = 0.0):
        super().__init__()
        self.input_dim = int(input_dim)
        self.hidden_size = int(hidden_size)
        self.encoder = nn.Sequential(
            nn.Dropout(p=dropout),
            nn.Linear(input_dim, 256),
            nn.LayerNorm(256),
            nn.Mish(inplace=True),
            nn.Linear(256, hidden_size),
            nn.LayerNorm(hidden_size),
            nn.Mish(inplace=True),
            nn.Linear(hidden_size, hidden_size),
        )
        self.mu_head = nn.Linear(hidden_size, hidden_size)
        self.logvar_head = nn.Linear(hidden_size, hidden_size)
        self.latent_norm = nn.LayerNorm(hidden_size)
        self.mask_predictor = nn.Linear(hidden_size, input_dim)
        self.decoder = nn.Sequential(
            nn.Linear(hidden_size + input_dim, decoder_hidden),
            nn.Mish(inplace=True),
            nn.Linear(decoder_hidden, input_dim),
        )
        self.nb_mean_decoder = nn.Sequential(nn.Linear(hidden_size, decoder_hidden), nn.Mish(inplace=True), nn.Linear(decoder_hidden, input_dim))
        self.nb_dropout_decoder = nn.Sequential(nn.Linear(hidden_size, decoder_hidden), nn.Mish(inplace=True), nn.Linear(decoder_hidden, input_dim))
        self.nb_log_theta = nn.Parameter(torch.zeros(input_dim))
        self.time_projector = nn.Sequential(nn.Linear(hidden_size, hidden_size), nn.Mish(inplace=True), nn.Linear(hidden_size, hidden_size))
        self.latent_denoiser = nn.Sequential(
            nn.Linear(hidden_size * 2, decoder_hidden),
            nn.LayerNorm(decoder_hidden),
            nn.Mish(inplace=True),
            nn.Linear(decoder_hidden, hidden_size),
        )

    def encode(self, x: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
        base = self.encoder(x)
        mu = self.latent_norm(base + 0.1 * self.mu_head(base))
        logvar = self.logvar_head(base).clamp(-6.0, 4.0)
        return base, mu, logvar

    def reparameterize(self, mu: torch.Tensor, logvar: torch.Tensor) -> torch.Tensor:
        if not self.training:
            return mu
        return mu + torch.randn_like(mu) * torch.exp(0.5 * logvar)

    def forward(self, x: torch.Tensor) -> dict[str, torch.Tensor]:
        _, mu, logvar = self.encode(x)
        z_sample = self.reparameterize(mu, logvar)
        z_for_decode = 0.8 * mu + 0.2 * z_sample
        mask_logits = self.mask_predictor(z_for_decode)
        recon = self.decoder(torch.cat([z_for_decode, mask_logits], dim=1))
        return {
            "embedding": mu,
            "z_sample": z_sample,
            "reconstruction": recon,
            "mask_logits": mask_logits,
            "mu": mu,
            "logvar": logvar,
            "nb_mean": F.softplus(self.nb_mean_decoder(mu)) + 1e-4,
            "nb_dropout_logits": self.nb_dropout_decoder(mu),
            "nb_theta": F.softplus(self.nb_log_theta) + 1e-4,
        }

    def predict_noise(self, z_noisy: torch.Tensor, t: torch.Tensor) -> torch.Tensor:
        t_emb = self.time_projector(sinusoidal_time_embedding(t, self.hidden_size))
        return self.latent_denoiser(torch.cat([z_noisy, t_emb], dim=1))
