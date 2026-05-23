import math

import torch
import torch.nn as nn


class LatentDiffusionPrior(nn.Module):
    def __init__(self, latent_dim: int, hidden_dim: int = 128, timesteps: int = 100):
        super().__init__()
        self.timesteps = timesteps
        self.denoiser = nn.Sequential(
            nn.Linear(latent_dim + 1, hidden_dim),
            nn.LayerNorm(hidden_dim),
            nn.GELU(),
            nn.Linear(hidden_dim, hidden_dim),
            nn.LayerNorm(hidden_dim),
            nn.GELU(),
            nn.Linear(hidden_dim, latent_dim),
        )

    def _time_embedding(self, t, device):
        if not torch.is_tensor(t):
            t = torch.tensor(t, device=device, dtype=torch.float32)
        if t.dim() == 0:
            t = t[None]
        return t.float().unsqueeze(-1)

    def q_sample(self, z0, t, noise=None):
        if noise is None:
            noise = torch.randn_like(z0)
        t = t.float().view(-1, 1)
        alpha = torch.cos((t + 0.008) / 1.008 * math.pi / 2) ** 2
        alpha = alpha.clamp(min=1e-4, max=1.0)
        return alpha.sqrt() * z0 + (1 - alpha).sqrt() * noise, noise

    def predict_noise(self, zt, t):
        t_emb = self._time_embedding(t, zt.device)
        return self.denoiser(torch.cat([zt, t_emb], dim=-1))

    def loss(self, z0):
        batch = z0.shape[0]
        t = torch.rand(batch, device=z0.device)
        zt, noise = self.q_sample(z0, t)
        pred = self.predict_noise(zt, t)
        return torch.mean((pred - noise) ** 2)

    @torch.no_grad()
    def denoise(self, z):
        for step in reversed(range(self.timesteps)):
            t = torch.full((z.shape[0],), float(step) / float(max(self.timesteps - 1, 1)), device=z.device)
            noise = self.predict_noise(z, t)
            z = z - noise / max(self.timesteps, 1)
        return z
