import math

import torch
import torch.nn as nn
import torch.nn.functional as F


def timestep_embedding(timesteps: torch.Tensor, dim: int) -> torch.Tensor:
    half = dim // 2
    freqs = torch.exp(
        -math.log(10000) * torch.arange(half, device=timesteps.device, dtype=torch.float32) / max(half - 1, 1)
    )
    args = timesteps[:, None].float() * freqs[None, :]
    emb = torch.cat([torch.sin(args), torch.cos(args)], dim=-1)
    if dim % 2:
        emb = F.pad(emb, (0, 1))
    return emb


class DiffusionBlock(nn.Module):
    def __init__(self, hidden_dim: int, time_dim: int, dropout: float):
        super().__init__()
        self.time_proj = nn.Sequential(
            nn.Linear(time_dim, hidden_dim),
            nn.Mish(inplace=True),
            nn.Linear(hidden_dim, hidden_dim),
        )
        self.net = nn.Sequential(
            nn.Linear(hidden_dim, hidden_dim * 2),
            nn.LayerNorm(hidden_dim * 2),
            nn.Mish(inplace=True),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim * 2, hidden_dim),
        )
        self.time_scale = nn.Parameter(torch.tensor(0.1))

    def forward(self, x: torch.Tensor, time_emb: torch.Tensor) -> torch.Tensor:
        h = x + self.time_scale * self.time_proj(time_emb)
        return h + self.net(h)


class LatentDenoiser(nn.Module):
    def __init__(
        self,
        latent_dim: int,
        hidden_dim: int = 256,
        num_layers: int = 3,
        time_dim: int = 128,
        dropout: float = 0.1,
    ):
        super().__init__()
        self.time_dim = time_dim
        self.time_mlp = nn.Sequential(
            nn.Linear(time_dim, hidden_dim),
            nn.Mish(inplace=True),
            nn.Linear(hidden_dim, time_dim),
        )
        self.input_proj = nn.Linear(latent_dim, hidden_dim)
        self.blocks = nn.ModuleList([DiffusionBlock(hidden_dim, time_dim, dropout) for _ in range(num_layers)])
        self.output_proj = nn.Linear(hidden_dim, latent_dim)

    def forward(self, z_t: torch.Tensor, timesteps: torch.Tensor) -> torch.Tensor:
        t_emb = timestep_embedding(timesteps, self.time_dim)
        t_emb = self.time_mlp(t_emb)
        h = self.input_proj(z_t)
        for block in self.blocks:
            h = block(h, t_emb)
        return self.output_proj(h)


class LatentDiffusionAE(nn.Module):
    """Sparse reconstruction autoencoder with a DDPM-style latent denoiser."""

    def __init__(
        self,
        num_genes: int,
        latent_dim: int = 32,
        hidden_dims=None,
        diffusion_steps: int = 100,
        dropout: float = 0.1,
    ):
        super().__init__()
        if hidden_dims is None:
            hidden_dims = [512, 256]

        self.num_genes = num_genes
        self.latent_dim = latent_dim
        self.diffusion_steps = diffusion_steps

        enc_layers = []
        in_dim = num_genes
        for hidden_dim in hidden_dims:
            enc_layers.extend(
                [
                    nn.Dropout(dropout),
                    nn.Linear(in_dim, hidden_dim),
                    nn.LayerNorm(hidden_dim),
                    nn.Mish(inplace=True),
                ]
            )
            in_dim = hidden_dim
        enc_layers.append(nn.Linear(in_dim, latent_dim))
        self.encoder = nn.Sequential(*enc_layers)

        dec_layers = []
        in_dim = latent_dim
        for hidden_dim in reversed(hidden_dims):
            dec_layers.extend(
                [
                    nn.Linear(in_dim, hidden_dim),
                    nn.LayerNorm(hidden_dim),
                    nn.Mish(inplace=True),
                    nn.Dropout(dropout),
                ]
            )
            in_dim = hidden_dim
        dec_layers.append(nn.Linear(in_dim, num_genes))
        self.decoder = nn.Sequential(*dec_layers)

        self.denoiser = LatentDenoiser(
            latent_dim=latent_dim,
            hidden_dim=max(128, latent_dim * 4),
            num_layers=3,
            time_dim=128,
            dropout=dropout,
        )

        betas = self._cosine_beta_schedule(diffusion_steps)
        alphas = 1.0 - betas
        alphas_cumprod = torch.cumprod(alphas, dim=0)
        self.register_buffer("betas", betas)
        self.register_buffer("alphas_cumprod", alphas_cumprod)
        self.register_buffer("sqrt_alphas_cumprod", torch.sqrt(alphas_cumprod))
        self.register_buffer("sqrt_one_minus_alphas_cumprod", torch.sqrt(1.0 - alphas_cumprod))

    def _cosine_beta_schedule(self, timesteps: int, s: float = 0.008) -> torch.Tensor:
        steps = timesteps + 1
        x = torch.linspace(0, timesteps, steps)
        alphas_cumprod = torch.cos(((x / timesteps) + s) / (1 + s) * torch.pi * 0.5) ** 2
        alphas_cumprod = alphas_cumprod / alphas_cumprod[0]
        betas = 1.0 - (alphas_cumprod[1:] / alphas_cumprod[:-1])
        return torch.clamp(betas, 0.0001, 0.02)

    def _extract(self, values: torch.Tensor, timesteps: torch.Tensor, shape: torch.Size) -> torch.Tensor:
        out = values.to(timesteps.device).gather(0, timesteps)
        return out.view(timesteps.shape[0], *((1,) * (len(shape) - 1)))

    def encode(self, x: torch.Tensor) -> torch.Tensor:
        return torch.clamp(self.encoder(x), min=-10.0, max=10.0)

    def decode(self, z: torch.Tensor) -> torch.Tensor:
        return self.decoder(z)

    def q_sample(self, z0: torch.Tensor, timesteps: torch.Tensor, noise: torch.Tensor = None) -> torch.Tensor:
        if noise is None:
            noise = torch.randn_like(z0)
        sqrt_alpha = self._extract(self.sqrt_alphas_cumprod, timesteps, z0.shape)
        sqrt_one_minus = self._extract(self.sqrt_one_minus_alphas_cumprod, timesteps, z0.shape)
        return sqrt_alpha * z0 + sqrt_one_minus * noise

    def masked_mse(self, pred: torch.Tensor, target: torch.Tensor, mask: torch.Tensor = None) -> torch.Tensor:
        loss = F.mse_loss(pred, target, reduction="none")
        if mask is None:
            return loss.mean()
        denom = mask.sum().clamp_min(1.0)
        return (loss * mask).sum() / denom

    def forward(
        self,
        x: torch.Tensor,
        mask: torch.Tensor = None,
        return_recon: bool = True,
        recon_from_denoised: bool = False,
    ) -> dict:
        z = self.encode(x)
        batch_size = x.shape[0]
        timesteps = torch.randint(0, self.diffusion_steps, (batch_size,), device=x.device)
        z_noisy = self.q_sample(z, timesteps)
        z_denoised = self.denoiser(z_noisy, timesteps)
        diffusion_loss = F.mse_loss(z_denoised, z)

        x_recon = None
        recon_loss = torch.tensor(0.0, device=x.device)
        if return_recon:
            z_for_recon = z_denoised if recon_from_denoised else z
            x_recon = self.decode(z_for_recon)
            recon_loss = self.masked_mse(x_recon, x, mask)

        return {
            "z": z,
            "z_denoised": z_denoised,
            "x_recon": x_recon,
            "losses": {
                "diffusion": diffusion_loss,
                "recon": recon_loss,
            },
        }

    @torch.no_grad()
    def denoise_embedding(self, z: torch.Tensor, start_frac: float = 0.35) -> torch.Tensor:
        """Deterministically refine an encoded cell latent instead of sampling from pure noise."""
        self.eval()
        start_t = int(round((self.diffusion_steps - 1) * start_frac))
        start_t = max(0, min(start_t, self.diffusion_steps - 1))
        t = torch.full((z.shape[0],), start_t, device=z.device, dtype=torch.long)
        z_t = self.q_sample(z, t, noise=torch.zeros_like(z))

        for step in reversed(range(start_t + 1)):
            step_t = torch.full((z.shape[0],), step, device=z.device, dtype=torch.long)
            pred_x0 = self.denoiser(z_t, step_t)
            if step == 0:
                z_t = pred_x0
                break
            prev_t = torch.full((z.shape[0],), step - 1, device=z.device, dtype=torch.long)
            alpha_t = self._extract(self.alphas_cumprod, step_t, z.shape)
            alpha_prev = self._extract(self.alphas_cumprod, prev_t, z.shape)
            eps = (z_t - torch.sqrt(alpha_t) * pred_x0) / torch.sqrt(1.0 - alpha_t + 1e-8)
            z_t = torch.sqrt(alpha_prev) * pred_x0 + torch.sqrt(1.0 - alpha_prev) * eps
            z_t = torch.clamp(z_t, -10.0, 10.0)
        return torch.clamp(z_t, -10.0, 10.0)

