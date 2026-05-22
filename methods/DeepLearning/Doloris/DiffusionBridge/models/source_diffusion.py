from __future__ import annotations

from dataclasses import dataclass
import math

import torch
import torch.nn.functional as F
from torch import nn


def _extract(coeff: torch.Tensor, t: torch.Tensor, x_shape: torch.Size) -> torch.Tensor:
    values = coeff.gather(0, t)
    return values.view(-1, *([1] * (len(x_shape) - 1)))


@dataclass
class BridgeOutput:
    shared_latent: torch.Tensor
    target_embedding: torch.Tensor
    support_anchor: torch.Tensor


class SinusoidalTimeEmbedding(nn.Module):
    def __init__(self, dim: int):
        super().__init__()
        self.dim = dim

    def forward(self, t: torch.Tensor) -> torch.Tensor:
        half_dim = self.dim // 2
        scale = math.log(10000.0) / max(half_dim - 1, 1)
        freq = torch.exp(torch.arange(half_dim, device=t.device, dtype=torch.float32) * -scale)
        angles = t.float().unsqueeze(1) * freq.unsqueeze(0)
        emb = torch.cat([angles.sin(), angles.cos()], dim=1)
        if self.dim % 2 == 1:
            emb = F.pad(emb, (0, 1))
        return emb


class ResidualMLPBlock(nn.Module):
    def __init__(self, dim: int, dropout: float = 0.0):
        super().__init__()
        self.block = nn.Sequential(
            nn.Linear(dim, dim),
            nn.SiLU(),
            nn.Dropout(dropout),
            nn.Linear(dim, dim),
        )
        self.norm = nn.LayerNorm(dim)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.norm(x + self.block(x))


class TimeConditionedDenoiser(nn.Module):
    def __init__(
        self,
        latent_dim: int,
        hidden_dim: int,
        time_dim: int,
        condition_dim: int = 0,
        depth: int = 3,
        dropout: float = 0.0,
    ):
        super().__init__()
        self.time_embedding = SinusoidalTimeEmbedding(time_dim)
        self.time_proj = nn.Sequential(nn.Linear(time_dim, hidden_dim), nn.SiLU(), nn.Linear(hidden_dim, hidden_dim))
        self.input_proj = nn.Linear(latent_dim, hidden_dim)
        self.condition_proj = nn.Linear(condition_dim, hidden_dim) if condition_dim > 0 else None
        self.blocks = nn.ModuleList([ResidualMLPBlock(hidden_dim, dropout=dropout) for _ in range(depth)])
        self.output = nn.Linear(hidden_dim, latent_dim)

    def forward(self, x_t: torch.Tensor, t: torch.Tensor, condition: torch.Tensor | None = None) -> torch.Tensor:
        hidden = self.input_proj(x_t) + self.time_proj(self.time_embedding(t))
        if self.condition_proj is not None and condition is not None:
            hidden = hidden + self.condition_proj(condition)
        for block in self.blocks:
            hidden = block(hidden)
        return self.output(hidden)


class LatentDomainDiffusion(nn.Module):
    def __init__(
        self,
        domain_dim: int,
        shared_dim: int,
        hidden_dim: int = 256,
        time_dim: int = 128,
        condition_dim: int = 0,
        num_steps: int = 50,
        beta_start: float = 1e-4,
        beta_end: float = 2e-2,
        dropout: float = 0.0,
    ):
        super().__init__()
        self.domain_dim = domain_dim
        self.shared_dim = shared_dim
        self.condition_dim = condition_dim
        self.num_steps = num_steps

        self.domain_encoder = nn.Sequential(
            nn.Linear(domain_dim, hidden_dim),
            nn.SiLU(),
            nn.Linear(hidden_dim, hidden_dim),
            nn.SiLU(),
            nn.Linear(hidden_dim, shared_dim),
        )
        self.domain_decoder = nn.Sequential(
            nn.Linear(shared_dim, hidden_dim),
            nn.SiLU(),
            nn.Linear(hidden_dim, hidden_dim),
            nn.SiLU(),
            nn.Linear(hidden_dim, domain_dim),
        )
        self.denoiser = TimeConditionedDenoiser(
            latent_dim=shared_dim,
            hidden_dim=hidden_dim,
            time_dim=time_dim,
            condition_dim=condition_dim,
            depth=3,
            dropout=dropout,
        )

        betas = torch.linspace(beta_start, beta_end, num_steps, dtype=torch.float32)
        alphas = 1.0 - betas
        alphas_cumprod = torch.cumprod(alphas, dim=0)
        alphas_cumprod_prev = torch.cat([torch.ones(1, dtype=torch.float32), alphas_cumprod[:-1]], dim=0)
        alphas_cumprod_next = torch.cat([alphas_cumprod[1:], torch.zeros(1, dtype=torch.float32)], dim=0)

        self.register_buffer('betas', betas)
        self.register_buffer('alphas_cumprod', alphas_cumprod)
        self.register_buffer('alphas_cumprod_prev', alphas_cumprod_prev)
        self.register_buffer('alphas_cumprod_next', alphas_cumprod_next)
        self.register_buffer('sqrt_alphas_cumprod', torch.sqrt(alphas_cumprod))
        self.register_buffer('sqrt_one_minus_alphas_cumprod', torch.sqrt(1.0 - alphas_cumprod))
        self.register_buffer('sqrt_recip_alphas_cumprod', torch.sqrt(1.0 / alphas_cumprod.clamp_min(1e-8)))
        self.register_buffer(
            'sqrt_recipm1_alphas_cumprod',
            torch.sqrt((1.0 / alphas_cumprod.clamp_min(1e-8)) - 1.0),
        )

    def encode_domain(self, x: torch.Tensor) -> torch.Tensor:
        return self.domain_encoder(x)

    def decode_domain(self, z: torch.Tensor) -> torch.Tensor:
        return self.domain_decoder(z)

    def q_sample(self, x_start: torch.Tensor, t: torch.Tensor, noise: torch.Tensor | None = None) -> torch.Tensor:
        if noise is None:
            noise = torch.randn_like(x_start)
        return _extract(self.sqrt_alphas_cumprod, t, x_start.shape) * x_start + _extract(
            self.sqrt_one_minus_alphas_cumprod, t, x_start.shape
        ) * noise

    def predict_x0_from_eps(self, x_t: torch.Tensor, t: torch.Tensor, eps: torch.Tensor) -> torch.Tensor:
        return _extract(self.sqrt_recip_alphas_cumprod, t, x_t.shape) * x_t - _extract(
            self.sqrt_recipm1_alphas_cumprod, t, x_t.shape
        ) * eps

    def predict_eps(self, x_t: torch.Tensor, t: torch.Tensor, condition: torch.Tensor | None = None) -> torch.Tensor:
        return self.denoiser(x_t, t, condition=condition)

    def denoise_shared(self, x_t: torch.Tensor, t: torch.Tensor, condition: torch.Tensor | None = None) -> torch.Tensor:
        eps = self.predict_eps(x_t, t, condition=condition)
        return self.predict_x0_from_eps(x_t, t, eps)

    def training_loss(
        self,
        domain_x: torch.Tensor,
        condition: torch.Tensor | None = None,
        recon_weight: float = 1.0,
        prior_weight: float = 1e-3,
        zero_weight: float = 0.25,
    ) -> dict[str, torch.Tensor]:
        shared_seed = self.encode_domain(domain_x)
        t = torch.randint(0, self.num_steps, (domain_x.shape[0],), device=domain_x.device, dtype=torch.long)
        noise = torch.randn_like(shared_seed)
        x_t = self.q_sample(shared_seed, t, noise=noise)
        pred_eps = self.predict_eps(x_t, t, condition=condition)
        pred_shared0 = self.predict_x0_from_eps(x_t, t, pred_eps)
        pred_domain = self.decode_domain(pred_shared0)

        diffusion_loss = F.mse_loss(pred_eps, noise)

        if domain_x.shape[-1] == pred_domain.shape[-1]:
            weights = torch.where(domain_x > 0, torch.ones_like(domain_x), torch.full_like(domain_x, zero_weight))
            recon_loss = ((pred_domain - domain_x).pow(2) * weights).mean()
        else:
            recon_loss = F.mse_loss(pred_domain, domain_x)

        prior_loss = shared_seed.pow(2).mean()
        loss = diffusion_loss + recon_weight * recon_loss + prior_weight * prior_loss
        return {
            'loss': loss,
            'diffusion_loss': diffusion_loss.detach(),
            'recon_loss': recon_loss.detach(),
            'prior_loss': prior_loss.detach(),
        }

    # DDIM inversion pushes a domain sample toward the shared isotropic Gaussian.
    # This mirrors DOLORIS: source/teacher domains meet in the same Gaussian endpoint.
    def ddim_reverse_step(self, x_t: torch.Tensor, t: torch.Tensor, condition: torch.Tensor | None = None) -> torch.Tensor:
        pred_x0 = self.denoise_shared(x_t, t, condition=condition)
        eps = (x_t - _extract(self.sqrt_alphas_cumprod, t, x_t.shape) * pred_x0) / _extract(
            self.sqrt_one_minus_alphas_cumprod, t, x_t.shape
        ).clamp_min(1e-8)
        alpha_bar_next = _extract(self.alphas_cumprod_next, t, x_t.shape)
        return pred_x0 * alpha_bar_next.sqrt() + (1.0 - alpha_bar_next).clamp_min(0.0).sqrt() * eps

    # DDIM sampling pulls the shared Gaussian back to a structured domain state.
    # We keep eta=0 so the bridge is deterministic and the latent correspondence is stable.
    def ddim_sample_step(self, x_t: torch.Tensor, t: torch.Tensor, condition: torch.Tensor | None = None) -> torch.Tensor:
        pred_x0 = self.denoise_shared(x_t, t, condition=condition)
        eps = (x_t - _extract(self.sqrt_alphas_cumprod, t, x_t.shape) * pred_x0) / _extract(
            self.sqrt_one_minus_alphas_cumprod, t, x_t.shape
        ).clamp_min(1e-8)
        alpha_bar_prev = _extract(self.alphas_cumprod_prev, t, x_t.shape)
        return pred_x0 * alpha_bar_prev.sqrt() + (1.0 - alpha_bar_prev).clamp_min(0.0).sqrt() * eps

    def ddim_reverse_sample_loop(self, image: torch.Tensor, condition: torch.Tensor | None = None) -> torch.Tensor:
        state = self.encode_domain(image)
        for step in range(self.num_steps):
            t = torch.full((image.shape[0],), step, device=image.device, dtype=torch.long)
            state = self.ddim_reverse_step(state, t, condition=condition)
        return state

    def ddim_sample_loop(self, noise: torch.Tensor, condition: torch.Tensor | None = None) -> torch.Tensor:
        state = noise
        for step in reversed(range(self.num_steps)):
            t = torch.full((noise.shape[0],), step, device=noise.device, dtype=torch.long)
            state = self.ddim_sample_step(state, t, condition=condition)
        return self.decode_domain(state)


class SourceDiffusion(LatentDomainDiffusion):
    def __init__(
        self,
        domain_dim: int,
        shared_dim: int = 64,
        hidden_dim: int = 256,
        time_dim: int = 128,
        num_steps: int = 50,
        dropout: float = 0.0,
    ):
        super().__init__(
            domain_dim=domain_dim,
            shared_dim=shared_dim,
            hidden_dim=hidden_dim,
            time_dim=time_dim,
            condition_dim=0,
            num_steps=num_steps,
            dropout=dropout,
        )


class TargetDiffusion(LatentDomainDiffusion):
    def __init__(
        self,
        domain_dim: int,
        shared_dim: int = 64,
        hidden_dim: int = 256,
        time_dim: int = 128,
        condition_dim: int | None = None,
        num_steps: int = 50,
        dropout: float = 0.0,
    ):
        super().__init__(
            domain_dim=domain_dim,
            shared_dim=shared_dim,
            hidden_dim=hidden_dim,
            time_dim=time_dim,
            condition_dim=shared_dim if condition_dim is None else condition_dim,
            num_steps=num_steps,
            dropout=dropout,
        )


class DiffusionBridge(nn.Module):
    def __init__(
        self,
        source: SourceDiffusion,
        target: TargetDiffusion,
        support_mask: nn.Module | None = None,
        support_hidden_dim: int = 256,
    ):
        super().__init__()
        self.source = source
        self.target = target
        self.support_mask = support_mask
        self.support_encoder = nn.Sequential(
            nn.Linear(source.domain_dim, support_hidden_dim),
            nn.SiLU(),
            nn.Linear(support_hidden_dim, support_hidden_dim),
            nn.SiLU(),
            nn.Linear(support_hidden_dim, target.shared_dim),
        )

    def encode_support(self, raw_x: torch.Tensor, raw_mask: torch.Tensor | None = None) -> torch.Tensor:
        masked_x = raw_x
        if self.support_mask is not None:
            masked_x = self.support_mask(masked_x, mask=raw_mask)
        return self.support_encoder(masked_x)

    def forward(self, raw_x: torch.Tensor, raw_mask: torch.Tensor | None = None) -> BridgeOutput:
        masked_x = raw_x
        if self.support_mask is not None:
            masked_x = self.support_mask(masked_x, mask=raw_mask)
        shared_latent = self.source.ddim_reverse_sample_loop(masked_x)
        support_anchor = self.encode_support(raw_x, raw_mask=raw_mask)
        target_embedding = self.target.ddim_sample_loop(shared_latent, condition=support_anchor)
        return BridgeOutput(
            shared_latent=shared_latent,
            target_embedding=target_embedding,
            support_anchor=support_anchor,
        )

    def ddim_reverse_sample_loop(self, x: torch.Tensor, raw_mask: torch.Tensor | None = None, *args, **kwargs) -> torch.Tensor:
        if self.support_mask is not None:
            x = self.support_mask(x, mask=raw_mask)
        return self.source.ddim_reverse_sample_loop(x)

    def ddim_sample_loop(
        self,
        z: torch.Tensor,
        raw_x: torch.Tensor | None = None,
        raw_mask: torch.Tensor | None = None,
        condition: torch.Tensor | None = None,
        *args,
        **kwargs,
    ) -> torch.Tensor:
        if condition is None:
            if raw_x is not None:
                condition = self.encode_support(raw_x, raw_mask=raw_mask)
            else:
                condition = torch.zeros(z.shape[0], self.target.condition_dim, device=z.device, dtype=z.dtype)
        return self.target.ddim_sample_loop(z, condition=condition)


class GaussianBridgePrior(nn.Module):
    def __init__(self, dim: int):
        super().__init__()
        self.proj = nn.Linear(dim, dim)

    def forward(self, z: torch.Tensor) -> torch.Tensor:
        return self.proj(z)
