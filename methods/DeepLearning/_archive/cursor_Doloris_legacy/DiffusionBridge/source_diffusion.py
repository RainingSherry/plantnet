"""
Dual-Domain Diffusion Bridge for Single-Cell Clustering

Core framework: Raw Sparse Count Domain -> Shared Gaussian Latent -> Cluster-Separable Denoised Domain

This module implements the source diffusion model that learns the raw sparse expression
distribution. Inspired by DOLORIS's DDIB (Dual Diffusion Implicit Bridges) paradigm.
"""
from __future__ import annotations

import math

import torch
import torch.nn as nn
import torch.nn.functional as F
from torch import nn

from nn import timestep_embedding


def _extract(coeff: torch.Tensor, t: torch.Tensor, x_shape: torch.Size) -> torch.Tensor:
    values = coeff.gather(0, t)
    return values.view(-1, *([1] * (len(x_shape) - 1)))


class MLPBlock(nn.Module):
    """MLP block with optional LayerNorm and SiLU activation."""

    def __init__(self, dim: int, depth: int = 3, dropout: float = 0.0, act: str = "SiLU"):
        super().__init__()
        layers = []
        for i in range(depth):
            layers.append(nn.Linear(dim, dim))
            if i < depth - 1:
                layers.append(nn.LayerNorm(dim) if act == "SiLU" else nn.BatchNorm1d(dim))
                layers.append(nn.SiLU() if act == "SiLU" else nn.ReLU())
                if dropout > 0:
                    layers.append(nn.Dropout(dropout))
        self.net = nn.Sequential(*layers)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.net(x)


class DenoiserBlock(nn.Module):
    """
    Time-conditioned denoiser block.
    Takes noisy latent x_t, timestep t, and optional condition, outputs prediction.
    """

    def __init__(
        self,
        input_dim: int,
        output_dim: int,
        time_embed_dim: int,
        cond_dim: int = 0,
        depth: int = 3,
        dropout: float = 0.0,
    ):
        super().__init__()
        self.time_mlp = nn.Sequential(
            nn.Linear(time_embed_dim, time_embed_dim * 2),
            nn.SiLU(),
            nn.Linear(time_embed_dim * 2, time_embed_dim),
            nn.SiLU(),
        )
        self.input_proj = nn.Linear(input_dim, time_embed_dim)
        self.cond_proj = nn.Linear(cond_dim, time_embed_dim) if cond_dim > 0 else None

        self.blocks = nn.ModuleList()
        for _ in range(depth):
            self.blocks.append(
                nn.Sequential(
                    nn.Linear(time_embed_dim, time_embed_dim),
                    nn.LayerNorm(time_embed_dim),
                    nn.SiLU(),
                    nn.Dropout(dropout),
                    nn.Linear(time_embed_dim, time_embed_dim),
                )
            )
        self.out = nn.Linear(time_embed_dim, output_dim)

    def forward(
        self, x_t: torch.Tensor, t_emb: torch.Tensor, cond: torch.Tensor | None = None
    ) -> torch.Tensor:
        h = self.input_proj(x_t) + self.time_mlp(t_emb)
        if self.cond_proj is not None and cond is not None:
            h = h + self.cond_proj(cond)
        for block in self.blocks:
            h = h + block(h)
        return self.out(h)


class LatentDomainDiffusion(nn.Module):
    """
    Base class for domain-specific diffusion models in the bridge framework.

    Learns to denoise samples from a domain (raw expression or target embedding)
    toward/away from a shared Gaussian latent space.
    """

    def __init__(
        self,
        domain_dim: int,
        shared_dim: int,
        hidden_dim: int = 256,
        time_embed_dim: int = 128,
        cond_dim: int = 0,
        num_steps: int = 50,
        beta_start: float = 1e-4,
        beta_end: float = 2e-2,
        dropout: float = 0.0,
    ):
        super().__init__()
        self.domain_dim = domain_dim
        self.shared_dim = shared_dim
        self.cond_dim = cond_dim
        self.num_steps = num_steps

        # Encode domain input to shared latent dimension
        self.domain_encoder = nn.Sequential(
            nn.Linear(domain_dim, hidden_dim),
            nn.LayerNorm(hidden_dim),
            nn.SiLU(),
            nn.Linear(hidden_dim, hidden_dim),
            nn.LayerNorm(hidden_dim),
            nn.SiLU(),
            nn.Linear(hidden_dim, shared_dim),
        )

        # Decode shared latent back to domain dimension
        self.domain_decoder = nn.Sequential(
            nn.Linear(shared_dim, hidden_dim),
            nn.LayerNorm(hidden_dim),
            nn.SiLU(),
            nn.Linear(hidden_dim, hidden_dim),
            nn.LayerNorm(hidden_dim),
            nn.SiLU(),
            nn.Linear(hidden_dim, domain_dim),
        )

        # Time-conditioned denoiser: predicts x_0 from (x_t, t)
        self.denoiser = DenoiserBlock(
            input_dim=shared_dim,
            output_dim=shared_dim,
            time_embed_dim=time_embed_dim,
            cond_dim=cond_dim,
            depth=3,
            dropout=dropout,
        )

        # Precompute diffusion schedule
        betas = torch.linspace(beta_start, beta_end, num_steps, dtype=torch.float32)
        alphas = 1.0 - betas
        alphas_cumprod = torch.cumprod(alphas, dim=0)
        alphas_cumprod_prev = torch.cat([torch.ones(1, dtype=torch.float32), alphas_cumprod[:-1]])
        alphas_cumprod_next = torch.cat([alphas_cumprod[1:], torch.zeros(1, dtype=torch.float32)])

        self.register_buffer("betas", betas)
        self.register_buffer("alphas_cumprod", alphas_cumprod)
        self.register_buffer("alphas_cumprod_prev", alphas_cumprod_prev)
        self.register_buffer("alphas_cumprod_next", alphas_cumprod_next)
        self.register_buffer("sqrt_alphas_cumprod", torch.sqrt(alphas_cumprod))
        self.register_buffer(
            "sqrt_one_minus_alphas_cumprod", torch.sqrt(1.0 - alphas_cumprod)
        )
        self.register_buffer(
            "sqrt_recip_alphas_cumprod", torch.sqrt(1.0 / alphas_cumprod.clamp_min(1e-8))
        )
        self.register_buffer(
            "sqrt_recipm1_alphas_cumprod",
            torch.sqrt((1.0 / alphas_cumprod.clamp_min(1e-8)) - 1.0),
        )

    def encode(self, x: torch.Tensor) -> torch.Tensor:
        """Project domain input to shared latent."""
        return self.domain_encoder(x)

    def decode(self, z: torch.Tensor) -> torch.Tensor:
        """Project shared latent back to domain."""
        return self.domain_decoder(z)

    def q_sample(self, x_start: torch.Tensor, t: torch.Tensor, noise: torch.Tensor | None = None) -> torch.Tensor:
        """Forward diffusion: add noise to x_start at timestep t."""
        if noise is None:
            noise = torch.randn_like(x_start)
        return (
            _extract(self.sqrt_alphas_cumprod, t, x_start.shape) * x_start
            + _extract(self.sqrt_one_minus_alphas_cumprod, t, x_start.shape) * noise
        )

    def predict_x0(self, x_t: torch.Tensor, t: torch.Tensor, eps: torch.Tensor) -> torch.Tensor:
        """Predict clean sample x_0 from noise prediction."""
        return (
            _extract(self.sqrt_recip_alphas_cumprod, t, x_t.shape) * x_t
            - _extract(self.sqrt_recipm1_alphas_cumprod, t, x_t.shape) * eps
        )

    def predict_noise(self, x_t: torch.Tensor, t: torch.Tensor, cond: torch.Tensor | None = None) -> torch.Tensor:
        """Predict noise at timestep t."""
        t_emb = timestep_embedding(t, self.denoiser.time_mlp[0].in_features)
        return self.denoiser(x_t, t_emb, cond=cond)

    def denoise_step(
        self, x_t: torch.Tensor, t: torch.Tensor, cond: torch.Tensor | None = None
    ) -> torch.Tensor:
        """One step of denoising: predict x_0 then reconstruct."""
        eps = self.predict_noise(x_t, t, cond=cond)
        return self.predict_x0(x_t, t, eps)

    def training_loss(
        self,
        domain_x: torch.Tensor,
        cond: torch.Tensor | None = None,
        recon_weight: float = 1.0,
        prior_weight: float = 1e-3,
        zero_weight: float = 0.25,
    ) -> dict[str, torch.Tensor]:
        """
        Compute training loss for the diffusion model.

        L = L_diffusion + L_reconstruction + L_prior

        - L_diffusion: MSE between predicted and true noise
        - L_reconstruction: weighted MSE of decoded domain vs original domain
        - L_prior: encourages encoded latent to be near isotropic Gaussian
        """
        shared_seed = self.encode(domain_x)
        B = domain_x.shape[0]
        t = torch.randint(0, self.num_steps, (B,), device=domain_x.device, dtype=torch.long)
        noise = torch.randn_like(shared_seed)
        x_t = self.q_sample(shared_seed, t, noise=noise)

        t_emb = timestep_embedding(t, self.denoiser.time_mlp[0].in_features)
        pred_eps = self.denoiser(x_t, t_emb, cond=cond)
        pred_x0 = self.predict_x0(x_t, t, pred_eps)
        pred_domain = self.decode(pred_x0)

        diffusion_loss = F.mse_loss(pred_eps, noise)

        if domain_x.shape[-1] == pred_domain.shape[-1]:
            weights = torch.where(
                domain_x > 0,
                torch.ones_like(domain_x),
                torch.full_like(domain_x, zero_weight),
            )
            recon_loss = ((pred_domain - domain_x).pow(2) * weights).mean()
        else:
            recon_loss = F.mse_loss(pred_domain, domain_x)

        prior_loss = shared_seed.pow(2).mean()

        loss = diffusion_loss + recon_weight * recon_loss + prior_weight * prior_loss
        return {
            "loss": loss,
            "diffusion_loss": diffusion_loss.detach(),
            "recon_loss": recon_loss.detach(),
            "prior_loss": prior_loss.detach(),
        }

    # -------------------------------------------------------------------------
    # DDIM inference (used by the bridge)
    # -------------------------------------------------------------------------

    def ddim_reverse_step(
        self, x_t: torch.Tensor, t: torch.Tensor, cond: torch.Tensor | None = None
    ) -> torch.Tensor:
        """
        DDIM reverse step (DDIM inversion).
        Inverts x_t one step toward the Gaussian prior.
        This is used by the source model to encode raw cells into the shared latent.
        """
        pred_x0 = self.denoise_step(x_t, t, cond=cond)
        eps = (x_t - _extract(self.sqrt_alphas_cumprod, t, x_t.shape) * pred_x0) / _extract(
            self.sqrt_one_minus_alphas_cumprod, t, x_t.shape
        ).clamp_min(1e-8)
        alpha_bar_next = _extract(self.alphas_cumprod_next, t, x_t.shape)
        return pred_x0 * alpha_bar_next.sqrt() + (1.0 - alpha_bar_next).clamp_min(0.0).sqrt() * eps

    def ddim_sample_step(
        self, x_t: torch.Tensor, t: torch.Tensor, cond: torch.Tensor | None = None
    ) -> torch.Tensor:
        """
        DDIM sample step (forward denoising).
        Denoises x_t one step toward the domain distribution.
        This is used by the target model to decode from the shared latent.
        """
        pred_x0 = self.denoise_step(x_t, t, cond=cond)
        eps = (x_t - _extract(self.sqrt_alphas_cumprod, t, x_t.shape) * pred_x0) / _extract(
            self.sqrt_one_minus_alphas_cumprod, t, x_t.shape
        ).clamp_min(1e-8)
        alpha_bar_prev = _extract(self.alphas_cumprod_prev, t, x_t.shape)
        return pred_x0 * alpha_bar_prev.sqrt() + (1.0 - alpha_bar_prev).clamp_min(0.0).sqrt() * eps

    def ddim_reverse_sample_loop(
        self, image: torch.Tensor, cond: torch.Tensor | None = None
    ) -> torch.Tensor:
        """
        Full DDIM reverse trajectory: raw domain -> shared Gaussian latent.
        Encodes a raw sample into the latent space via reverse ODE.
        """
        state = self.encode(image)
        for step in range(self.num_steps):
            t = torch.full((image.shape[0],), step, device=image.device, dtype=torch.long)
            state = self.ddim_reverse_step(state, t, cond=cond)
        return state

    def ddim_sample_loop(
        self, noise: torch.Tensor, cond: torch.Tensor | None = None
    ) -> torch.Tensor:
        """
        Full DDIM sample trajectory: shared Gaussian latent -> domain.
        Decodes from the latent space to the target domain via forward ODE.
        """
        state = noise
        for step in reversed(range(self.num_steps)):
            t = torch.full((noise.shape[0],), step, device=noise.device, dtype=torch.long)
            state = self.ddim_sample_step(state, t, cond=cond)
        return self.decode(state)


class SourceDiffusion(LatentDomainDiffusion):
    """
    Source diffusion model: learns raw sparse expression distribution.

    This model takes raw gene expression vectors (log1p-normalized counts)
    and learns to project them into the shared Gaussian latent space via
    DDIM inversion. The latent should capture the essential structure of
    raw sparse data while being isotropically Gaussian.
    """

    def __init__(
        self,
        domain_dim: int,
        shared_dim: int = 64,
        hidden_dim: int = 256,
        time_embed_dim: int = 128,
        num_steps: int = 50,
        dropout: float = 0.0,
    ):
        super().__init__(
            domain_dim=domain_dim,
            shared_dim=shared_dim,
            hidden_dim=hidden_dim,
            time_embed_dim=time_embed_dim,
            cond_dim=0,
            num_steps=num_steps,
            dropout=dropout,
        )


class TargetDiffusion(LatentDomainDiffusion):
    """
    Target diffusion model: learns cluster-friendly denoised representation.

    This model takes denoised target embeddings (from teacher methods like
    PCA+graph smoothing or PhytoCluster) and learns to sample from the
    cluster-separable manifold. It uses the shared latent as input and
    optionally conditions on support anchors.
    """

    def __init__(
        self,
        domain_dim: int,
        shared_dim: int = 64,
        hidden_dim: int = 256,
        time_embed_dim: int = 128,
        cond_dim: int | None = None,
        num_steps: int = 50,
        dropout: float = 0.0,
    ):
        super().__init__(
            domain_dim=domain_dim,
            shared_dim=shared_dim,
            hidden_dim=hidden_dim,
            time_embed_dim=time_embed_dim,
            cond_dim=shared_dim if cond_dim is None else cond_dim,
            num_steps=num_steps,
            dropout=dropout,
        )
