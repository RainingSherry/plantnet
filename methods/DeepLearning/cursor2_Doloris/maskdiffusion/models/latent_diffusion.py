"""Latent-space denoising diffusion autoencoder (DDPM) for single-cell analysis.

This module implements a DDPM-style model operating in the latent space of cell embeddings:
  1. Encoder: X -> z (gene expression to latent)
  2. Denoiser: z_t -> z_0 prediction (noise removal in latent space)
  3. Decoder: z -> X_hat (latent back to gene expression)

Key design decisions:
  - Cosine noise schedule (Nichol & Dhariwal 2021)
  - Sinusoidal timestep embedding (standard in DDPM/ViT)
  - FiLM-conditioned denoiser blocks (stable and expressive)
  - Both direct and diffusion embeddings available for clustering
"""

import math
import torch
import torch.nn as nn
import torch.nn.functional as F

from .encoder import SparseEncoder
from .decoder import GeneDecoder


# ── Time embedding ──────────────────────────────────────────────────────────────


class SinusoidalTimeEmbedding(nn.Module):
    """Sinusoidal positional embedding for diffusion timesteps.

    Maps t in [0, 1] -> (dim,) vector via alternating sin/cos at geometrically
    increasing frequencies. The standard approach from "Attention is All You Need"
    and used in virtually all DDPM implementations.
    """

    def __init__(self, dim: int):
        super().__init__()
        assert dim % 2 == 0, f"time embedding dim must be even, got {dim}"
        self.dim = dim

    def forward(self, t: torch.Tensor) -> torch.Tensor:
        device = t.device
        half = self.dim // 2
        freqs = torch.exp(
            -math.log(math.pi * 2) * torch.arange(half, device=device, dtype=t.dtype) / half
        )
        args = t[:, None].to(dtype=freqs.dtype) * freqs[None, :]
        emb = torch.cat([args.sin(), args.cos()], dim=-1)
        return emb


# ── FiLM-conditioned denoiser block ─────────────────────────────────────────────


class DenoiserBlock(nn.Module):
    """One FiLM-conditioned residual block of the denoiser.

    Uses FiLM (Feature-wise Linear Modulation) conditioning from Perez et al. 2018.
    Instead of concatenating t_emb at every layer (causes shape issues), we:
      1. Apply LayerNorm to the hidden state
      2. Project t_emb to a shift and scale vector
      3. Modulate: h_norm = h_norm * (1 + shift) + bias
      4. Apply GELU -> linear -> residual add

    This is the architecture used in modern DDPM implementations (e.g., ADM, Denoising
    Diffusion Implicit Models) and avoids the residual-shape mismatch from concatenation.
    """

    def __init__(self, hidden_dim: int, time_dim: int):
        super().__init__()
        self.norm = nn.LayerNorm(hidden_dim)
        self.lin1 = nn.Linear(hidden_dim, hidden_dim)
        self.lin2 = nn.Linear(hidden_dim, hidden_dim)

        # FiLM conditioning: time_emb -> (shift, scale)
        self.film = nn.Sequential(
            nn.Linear(time_dim, hidden_dim),
            nn.Sigmoid(),  # scale in (0, 1)
        )
        self.film_bias = nn.Linear(time_dim, hidden_dim)  # shift

        self.dropout = nn.Dropout(0.1)

    def forward(self, h: torch.Tensor, t_emb: torch.Tensor) -> torch.Tensor:
        h_norm = self.norm(h)
        h_norm = h_norm * (0.5 + 0.5 * self.film(t_emb)) + self.film_bias(t_emb)
        h_norm = F.gelu(h_norm)
        h_norm = self.dropout(self.lin1(h_norm))
        h_norm = F.gelu(self.lin2(h_norm))
        return h + h_norm  # residual connection in hidden_dim space


# ── Full denoiser network ───────────────────────────────────────────────────────


class Denoiser(nn.Module):
    """FiLM-conditioned denoiser for DDPM in latent space.

    Architecture:
        [z_t, t_emb] -> Linear(latent_dim + time_dim, hidden_dim)
                     -> N x DenoiserBlock(hidden_dim, time_dim)
                     -> Linear(hidden_dim, latent_dim) -> noise prediction
    """

    def __init__(
        self,
        latent_dim: int,
        time_dim: int = 32,
        hidden_dim: int = 256,
        n_blocks: int = 3,
    ):
        super().__init__()
        self.input_proj = nn.Sequential(
            nn.Linear(latent_dim + time_dim, hidden_dim),
            nn.LayerNorm(hidden_dim),
            nn.GELU(),
        )
        self.blocks = nn.ModuleList(
            [DenoiserBlock(hidden_dim, time_dim) for _ in range(n_blocks)]
        )
        self.output_proj = nn.Linear(hidden_dim, latent_dim)

    def forward(self, z_t: torch.Tensor, t_emb: torch.Tensor) -> torch.Tensor:
        h = self.input_proj(torch.cat([z_t, t_emb], dim=-1))
        for block in self.blocks:
            h = block(h, t_emb)
        return self.output_proj(h)


# ── Main LatentDiffusionAE ─────────────────────────────────────────────────────


class LatentDiffusionAE(nn.Module):
    """Latent-space denoising diffusion autoencoder.

    Data flow:
        X -> Encoder -> z (clean latent)
             |
             +-> q_sample(z, t) = sqrt(α_bar_t)·z + sqrt(1-α_bar_t)·ε  (forward)
             |
             +-> Denoiser(z_t, t) = ε̂  (predict noise, for training loss)
             |
             +-> p_sample_loop(z_t)  (reverse, for inference)
             |
        Decoder(z_denoised) -> X̂  (reconstruction)

    Both direct embedding (no diffusion) AND diffusion embedding are available.
    """

    def __init__(
        self,
        n_genes: int,
        latent_dim: int = 32,
        hidden_dim: int = 256,
        diffusion_hidden_dim: int = 256,
        diffusion_steps: int = 100,
        dropout: float = 0.1,
    ):
        super().__init__()
        self.latent_dim = latent_dim
        self.diffusion_steps = diffusion_steps

        self.encoder = SparseEncoder(
            n_genes=n_genes,
            latent_dim=latent_dim,
            hidden_dim=hidden_dim,
            dropout=dropout,
        )
        self.decoder = GeneDecoder(
            latent_dim=latent_dim,
            n_genes=n_genes,
            hidden_dim=hidden_dim,
            dropout=dropout,
        )

        # Time embedding dimension
        time_dim = max(latent_dim // 2, 32)
        self.denoiser = Denoiser(
            latent_dim=latent_dim,
            time_dim=time_dim,
            hidden_dim=diffusion_hidden_dim,
            n_blocks=3,
        )
        self.time_embed = SinusoidalTimeEmbedding(time_dim)

        self.register_buffer("_dummy", torch.tensor(0.0))
        self._init_diffusion_schedule()

    def _init_diffusion_schedule(self):
        """Compute cosine schedule (Nichol & Dhariwal 2021).
        α_bar[t] = cos²((t/T + s) / (1+s) · π/2)  with s = 0.008
        β[t] = 1 - α[t] / α_bar[t-1]  clipped to [1e-4, 0.999]
        """
        T = self.diffusion_steps
        s = 0.008
        t = torch.linspace(0, T, T + 1)
        sched = torch.cos(((t / T + s) / (1 + s)) * (math.pi / 2)) ** 2
        alphas_cumprod = sched / sched[0]
        alphas = alphas_cumprod[1:] / alphas_cumprod[:-1]
        betas = (1 - alphas).clamp(1e-4, 0.999).cpu()

        self.register_buffer("betas", betas, persistent=False)
        self.register_buffer("alphas_cumprod", alphas_cumprod[:-1], persistent=False)
        self.register_buffer(
            "sqrt_alphas_cumprod", alphas_cumprod[:-1].sqrt(), persistent=False
        )
        self.register_buffer(
            "sqrt_one_minus_alphas_cumprod", (1 - alphas_cumprod[:-1]).sqrt(), persistent=False
        )
        self.register_buffer("sqrt_recip_alphas", (1 / alphas).sqrt(), persistent=False)

    @property
    def device(self) -> torch.device:
        return self._dummy.device

    def encode(self, x: torch.Tensor) -> torch.Tensor:
        """Encode input to latent space, clamped to prevent extreme values."""
        z = self.encoder(x)
        return torch.clamp(z, -10.0, 10.0)

    def decode(self, z: torch.Tensor) -> torch.Tensor:
        """Decode latent to gene expression reconstruction."""
        return self.decoder(z)

    # ── Forward diffusion (corrupt z0 -> zt) ───────────────────────────────────

    def q_sample(self, z0: torch.Tensor, t: torch.Tensor, noise: torch.Tensor = None) -> tuple:
        """Forward diffusion: corrupt clean z0 at timestep t.
        z_t = sqrt(α_bar_t) · z_0 + sqrt(1 - α_bar_t) · ε
        Returns (z_t, noise).
        """
        if noise is None:
            noise = torch.randn_like(z0)
        sqrt_abar = self._interp(self.sqrt_alphas_cumprod, t)[:, None]
        sqrt_oabar = self._interp(self.sqrt_one_minus_alphas_cumprod, t)[:, None]
        z_t = sqrt_abar * z0 + sqrt_oabar * noise
        return z_t, noise

    @staticmethod
    def _interp(buf: torch.Tensor, t: torch.Tensor) -> torch.Tensor:
        """Linearly interpolate schedule array at fractional timestep t."""
        T = len(buf) - 1
        t_clipped = t.clamp(0.0, 1.0 - 1e-6)
        idx = t_clipped * T
        lo = idx.long()
        hi = (lo + 1).clamp(max=len(buf) - 1)
        frac = idx - lo.float()
        return buf[lo] * (1 - frac) + buf[hi] * frac

    # ── Reverse process ────────────────────────────────────────────────────────

    def predict_noise(self, z_t: torch.Tensor, t: torch.Tensor) -> torch.Tensor:
        """Predict noise ε_θ(z_t, t) using the denoiser network."""
        t_emb = self.time_embed(t)
        return self.denoiser(z_t, t_emb)

    def p_sample(self, z_t: torch.Tensor, t: torch.Tensor) -> torch.Tensor:
        """Single reverse diffusion step (DDPM formula).

        Given z_t and timestep t, predicts z_{t-1}:
            z_{t-1} = μ_θ(z_t, t) + σ_t · ε,  ε ~ N(0, I)

        where μ_θ(z_t, t) = (z_t - √(β_t/ᾱ_t)·ε̂) / √(ᾱ_{t-1})
        and ε̂ = Denoiser(z_t, t).
        """
        noise_pred = self.predict_noise(z_t, t)
        t_next = t - 1.0 / (self.diffusion_steps - 1)
        t_next = torch.clamp(t_next, 0.0, 1.0)

        sqrt_abar = self._interp(self.sqrt_alphas_cumprod, t)[:, None]
        sqrt_abar_next = self._interp(self.sqrt_alphas_cumprod, t_next)[:, None]

        alpha_t = sqrt_abar ** 2
        alpha_t_next = sqrt_abar_next ** 2
        beta_t = (1 - alpha_t / (alpha_t_next + 1e-8)).clamp(1e-4, 0.999)

        # DDPM mean: μ = (z_t - √(β_t/ᾱ_t)·ε̂) / √(ᾱ_{t-1})
        mean = (z_t - beta_t.sqrt() * noise_pred / (alpha_t.sqrt() + 1e-8)) / (
            sqrt_abar_next + 1e-8
        )
        z_prev = mean + beta_t.sqrt() * torch.randn_like(z_t)
        return torch.clamp(z_prev, -10.0, 10.0)

    @torch.no_grad()
    def p_sample_loop(self, z_t: torch.Tensor, n_steps: int = None) -> torch.Tensor:
        """Full reverse diffusion: z_T -> z_0 by running n_steps reverse steps."""
        n_steps = n_steps or self.diffusion_steps
        z = z_t.clone()
        for step in reversed(range(n_steps)):
            t_val = step / max(n_steps - 1, 1)
            t = torch.full((z.shape[0],), t_val, device=z.device, dtype=z.dtype)
            z = self.p_sample(z, t)
        return z

    @torch.no_grad()
    def denoise_embedding(self, z0: torch.Tensor, start_frac: float = 0.35) -> torch.Tensor:
        """Refine an encoded latent from a partial noising point.

        The previous cursor2 implementation fed a clean encoder latent directly
        into the full reverse chain, which treats it like pure z_T and often
        destroys cluster geometry. This deterministic DDIM-like refinement keeps
        the cell identity anchored in z0 while still testing whether the denoiser
        improves the latent space.
        """
        n_steps = max(2, int(round(self.diffusion_steps * start_frac)))
        start_t = torch.full((z0.shape[0],), start_frac, device=z0.device, dtype=z0.dtype)
        z, _ = self.q_sample(z0, start_t, noise=torch.zeros_like(z0))
        for step in reversed(range(n_steps)):
            t_val = step / max(self.diffusion_steps - 1, 1)
            t = torch.full((z.shape[0],), t_val, device=z.device, dtype=z.dtype)
            noise_pred = self.predict_noise(z, t)
            sqrt_abar = self._interp(self.sqrt_alphas_cumprod, t)[:, None]
            sqrt_oabar = self._interp(self.sqrt_one_minus_alphas_cumprod, t)[:, None]
            pred_x0 = (z - sqrt_oabar * noise_pred) / (sqrt_abar + 1e-8)
            if step == 0:
                z = pred_x0
            else:
                t_prev = torch.full(
                    (z.shape[0],),
                    (step - 1) / max(self.diffusion_steps - 1, 1),
                    device=z.device,
                    dtype=z.dtype,
                )
                sqrt_abar_prev = self._interp(self.sqrt_alphas_cumprod, t_prev)[:, None]
                sqrt_oabar_prev = self._interp(self.sqrt_one_minus_alphas_cumprod, t_prev)[:, None]
                z = sqrt_abar_prev * pred_x0 + sqrt_oabar_prev * noise_pred
            z = torch.clamp(z, -10.0, 10.0)
        return z

    # ── Training loss ─────────────────────────────────────────────────────────

    def diffusion_loss(self, z0: torch.Tensor) -> torch.Tensor:
        """DDPM training objective: MSE between predicted noise and true noise."""
        batch = z0.shape[0]
        t = torch.rand(batch, device=z0.device, dtype=z0.dtype)
        z_t, noise = self.q_sample(z0, t)
        noise_pred = self.predict_noise(z_t, t)
        return torch.mean((noise_pred - noise) ** 2)

    def recon_loss(
        self,
        x: torch.Tensor,
        z: torch.Tensor,
        mask: torch.Tensor = None,
    ) -> torch.Tensor:
        """Reconstruction loss (MSE) on gene expression.
        If mask is provided, only compute loss on masked (active) positions.
        mask shape: (batch, n_genes), values in [0, 1].
        """
        x_hat = self.decode(z)
        if mask is not None:
            diff = (x_hat - x) ** 2
            loss = (diff * mask).sum() / (mask.sum() + 1e-8)
        else:
            loss = torch.mean((x_hat - x) ** 2)
        return loss

    # ── Full forward pass ───────────────────────────────────────────────────

    def forward(
        self,
        x: torch.Tensor,
        mask: torch.Tensor = None,
        return_recon: bool = True,
        sample_diffusion: bool = False,
    ) -> dict:
        """Full forward pass through encoder + (optional) diffusion + decoder.

        Args:
            x: Input gene expression, shape (batch, n_genes).
            mask: Soft support mask (batch, n_genes), values in [0, 1].
            return_recon: If True, also return reconstructed gene expression.
            sample_diffusion: If True, run full diffusion sampling for z_denoised.
        Returns:
            dict with: z, z_denoised, x_recon, diffusion_loss, recon_loss
        """
        z = self.encode(x)
        z_denoised = z.clone()

        # Diffusion loss on latent space (always computed during training)
        loss_diffusion = self.diffusion_loss(z)

        # Reverse diffusion (only at inference when sample_diffusion=True)
        if sample_diffusion:
            z_denoised = self.p_sample_loop(z)

        # Reconstruction
        x_recon = None
        loss_recon = torch.tensor(0.0, device=x.device)
        if return_recon:
            x_recon = self.decode(z_denoised)
            loss_recon = self.recon_loss(x, z_denoised, mask)

        return {
            "z": z,
            "z_denoised": z_denoised,
            "x_recon": x_recon,
            "diffusion_loss": loss_diffusion,
            "recon_loss": loss_recon,
        }

    # ── Embedding extraction ──────────────────────────────────────────────────

    def get_direct_embedding(self, x: torch.Tensor) -> torch.Tensor:
        """Return raw encoder output (no diffusion sampling)."""
        z = self.encode(x)
        return torch.clamp(z, -10.0, 10.0)

    @torch.no_grad()
    def get_diffusion_embedding(self, x: torch.Tensor) -> torch.Tensor:
        """Return diffusion-sampled embedding (full reverse process)."""
        z = self.encode(x)
        return self.denoise_embedding(z)
