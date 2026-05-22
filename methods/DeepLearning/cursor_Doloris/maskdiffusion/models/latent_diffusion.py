# -*- coding: utf-8 -*-
"""
LatentDiffusionAE: Latent Space Diffusion Autoencoder for scRNA-seq Clustering

核心思想：
1. Encoder 将高维表达 X 映射到低维潜空间 z
2. Diffusion prior 在 z 上进行去噪扩散
3. Decoder 从 z 重构表达值（仅在 active genes 上）
4. 与 SupportMaskNet 配合：mask 指导 decoder 只关注激活基因

关键设计：
- 在低维空间做扩散（而非高维基因空间），大幅降低计算成本
- 结合 DOLORIS 的 sparsity masking：decoder 输出只在激活基因上有意义
- 条件信息可以是 cell type embedding 或 cluster embedding
"""

import torch
import torch.nn as nn
import torch.nn.functional as F
import math
import numpy as np


def timestep_embedding(timesteps: torch.Tensor, dim: int) -> torch.Tensor:
    """
    Create sinusoidal timestep embeddings.

    Args:
        timesteps: 1D tensor of timesteps
        dim: dimension of the embedding

    Returns:
        (len(timesteps), dim) embedding tensor
    """
    half_dim = dim // 2
    emb = math.log(10000) / (half_dim - 1)
    emb = torch.exp(torch.arange(half_dim, device=timesteps.device) * -emb)
    emb = timesteps[:, None].float() * emb[None, :]
    emb = torch.cat([torch.sin(emb), torch.cos(emb)], dim=-1)
    return emb


class DiffusionBlock(nn.Module):
    """
    Single diffusion transformation block with time conditioning.

    Architecture:
        Input + TimeEmb → MLP → Residual Connection
    """

    def __init__(
        self,
        hidden_dim: int,
        time_embed_dim: int,
        dropout: float = 0.1,
    ):
        super().__init__()
        self.time_mlp = nn.Sequential(
            nn.Linear(time_embed_dim, hidden_dim),
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

        self.residual_scale = nn.Parameter(torch.tensor(0.0))

    def forward(self, x: torch.Tensor, time_emb: torch.Tensor) -> torch.Tensor:
        """
        Args:
            x: (batch, hidden_dim)
            time_emb: (batch, time_embed_dim)
        """
        h = x + self.residual_scale * self.time_mlp(time_emb)
        h = h + self.net(h)
        return h


class ConditionalMLP(nn.Module):
    """
    MLP with conditioning on auxiliary information (e.g., cluster labels, cell type).

    Architecture:
        concat([x, cond]) → Linear → layers → output
    """

    def __init__(
        self,
        input_dim: int,
        output_dim: int,
        cond_dim: int,
        hidden_dims: list = None,
        dropout: float = 0.1,
    ):
        super().__init__()
        if hidden_dims is None:
            hidden_dims = [256, 128]

        layers = []
        in_dim = input_dim + cond_dim
        for hidden_dim in hidden_dims:
            layers.extend([
                nn.Linear(in_dim, hidden_dim),
                nn.LayerNorm(hidden_dim),
                nn.Mish(inplace=True),
                nn.Dropout(dropout),
            ])
            in_dim = hidden_dim

        layers.append(nn.Linear(in_dim, output_dim))
        self.net = nn.Sequential(*layers)

    def forward(self, x: torch.Tensor, cond: torch.Tensor) -> torch.Tensor:
        """
        Args:
            x: (batch, input_dim)
            cond: (batch, cond_dim)
        """
        return self.net(torch.cat([x, cond], dim=-1))


class LatentDenoiser(nn.Module):
    """
    Denoiser network for latent space diffusion.

    Predicts x0 from noisy latent z_t and timestep t.
    """

    def __init__(
        self,
        latent_dim: int,
        hidden_dim: int = 256,
        num_layers: int = 3,
        time_embed_dim: int = 128,
        dropout: float = 0.1,
        cond_dim: int = 0,
    ):
        super().__init__()
        self.latent_dim = latent_dim
        self.time_embed_dim = time_embed_dim

        # Time embedding
        self.time_embed = nn.Sequential(
            nn.Linear(time_embed_dim, hidden_dim),
            nn.Mish(inplace=True),
            nn.Linear(hidden_dim, hidden_dim),
        )

        # Input projection
        self.input_proj = nn.Linear(latent_dim, hidden_dim)

        # Diffusion blocks
        self.blocks = nn.ModuleList([
            DiffusionBlock(hidden_dim, time_embed_dim, dropout)
            for _ in range(num_layers)
        ])

        # Output projection
        self.output_proj = nn.Linear(hidden_dim, latent_dim)

        # Optional conditioning
        self.has_cond = cond_dim > 0
        if self.has_cond:
            self.cond_embed = nn.Linear(cond_dim, hidden_dim)

    def forward(
        self,
        z_t: torch.Tensor,
        timesteps: torch.Tensor,
        cond: torch.Tensor = None,
    ) -> torch.Tensor:
        """
        Args:
            z_t: Noisy latent at timestep t (batch, latent_dim)
            timesteps: Timestep values (batch,)
            cond: Optional conditioning (batch, cond_dim)

        Returns:
            Predicted clean latent x0 (batch, latent_dim)
        """
        # Time embedding
        t_emb = timestep_embedding(timesteps, self.time_embed_dim)
        t_emb = self.time_embed(t_emb)  # (batch, hidden_dim)

        # Input projection
        h = self.input_proj(z_t)  # (batch, hidden_dim)

        # Add conditioning
        if self.has_cond and cond is not None:
            cond_emb = self.cond_embed(cond)
            h = h + cond_emb

        # Apply diffusion blocks
        for block in self.blocks:
            h = block(h, t_emb)

        # Output
        x0_pred = self.output_proj(h)
        return x0_pred


class LatentDiffusionAE(nn.Module):
    """
    Latent Space Diffusion Autoencoder for scRNA-seq clustering.

    Pipeline:
        X → Encoder → z → Diffuse → z_denoised → Decoder → X_hat

    Losses:
        1. Diffusion loss: MSE between predicted x0 and true z
        2. Reconstruction loss: MSE between decoded X and input X (masked)
        3. Optional: Cluster alignment loss

    Args:
        num_genes: Number of genes (input dimension)
        latent_dim: Dimension of latent space
        hidden_dims: Encoder hidden layer dimensions
        diffusion_steps: Number of diffusion timesteps
        diffusion_type: 'ddpm' or 'ddim'
    """

    def __init__(
        self,
        num_genes: int,
        latent_dim: int = 32,
        hidden_dims: list = None,
        diffusion_steps: int = 100,
        diffusion_type: str = 'ddpm',
        dropout: float = 0.1,
        cond_dim: int = 0,
    ):
        super().__init__()
        if hidden_dims is None:
            hidden_dims = [512, 256]

        self.num_genes = num_genes
        self.latent_dim = latent_dim
        self.diffusion_steps = diffusion_steps
        self.diffusion_type = diffusion_type

        # Encoder: X → z
        encoder_layers = []
        in_dim = num_genes
        for hidden_dim in hidden_dims:
            encoder_layers.extend([
                nn.Linear(in_dim, hidden_dim),
                nn.LayerNorm(hidden_dim),
                nn.Mish(inplace=True),
                nn.Dropout(dropout),
            ])
            in_dim = hidden_dim
        encoder_layers.append(nn.Linear(in_dim, latent_dim))
        self.encoder = nn.Sequential(*encoder_layers)

        # Decoder: z → X
        decoder_layers = []
        in_dim = latent_dim
        for hidden_dim in hidden_dims[::-1]:
            decoder_layers.extend([
                nn.Linear(in_dim, hidden_dim),
                nn.LayerNorm(hidden_dim),
                nn.Mish(inplace=True),
                nn.Dropout(dropout),
            ])
            in_dim = hidden_dim
        decoder_layers.append(nn.Linear(in_dim, num_genes))
        self.decoder = nn.Sequential(*decoder_layers)

        # Latent denoiser
        self.denoiser = LatentDenoiser(
            latent_dim=latent_dim,
            hidden_dim=max(latent_dim * 2, 128),
            num_layers=3,
            time_embed_dim=128,
            dropout=dropout,
            cond_dim=cond_dim,
        )

        # Diffusion schedule
        self.register_buffer(
            'betas',
            self._cosine_beta_schedule(diffusion_steps),
        )
        alphas = 1.0 - self.betas
        self.register_buffer('alphas_cumprod', torch.cumprod(alphas, dim=0))
        self.register_buffer('sqrt_alphas_cumprod', torch.sqrt(self.alphas_cumprod))
        self.register_buffer('sqrt_one_minus_alphas_cumprod', torch.sqrt(1.0 - self.alphas_cumprod))

    def _cosine_beta_schedule(self, timesteps: int, s: float = 0.008) -> torch.Tensor:
        """
        Cosine schedule as proposed in https://arxiv.org/abs/2102.09672.
        """
        steps = timesteps + 1
        x = torch.linspace(0, timesteps, steps)
        alphas_cumprod = torch.cos(((x / timesteps) + s) / (1 + s) * torch.pi * 0.5) ** 2
        alphas_cumprod = alphas_cumprod / alphas_cumprod[0]
        betas = 1 - (alphas_cumprod[1:] / alphas_cumprod[:-1])
        return torch.clip(betas, 0.0001, 0.02)

    def q_sample(self, z0: torch.Tensor, t: torch.Tensor, noise: torch.Tensor = None) -> torch.Tensor:
        """
        Forward diffusion: add noise to z0.

        Args:
            z0: Clean latent (batch, latent_dim)
            t: Timestep (batch,)
            noise: Optional noise tensor

        Returns:
            Noisy latent z_t (batch, latent_dim)
        """
        if noise is None:
            noise = torch.randn_like(z0)

        sqrt_alphas_cumprod_t = self._extract(self.sqrt_alphas_cumprod, t, z0.shape)
        sqrt_one_minus_alphas_cumprod_t = self._extract(self.sqrt_one_minus_alphas_cumprod, t, z0.shape)

        return sqrt_alphas_cumprod_t * z0 + sqrt_one_minus_alphas_cumprod_t * noise

    def _extract(self, a: torch.Tensor, t: torch.Tensor, x_shape: tuple) -> torch.Tensor:
        """Extract values from a at indices t and reshape for broadcasting."""
        batch_size = t.shape[0]
        out = a.to(t.device).gather(0, t)
        return out.view(batch_size, *((1,) * (len(x_shape) - 1)))

    def p_mean_variance(
        self,
        z_t: torch.Tensor,
        t: torch.Tensor,
        cond: torch.Tensor = None,
    ) -> tuple:
        """
        Compute mean and variance for reverse step.

        Args:
            z_t: Noisy latent (batch, latent_dim)
            t: Timestep (batch,)
            cond: Optional conditioning

        Returns:
            (mean, variance)
        """
        # Predict x0
        x0_pred = self.denoiser(z_t, t, cond)

        # Compute posterior mean
        t_prev = (t - 1).clamp(0, self.diffusion_steps - 1)
        
        # Safe extraction for alpha_t - handle dimension mismatch
        alpha_t = self._extract(self.alphas_cumprod, t, z_t.shape)
        alpha_t_prev = self._extract(self.alphas_cumprod, t_prev, z_t.shape)
        
        # Compute posterior mean with numerical stability
        sqrt_alphas_cumprod_t = torch.sqrt(alpha_t + 1e-8)
        sqrt_one_minus_alpha_t = torch.sqrt(1 - alpha_t + 1e-8)
        
        # x0_pred from denoiser is the predicted clean latent
        pred_x0 = x0_pred
        
        # Direction pointing to x_t
        pred_eps = (z_t - sqrt_alphas_cumprod_t * pred_x0) / (sqrt_one_minus_alpha_t + 1e-8)
        
        # x_{t-1} = sqrt(alpha_{t-1}) * x0 + sqrt(1 - alpha_{t-1}) * eps
        sqrt_alpha_prev = torch.sqrt(alpha_t_prev + 1e-8)
        sqrt_one_minus_alpha_prev = torch.sqrt(1 - alpha_t_prev + 1e-8)
        
        model_mean = sqrt_alpha_prev * pred_x0 + sqrt_one_minus_alpha_prev * pred_eps
        
        variance = self._extract(self.betas, t, z_t.shape)

        return model_mean, variance, x0_pred

    @torch.no_grad()
    def p_sample(
        self,
        z_t: torch.Tensor,
        t: int,
        cond: torch.Tensor = None,
    ) -> torch.Tensor:
        """
        Single reverse diffusion step (DDPM sampling).

        Args:
            z_t: Current noisy latent
            t: Current timestep
            cond: Optional conditioning

        Returns:
            z_{t-1}
        """
        t_tensor = torch.full((z_t.shape[0],), t, device=z_t.device, dtype=torch.long)
        mean, variance, _ = self.p_mean_variance(z_t, t_tensor, cond)

        if t == 0:
            return mean

        noise = torch.randn_like(z_t)
        sample = mean + torch.sqrt(variance + 1e-8) * noise
        # Clamp to prevent extreme values
        sample = torch.clamp(sample, min=-10, max=10)
        return sample

    @torch.no_grad()
    def p_sample_loop(
        self,
        shape: tuple,
        cond: torch.Tensor = None,
        ddim_steps: int = 50,
    ) -> torch.Tensor:
        """
        Full reverse diffusion sampling.

        Args:
            shape: (batch_size, latent_dim)
            cond: Optional conditioning
            ddim_steps: Number of DDIM steps (if using DDIM)

        Returns:
            Clean latent z0 (batch, latent_dim)
        """
        device = next(self.parameters()).device
        z_t = torch.randn(shape, device=device)

        # Use deterministic DDPM sampling for simplicity
        for t in reversed(range(self.diffusion_steps)):
            t_tensor = torch.full((shape[0],), t, device=device, dtype=torch.long)
            
            # Predict x0
            x0_pred = self.denoiser(z_t, t_tensor, cond)
            
            # Get parameters for this timestep
            alpha_t = self.alphas_cumprod[t].to(device)
            beta_t = self.betas[t].to(device)
            
            if t > 0:
                # Sample from posterior
                pred_mean = x0_pred * torch.sqrt(alpha_t)
                pred_mean = pred_mean + torch.sqrt(1 - alpha_t) * ((z_t - torch.sqrt(alpha_t) * x0_pred) / torch.sqrt(1 - alpha_t + 1e-8))
                
                # Add noise
                noise = torch.randn_like(z_t)
                z_t = pred_mean + torch.sqrt(beta_t + 1e-8) * noise
            else:
                z_t = x0_pred
            
            # Clamp to prevent extreme values
            z_t = torch.clamp(z_t, min=-10, max=10)

        return z_t

    def encode(self, x: torch.Tensor) -> torch.Tensor:
        """Encode input to latent space."""
        z = self.encoder(x)
        # Clamp to prevent NaN
        z = torch.clamp(z, min=-10, max=10)
        return z

    def decode(self, z: torch.Tensor) -> torch.Tensor:
        """Decode latent to expression space."""
        return self.decoder(z)

    def forward(
        self,
        x: torch.Tensor,
        mask: torch.Tensor = None,
        return_recon: bool = True,
        sample_diffusion: bool = False,
        cond: torch.Tensor = None,
    ) -> dict:
        """
        Full forward pass.

        Args:
            x: Input expression (batch, n_genes)
            mask: Optional binary mask for active genes (batch, n_genes)
            return_recon: Whether to return reconstruction
            sample_diffusion: Whether to run full diffusion sampling
            cond: Optional conditioning

        Returns:
            dict with keys: z, z_denoised, x_recon, losses
        """
        # Encode
        z = self.encode(x)

        # Diffusion (training mode)
        if self.training or not sample_diffusion:
            # Sample random timestep
            batch_size = x.shape[0]
            t = torch.randint(0, self.diffusion_steps, (batch_size,), device=x.device)

            # Add noise
            z_noisy = self.q_sample(z, t)

            # Denoise
            z_denoised = self.denoiser(z_noisy, t, cond)
        else:
            # Run full diffusion sampling
            z_denoised = self.p_sample_loop(
                shape=(x.shape[0], self.latent_dim),
                cond=cond,
            )

        # Compute losses
        losses = {}

        # Diffusion loss: predict z from noisy z
        if self.training:
            t = torch.randint(0, self.diffusion_steps, (batch_size,), device=x.device)
            z_noisy = self.q_sample(z, t)
            z_pred = self.denoiser(z_noisy, t, cond)
            losses['diffusion'] = F.mse_loss(z_pred, z)
        else:
            losses['diffusion'] = torch.tensor(0.0, device=x.device)

        # Reconstruction loss
        if return_recon:
            x_recon = self.decode(z_denoised)
            if mask is not None:
                # Only compute loss on active genes
                losses['recon'] = (F.mse_loss(x_recon, x, reduction='none') * mask).sum() / (mask.sum() + 1e-8)
            else:
                losses['recon'] = F.mse_loss(x_recon, x)
        else:
            x_recon = None

        return {
            'z': z,
            'z_denoised': z_denoised,
            'x_recon': x_recon,
            'losses': losses,
        }

    def get_clustering_embedding(
        self,
        x: torch.Tensor,
        mask: torch.Tensor = None,
        cond: torch.Tensor = None,
    ) -> torch.Tensor:
        """
        Get denoised embedding for clustering.

        This is the main entry point for clustering after training.
        """
        self.eval()
        with torch.no_grad():
            return self.forward(x, mask=mask, sample_diffusion=True, cond=cond)['z_denoised']

    def get_direct_embedding(self, x: torch.Tensor) -> torch.Tensor:
        """
        Get direct encoder embedding without diffusion sampling.
        Faster and suitable for clustering.

        Args:
            x: Input expression (batch, n_genes)

        Returns:
            Latent embedding (batch, latent_dim)
        """
        self.eval()
        with torch.no_grad():
            z = self.encode(x)
            return torch.clamp(z, min=-10, max=10)
