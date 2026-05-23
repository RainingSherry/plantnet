# -*- coding: utf-8 -*-
"""
Encoder and Decoder modules for scRNA-seq data.

Provides flexible encoder/decoder architectures for the maskdiffusion pipeline.
"""

import torch
import torch.nn as nn
import torch.nn.functional as F


class SparseEncoder(nn.Module):
    """
    Sparse-aware encoder for scRNA-seq data.

    Key features:
    - LayerNorm instead of BatchNorm (better for small batch sizes)
    - Mish activation (smoother than ReLU)
    - Optional support-aware input weighting
    """

    def __init__(
        self,
        input_dim: int,
        latent_dim: int,
        hidden_dims: list = None,
        dropout: float = 0.1,
        use_batch_norm: bool = False,
    ):
        super().__init__()
        if hidden_dims is None:
            hidden_dims = [512, 256, 128]

        self.input_dim = input_dim
        self.latent_dim = latent_dim

        # Build encoder layers
        layers = []
        in_dim = input_dim
        for i, hidden_dim in enumerate(hidden_dims):
            layers.append(nn.Dropout(p=dropout))
            layers.append(nn.Linear(in_dim, hidden_dim))
            if use_batch_norm:
                layers.append(nn.BatchNorm1d(hidden_dim))
            else:
                layers.append(nn.LayerNorm(hidden_dim))
            layers.append(nn.Mish(inplace=True))
            in_dim = hidden_dim

        # Latent projection
        layers.append(nn.Linear(in_dim, latent_dim))

        self.encoder = nn.Sequential(*layers)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        """
        Args:
            x: Input expression (batch, n_genes)

        Returns:
            Latent representation (batch, latent_dim)
        """
        return self.encoder(x)


class SupportAwareEncoder(nn.Module):
    """
    Encoder that takes both expression and support (mask) as input.

    Architecture:
        X + Support_M → MLP → z
    """

    def __init__(
        self,
        input_dim: int,
        latent_dim: int,
        hidden_dims: list = None,
        dropout: float = 0.1,
    ):
        super().__init__()
        if hidden_dims is None:
            hidden_dims = [512, 256, 128]

        self.input_dim = input_dim
        self.latent_dim = latent_dim

        # Expression encoder
        self.expr_encoder = nn.Sequential(
            nn.Dropout(p=dropout),
            nn.Linear(input_dim, hidden_dims[0]),
            nn.LayerNorm(hidden_dims[0]),
            nn.Mish(inplace=True),
        )

        # Support encoder (predicts which genes are active)
        self.support_encoder = nn.Sequential(
            nn.Dropout(p=dropout),
            nn.Linear(input_dim, hidden_dims[0] // 2),
            nn.LayerNorm(hidden_dims[0] // 2),
            nn.Mish(inplace=True),
        )

        # Fusion layers
        fusion_dims = [hidden_dims[0] + hidden_dims[0] // 2] + hidden_dims[1:]
        in_dim = fusion_dims[0]
        fusion_layers = []
        for hidden_dim in fusion_dims[1:]:
            fusion_layers.extend([
                nn.Linear(in_dim, hidden_dim),
                nn.LayerNorm(hidden_dim),
                nn.Mish(inplace=True),
                nn.Dropout(p=dropout),
            ])
            in_dim = hidden_dim

        self.fusion = nn.Sequential(*fusion_layers)
        self.latent_proj = nn.Linear(in_dim, latent_dim)

    def forward(
        self,
        x: torch.Tensor,
        support: torch.Tensor = None,
    ) -> torch.Tensor:
        """
        Args:
            x: Input expression (batch, n_genes)
            support: Support mask (batch, n_genes), 1=active, 0=zero
                     If None, computed as (x > 0).float()

        Returns:
            Latent representation (batch, latent_dim)
        """
        if support is None:
            support = (x > 0).float()

        # Encode expression
        h_expr = self.expr_encoder(x)

        # Encode support
        h_support = self.support_encoder(support)

        # Fuse and project
        h = torch.cat([h_expr, h_support], dim=-1)
        h = self.fusion(h)
        z = self.latent_proj(h)

        return z


class VariationalEncoder(nn.Module):
    """
    Variational Autoencoder-style encoder with reparameterization trick.

    Outputs mean and log variance for latent distribution.
    """

    def __init__(
        self,
        input_dim: int,
        latent_dim: int,
        hidden_dims: list = None,
        dropout: float = 0.1,
    ):
        super().__init__()
        if hidden_dims is None:
            hidden_dims = [512, 256]

        self.input_dim = input_dim
        self.latent_dim = latent_dim

        # Build encoder
        layers = []
        in_dim = input_dim
        for hidden_dim in hidden_dims:
            layers.extend([
                nn.Dropout(p=dropout),
                nn.Linear(in_dim, hidden_dim),
                nn.LayerNorm(hidden_dim),
                nn.Mish(inplace=True),
            ])
            in_dim = hidden_dim

        self.encoder = nn.Sequential(*layers)

        # Latent mean and variance
        self.fc_mu = nn.Linear(in_dim, latent_dim)
        self.fc_logvar = nn.Linear(in_dim, latent_dim)

    def reparameterize(self, mu: torch.Tensor, logvar: torch.Tensor) -> torch.Tensor:
        """Reparameterization trick for VAE."""
        if self.training:
            std = torch.exp(0.5 * logvar)
            eps = torch.randn_like(std)
            return mu + eps * std
        return mu

    def forward(self, x: torch.Tensor) -> tuple:
        """
        Args:
            x: Input expression (batch, n_genes)

        Returns:
            (z, mu, logvar)
        """
        h = self.encoder(x)
        mu = self.fc_mu(h)
        logvar = self.fc_logvar(h)
        z = self.reparameterize(mu, logvar)
        return z, mu, logvar


class Decoder(nn.Module):
    """
    Decoder for reconstructing expression from latent.

    Supports masked reconstruction (only on active genes).
    """

    def __init__(
        self,
        latent_dim: int,
        output_dim: int,
        hidden_dims: list = None,
        dropout: float = 0.1,
    ):
        super().__init__()
        if hidden_dims is None:
            hidden_dims = [128, 256, 512]

        self.latent_dim = latent_dim
        self.output_dim = output_dim

        # Build decoder layers
        layers = []
        in_dim = latent_dim
        for hidden_dim in hidden_dims:
            layers.extend([
                nn.Dropout(p=dropout),
                nn.Linear(in_dim, hidden_dim),
                nn.LayerNorm(hidden_dim),
                nn.Mish(inplace=True),
            ])
            in_dim = hidden_dim

        layers.append(nn.Linear(in_dim, output_dim))

        self.decoder = nn.Sequential(*layers)

    def forward(
        self,
        z: torch.Tensor,
        mask: torch.Tensor = None,
    ) -> torch.Tensor:
        """
        Args:
            z: Latent representation (batch, latent_dim)
            mask: Optional mask for active genes (batch, n_genes)

        Returns:
            Reconstructed expression (batch, n_genes)
        """
        x_recon = self.decoder(z)

        if mask is not None:
            # Zero out non-active genes
            x_recon = x_recon * mask

        return x_recon


class GeneWiseDecoder(nn.Module):
    """
    Gene-wise decoder that predicts each gene independently.

    Used for gene-level tasks like marker gene prediction.
    """

    def __init__(
        self,
        latent_dim: int,
        num_genes: int,
        hidden_dim: int = 128,
        dropout: float = 0.1,
    ):
        super().__init__()
        self.num_genes = num_genes

        # Per-gene prediction heads
        self.gene_predictors = nn.ModuleList([
            nn.Sequential(
                nn.Linear(latent_dim, hidden_dim),
                nn.LayerNorm(hidden_dim),
                nn.Mish(inplace=True),
                nn.Dropout(dropout),
                nn.Linear(hidden_dim, 1),
            )
            for _ in range(num_genes)
        ])

    def forward(self, z: torch.Tensor) -> torch.Tensor:
        """
        Args:
            z: Latent representation (batch, latent_dim)

        Returns:
            Gene expression predictions (batch, num_genes)
        """
        predictions = []
        for predictor in self.gene_predictors:
            predictions.append(predictor(z))
        return torch.cat(predictions, dim=-1)
