# -*- coding: utf-8 -*-
"""
Decoder modules for scRNA-seq reconstruction.

Provides various decoder architectures for reconstructing expression from latent space.
"""

import torch
import torch.nn as nn
import torch.nn.functional as F


class MaskedDecoder(nn.Module):
    """
    Decoder that reconstructs expression only on active genes.

    Key design:
    - Takes latent z and gene activation mask as input
    - Only outputs meaningful values for active genes
    - For inactive genes, outputs zero
    """

    def __init__(
        self,
        latent_dim: int,
        num_genes: int,
        hidden_dims: list = None,
        dropout: float = 0.1,
    ):
        super().__init__()
        if hidden_dims is None:
            hidden_dims = [128, 256, 512]

        self.latent_dim = latent_dim
        self.num_genes = num_genes

        # Build decoder
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

        layers.append(nn.Linear(in_dim, num_genes))
        self.decoder = nn.Sequential(*layers)

    def forward(
        self,
        z: torch.Tensor,
        gene_mask: torch.Tensor = None,
    ) -> torch.Tensor:
        """
        Args:
            z: Latent representation (batch, latent_dim)
            gene_mask: Binary mask for active genes (batch, num_genes)
                      If provided, output is multiplied by mask

        Returns:
            Reconstructed expression (batch, num_genes)
        """
        x_recon = self.decoder(z)

        if gene_mask is not None:
            # Zero out inactive genes
            x_recon = x_recon * gene_mask

        return x_recon


class ConditionalDecoder(nn.Module):
    """
    Decoder with conditioning on auxiliary information.

    Supports conditioning on:
    - Cluster labels
    - Cell type embeddings
    - Batch information
    """

    def __init__(
        self,
        latent_dim: int,
        num_genes: int,
        cond_dim: int = 0,
        hidden_dims: list = None,
        dropout: float = 0.1,
    ):
        super().__init__()
        if hidden_dims is None:
            hidden_dims = [128, 256, 512]

        self.latent_dim = latent_dim
        self.num_genes = num_genes
        self.cond_dim = cond_dim

        # Latent projection
        self.z_proj = nn.Sequential(
            nn.Linear(latent_dim, hidden_dims[0]),
            nn.LayerNorm(hidden_dims[0]),
            nn.Mish(inplace=True),
        )

        # Conditioning projection
        if cond_dim > 0:
            self.cond_proj = nn.Sequential(
                nn.Linear(cond_dim, hidden_dims[0]),
                nn.LayerNorm(hidden_dims[0]),
                nn.Mish(inplace=True),
            )

        # Decoder layers
        layers = []
        in_dim = hidden_dims[0] * 2 if cond_dim > 0 else hidden_dims[0]
        for hidden_dim in hidden_dims[1:]:
            layers.extend([
                nn.Dropout(p=dropout),
                nn.Linear(in_dim, hidden_dim),
                nn.LayerNorm(hidden_dim),
                nn.Mish(inplace=True),
            ])
            in_dim = hidden_dim

        layers.append(nn.Linear(in_dim, num_genes))
        self.decoder = nn.Sequential(*layers)

    def forward(
        self,
        z: torch.Tensor,
        cond: torch.Tensor = None,
        gene_mask: torch.Tensor = None,
    ) -> torch.Tensor:
        """
        Args:
            z: Latent representation (batch, latent_dim)
            cond: Conditioning vector (batch, cond_dim)
            gene_mask: Binary mask for active genes

        Returns:
            Reconstructed expression (batch, num_genes)
        """
        h = self.z_proj(z)

        if self.cond_dim > 0 and cond is not None:
            cond_emb = self.cond_proj(cond)
            h = torch.cat([h, cond_emb], dim=-1)

        x_recon = self.decoder(h)

        if gene_mask is not None:
            x_recon = x_recon * gene_mask

        return x_recon


class ZeroInflatedDecoder(nn.Module):
    """
    Decoder that models zero-inflated expression.

    Predicts:
    1. Probability of gene being expressed (Bernoulli)
    2. Expression level if expressed (Gaussian/Negative Binomial)
    """

    def __init__(
        self,
        latent_dim: int,
        num_genes: int,
        hidden_dims: list = None,
        dropout: float = 0.1,
    ):
        super().__init__()
        if hidden_dims is None:
            hidden_dims = [256, 128]

        self.latent_dim = latent_dim
        self.num_genes = num_genes

        # Shared encoder
        self.encoder = nn.Sequential(
            nn.Dropout(p=dropout),
            nn.Linear(latent_dim, hidden_dims[0]),
            nn.LayerNorm(hidden_dims[0]),
            nn.Mish(inplace=True),
        )

        # Bernoulli head (probability of expression)
        self.pi_head = nn.Sequential(
            nn.Linear(hidden_dims[0], hidden_dims[-1]),
            nn.LayerNorm(hidden_dims[-1]),
            nn.Mish(inplace=True),
            nn.Linear(hidden_dims[-1], num_genes),
        )

        # Gaussian mean head
        self.mu_head = nn.Sequential(
            nn.Linear(hidden_dims[0], hidden_dims[-1]),
            nn.LayerNorm(hidden_dims[-1]),
            nn.Mish(inplace=True),
            nn.Linear(hidden_dims[-1], num_genes),
        )

        # Log variance head for Gaussian
        self.logvar_head = nn.Sequential(
            nn.Linear(hidden_dims[0], hidden_dims[-1]),
            nn.LayerNorm(hidden_dims[-1]),
            nn.Mish(inplace=True),
            nn.Linear(hidden_dims[-1], num_genes),
        )

    def forward(
        self,
        z: torch.Tensor,
        sample: bool = False,
    ) -> dict:
        """
        Args:
            z: Latent representation (batch, latent_dim)
            sample: If True, sample from the distribution

        Returns:
            dict with keys:
                - pi: probability of expression (batch, num_genes)
                - mu: expression mean (batch, num_genes)
                - logvar: log variance (batch, num_genes)
                - expression: sampled expression (batch, num_genes)
        """
        h = self.encoder(z)

        pi = torch.sigmoid(self.pi_head(h))
        mu = self.mu_head(h)
        logvar = self.logvar_head(h)

        # Sample expression
        if sample:
            std = torch.exp(0.5 * logvar)
            eps = torch.randn_like(std)
            expression = mu + eps * std
            expression = torch.clamp(expression, min=0)
        else:
            expression = pi * mu

        return {
            'pi': pi,
            'mu': mu,
            'logvar': logvar,
            'expression': expression,
        }

    def get_loss(
        self,
        z: torch.Tensor,
        x: torch.Tensor,
        mask: torch.Tensor = None,
    ) -> dict:
        """
        Compute zero-inflated loss.

        Loss = BCE(pi, mask) + Gaussian_NLL(mu, logvar, x | mask=1)
        """
        if mask is None:
            mask = (x > 0).float()

        output = self.forward(z, sample=False)

        # BCE loss for expression probability
        bce_loss = F.binary_cross_entropy(output['pi'], mask, reduction='mean')

        # Negative log likelihood for expressed genes
        mask_expanded = mask.unsqueeze(-1) if output['mu'].dim() > 2 else mask
        if mask_expanded.dim() == 1:
            mask_expanded = mask_expanded.unsqueeze(-1)

        # Compute NLL only for expressed genes
        expressed_mask = mask_expanded > 0
        if expressed_mask.sum() > 0:
            mu_exp = output['mu'][expressed_mask]
            logvar_exp = output['logvar'][expressed_mask]
            x_exp = x[mask_expanded > 0]

            var = torch.exp(logvar_exp) + 1e-6
            nll = 0.5 * (torch.log(var) + (x_exp - mu_exp) ** 2 / var + torch.log(2 * torch.tensor(np.pi)))
            gaussian_loss = nll.mean()
        else:
            gaussian_loss = torch.tensor(0.0, device=z.device)

        total_loss = bce_loss + gaussian_loss

        return {
            'loss': total_loss,
            'bce': bce_loss,
            'gaussian': gaussian_loss,
        }
