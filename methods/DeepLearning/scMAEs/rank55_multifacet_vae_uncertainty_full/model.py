from __future__ import annotations

import torch
from torch import nn
import torch.nn.functional as F


class MultiFacetVAEscMAE(nn.Module):
    """scMAE with MFCVAE-style multiple latent facets and MoG priors."""

    def __init__(
        self,
        input_dim: int,
        hidden_size: int = 128,
        decoder_hidden: int = 128,
        n_facets: int = 2,
        facet_dim: int = 32,
        n_components: int = 8,
        dropout: float = 0.0,
    ):
        super().__init__()
        self.input_dim = int(input_dim)
        self.hidden_size = int(hidden_size)
        self.n_facets = int(n_facets)
        self.facet_dim = int(facet_dim)
        self.n_components = int(n_components)
        self.encoder = nn.Sequential(
            nn.Dropout(dropout),
            nn.Linear(input_dim, 256),
            nn.LayerNorm(256),
            nn.GELU(),
            nn.Linear(256, hidden_size),
            nn.LayerNorm(hidden_size),
            nn.GELU(),
        )
        self.mu_heads = nn.ModuleList([nn.Linear(hidden_size, facet_dim) for _ in range(n_facets)])
        self.logvar_heads = nn.ModuleList([nn.Linear(hidden_size, facet_dim) for _ in range(n_facets)])
        total_z = n_facets * facet_dim
        self.embedding_head = nn.Sequential(nn.Linear(total_z, hidden_size), nn.LayerNorm(hidden_size), nn.GELU())
        self.mask_predictor = nn.Linear(hidden_size, input_dim)
        self.decoder = nn.Sequential(nn.Linear(hidden_size + input_dim, decoder_hidden), nn.GELU(), nn.Linear(decoder_hidden, input_dim))
        self.nb_mean_decoder = nn.Sequential(nn.Linear(hidden_size, decoder_hidden), nn.GELU(), nn.Linear(decoder_hidden, input_dim))
        self.nb_log_theta = nn.Parameter(torch.zeros(input_dim))
        self.mog_logits = nn.Parameter(torch.zeros(n_facets, n_components))
        self.mog_mu = nn.Parameter(torch.randn(n_facets, n_components, facet_dim) * 0.02)
        self.mog_logvar = nn.Parameter(torch.zeros(n_facets, n_components, facet_dim))

    def encode_facets(self, x: torch.Tensor) -> tuple[list[torch.Tensor], list[torch.Tensor]]:
        h = self.encoder(x)
        mus = [head(h) for head in self.mu_heads]
        logvars = [head(h).clamp(-8.0, 6.0) for head in self.logvar_heads]
        return mus, logvars

    def reparameterize(self, mus: list[torch.Tensor], logvars: list[torch.Tensor]) -> list[torch.Tensor]:
        if not self.training:
            return mus
        return [mu + torch.randn_like(mu) * torch.exp(0.5 * logvar) for mu, logvar in zip(mus, logvars)]

    def mixture_responsibilities(self, mus: list[torch.Tensor]) -> list[torch.Tensor]:
        resp = []
        for j, mu in enumerate(mus):
            diff = mu[:, None, :] - self.mog_mu[j][None, :, :]
            logvar = self.mog_logvar[j][None, :, :].clamp(-8.0, 6.0)
            log_prob = -0.5 * (diff.pow(2) / torch.exp(logvar) + logvar + torch.log(torch.tensor(2.0 * torch.pi, device=mu.device))).sum(dim=-1)
            log_prob = log_prob + F.log_softmax(self.mog_logits[j], dim=0).view(1, -1)
            resp.append(torch.softmax(log_prob, dim=1))
        return resp

    def forward(self, x: torch.Tensor) -> dict[str, torch.Tensor | list[torch.Tensor]]:
        mus, logvars = self.encode_facets(x)
        zs = self.reparameterize(mus, logvars)
        z_cat = torch.cat(zs, dim=1)
        emb = self.embedding_head(z_cat)
        mask_logits = self.mask_predictor(emb)
        recon = self.decoder(torch.cat([emb, mask_logits], dim=1))
        nb_mean = F.softplus(self.nb_mean_decoder(emb)) + 1e-4
        return {
            "embedding": emb,
            "reconstruction": recon,
            "mask_logits": mask_logits,
            "mus": mus,
            "logvars": logvars,
            "zs": zs,
            "responsibilities": self.mixture_responsibilities(mus),
            "nb_mean": nb_mean,
            "nb_theta": F.softplus(self.nb_log_theta) + 1e-4,
        }
