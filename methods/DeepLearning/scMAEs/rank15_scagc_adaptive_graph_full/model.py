from __future__ import annotations

import math

import torch
import torch.nn as nn
import torch.nn.functional as F


def normalize_adjacency(adjacency: torch.Tensor) -> torch.Tensor:
    if adjacency.ndim != 2 or adjacency.shape[0] != adjacency.shape[1]:
        raise ValueError(f"adjacency must be [cells, cells], got {tuple(adjacency.shape)}")
    eye = torch.eye(adjacency.shape[0], device=adjacency.device, dtype=adjacency.dtype)
    a_hat = adjacency + eye
    degree = a_hat.sum(dim=1).clamp_min(1e-6)
    d_inv_sqrt = torch.rsqrt(degree)
    return d_inv_sqrt[:, None] * a_hat * d_inv_sqrt[None, :]


class TAGCNLayer(nn.Module):
    """Topology adaptive graph convolution with 0..K hop kernels."""

    def __init__(self, in_features: int, out_features: int, hop_order: int, dropout: float) -> None:
        super().__init__()
        if in_features <= 0 or out_features <= 0 or hop_order < 0:
            raise ValueError("in_features/out_features must be positive and hop_order must be non-negative")
        self.hop_order = int(hop_order)
        self.kernels = nn.ModuleList([nn.Linear(in_features, out_features, bias=False) for _ in range(self.hop_order + 1)])
        self.bias = nn.Parameter(torch.zeros(out_features))
        self.norm = nn.LayerNorm(out_features)
        self.dropout = nn.Dropout(dropout)

    def forward(self, x: torch.Tensor, adjacency: torch.Tensor) -> torch.Tensor:
        if x.ndim != 2:
            raise ValueError(f"x must be [cells, features], got {tuple(x.shape)}")
        norm_adj = normalize_adjacency(adjacency)
        aggregated = x
        out = self.kernels[0](aggregated)
        for order in range(1, self.hop_order + 1):
            aggregated = norm_adj @ aggregated
            out = out + self.kernels[order](aggregated)
        out = self.norm(out + self.bias)
        return self.dropout(F.relu(out))


class AdaptiveGraphSampler(nn.Module):
    """RBF similarity plus Gumbel-TopK straight-through adaptive graph sampler."""

    def __init__(self, top_k: int, sigma: float, temperature: float) -> None:
        super().__init__()
        if top_k <= 0 or sigma <= 0.0 or temperature <= 0.0:
            raise ValueError("top_k, sigma, and temperature must be positive")
        self.top_k = int(top_k)
        self.sigma = float(sigma)
        self.temperature = float(temperature)

    def rbf_similarity(self, z: torch.Tensor) -> torch.Tensor:
        if z.ndim != 2:
            raise ValueError(f"z must be [cells, hidden], got {tuple(z.shape)}")
        dist_sq = torch.cdist(z, z).square()
        sim = torch.exp(-dist_sq / (2.0 * self.sigma * self.sigma))
        diag = torch.eye(sim.shape[0], device=sim.device, dtype=torch.bool)
        return sim.masked_fill(diag, -1e9)

    def forward(self, z: torch.Tensor, sample_gumbel: bool) -> torch.Tensor:
        cells = z.shape[0]
        if cells <= 1:
            return z.new_zeros(cells, cells)
        k = min(self.top_k, cells - 1)
        logits = self.rbf_similarity(z)
        if sample_gumbel:
            uniform = torch.rand_like(logits).clamp_(1e-6, 1.0 - 1e-6)
            logits = logits - torch.log(-torch.log(uniform))
        soft = F.softmax(logits / self.temperature, dim=1)
        topk = torch.topk(logits, k=k, dim=1).indices
        hard = torch.zeros_like(soft).scatter_(1, topk, 1.0)
        directed = (hard - soft).detach() + soft
        symmetric = torch.maximum(directed, directed.T)
        eye = torch.eye(cells, device=z.device, dtype=torch.bool)
        return symmetric.masked_fill(eye, 0.0)


class ScAGCAdaptiveGraphScMAE(nn.Module):
    """scAGC-style adaptive graph autoencoder with scMAE reconstruction branch."""

    def __init__(
        self,
        num_genes: int,
        hidden_size: int = 128,
        latent_size: int = 64,
        hop_order: int = 2,
        graph_top_k: int = 15,
        graph_sigma: float = 1.0,
        graph_temperature: float = 0.7,
        dropout: float = 0.1,
    ) -> None:
        super().__init__()
        if num_genes <= 0 or hidden_size <= 0 or latent_size <= 0:
            raise ValueError("num_genes, hidden_size, and latent_size must be positive")
        self.num_genes = int(num_genes)
        self.hidden_size = int(hidden_size)
        self.latent_size = int(latent_size)
        self.input_projection = nn.Sequential(
            nn.Linear(num_genes, hidden_size),
            nn.LayerNorm(hidden_size),
            nn.ReLU(),
            nn.Dropout(dropout),
        )
        self.tagcn1 = TAGCNLayer(hidden_size, hidden_size, hop_order, dropout)
        self.tagcn2 = TAGCNLayer(hidden_size, latent_size, hop_order, dropout)
        self.graph_sampler = AdaptiveGraphSampler(graph_top_k, graph_sigma, graph_temperature)
        self.adjacency_decoder_scale = 1.0 / math.sqrt(float(latent_size))
        self.zinb_mu = nn.Linear(latent_size, num_genes)
        self.zinb_theta = nn.Linear(latent_size, num_genes)
        self.zinb_pi = nn.Linear(latent_size, num_genes)
        self.mask_predictor = nn.Linear(latent_size, num_genes)
        self.scaled_decoder = nn.Sequential(
            nn.Linear(latent_size + num_genes, hidden_size),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_size, num_genes),
        )

    def encode_with_graph(self, x: torch.Tensor, adjacency: torch.Tensor) -> torch.Tensor:
        if x.ndim != 2 or x.shape[1] != self.num_genes:
            raise ValueError(f"x must be [cells, {self.num_genes}], got {tuple(x.shape)}")
        if adjacency.shape != (x.shape[0], x.shape[0]):
            raise ValueError(f"adjacency must be [cells, cells], got {tuple(adjacency.shape)} for {x.shape[0]} cells")
        h = self.input_projection(x)
        h = self.tagcn1(h, adjacency)
        return self.tagcn2(h, adjacency)

    def decode(self, z: torch.Tensor) -> dict[str, torch.Tensor]:
        if z.ndim != 2 or z.shape[1] != self.latent_size:
            raise ValueError(f"z must be [cells, {self.latent_size}], got {tuple(z.shape)}")
        mask_logits = self.mask_predictor(z)
        return {
            "adjacency_reconstruction": torch.sigmoid((z @ z.T) * self.adjacency_decoder_scale),
            "zinb_mu": F.softplus(self.zinb_mu(z)).clamp_min(1e-5),
            "zinb_theta": F.softplus(self.zinb_theta(z)).clamp_min(1e-5),
            "zinb_pi_logits": self.zinb_pi(z),
            "mask_logits": mask_logits,
            "reconstruction": self.scaled_decoder(torch.cat([z, mask_logits], dim=1)),
        }

    def forward(self, x: torch.Tensor, initial_adjacency: torch.Tensor) -> dict[str, torch.Tensor]:
        z_previous = self.encode_with_graph(x, initial_adjacency)
        adaptive_adjacency = self.graph_sampler(z_previous, sample_gumbel=self.training)
        z_current = self.encode_with_graph(x, adaptive_adjacency)
        decoded = self.decode(z_current)
        decoded.update(
            {
                "embedding": z_current,
                "previous_embedding": z_previous,
                "adaptive_adjacency": adaptive_adjacency,
                "initial_adjacency": initial_adjacency,
            }
        )
        return decoded

    @torch.no_grad()
    def feature(self, x: torch.Tensor, initial_adjacency: torch.Tensor) -> torch.Tensor:
        z_previous = self.encode_with_graph(x, initial_adjacency)
        adaptive_adjacency = self.graph_sampler(z_previous, sample_gumbel=False)
        return self.encode_with_graph(x, adaptive_adjacency)
