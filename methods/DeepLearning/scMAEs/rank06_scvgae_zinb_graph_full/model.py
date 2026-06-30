from __future__ import annotations

import torch
import torch.nn as nn
import torch.nn.functional as F


def build_knn_adjacency(x: torch.Tensor, k: int, edge_dropout: float = 0.0) -> torch.Tensor:
    if x.ndim != 2:
        raise ValueError(f"x must be [batch, genes], got {tuple(x.shape)}")
    n_cells = x.shape[0]
    if n_cells == 0:
        raise ValueError("empty batch")
    k_eff = min(max(1, int(k)), max(1, n_cells - 1))
    with torch.no_grad():
        x_norm = F.normalize(x, dim=1)
        sim = x_norm @ x_norm.t()
        sim.fill_diagonal_(-float("inf"))
        nn_idx = sim.topk(k_eff, dim=1).indices
        adj = torch.zeros((n_cells, n_cells), dtype=x.dtype, device=x.device)
        rows = torch.arange(n_cells, device=x.device).view(-1, 1).expand_as(nn_idx)
        adj[rows, nn_idx] = 1.0
        adj = torch.maximum(adj, adj.t())
        if edge_dropout > 0.0:
            keep = (torch.rand_like(adj) >= float(edge_dropout)).float()
            keep = torch.triu(keep, diagonal=1)
            keep = keep + keep.t()
            adj = adj * keep
        adj.fill_diagonal_(1.0)
        degree = adj.sum(dim=1).clamp_min(1.0)
        adj = degree.rsqrt().view(-1, 1) * adj * degree.rsqrt().view(1, -1)
    return adj


class GraphNorm(nn.Module):
    def __init__(self, hidden_size: int) -> None:
        super().__init__()
        self.weight = nn.Parameter(torch.ones(hidden_size))
        self.bias = nn.Parameter(torch.zeros(hidden_size))

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        mean = x.mean(dim=0, keepdim=True)
        var = (x - mean).square().mean(dim=0, keepdim=True)
        return (x - mean) / torch.sqrt(var + 1e-5) * self.weight + self.bias


class GraphConvolution(nn.Module):
    def __init__(self, in_features: int, out_features: int, activation: bool = True) -> None:
        super().__init__()
        self.linear = nn.Linear(in_features, out_features)
        self.norm = GraphNorm(out_features)
        self.activation = activation

    def forward(self, x: torch.Tensor, adj: torch.Tensor) -> torch.Tensor:
        if adj.ndim != 2 or adj.shape[0] != adj.shape[1] or adj.shape[0] != x.shape[0]:
            raise ValueError(f"adj must be [batch, batch], got {tuple(adj.shape)} for x {tuple(x.shape)}")
        h = adj @ self.linear(x)
        h = self.norm(h)
        return F.relu(h) if self.activation else h


class ScVGAEZINBScMAE(nn.Module):
    """ZINB-based variational graph autoencoder adapted to scMAE.

    The model follows the scVGAE method structure from the paper: a GCN encoder
    over a cell-cell similarity graph, variational latent variables, graph
    branches that predict ZINB mean/dropout/dispersion, and a reconstruction
    decoder. The masked reconstruction branch keeps scMAE's corruption objective.
    """

    def __init__(
        self,
        num_genes: int,
        hidden_size: int = 128,
        latent_size: int = 64,
        dropout: float = 0.2,
        knn_k: int = 15,
        edge_dropout: float = 0.05,
    ) -> None:
        super().__init__()
        if num_genes <= 0:
            raise ValueError("num_genes must be positive")
        self.num_genes = int(num_genes)
        self.hidden_size = int(hidden_size)
        self.latent_size = int(latent_size)
        self.knn_k = int(knn_k)
        self.edge_dropout = float(edge_dropout)
        self.dropout = nn.Dropout(dropout)

        self.encoder = GraphConvolution(num_genes, hidden_size, activation=True)
        self.z_mean = GraphConvolution(hidden_size, latent_size, activation=False)
        self.z_logvar = GraphConvolution(hidden_size, latent_size, activation=False)

        self.mean_head = GraphConvolution(latent_size, num_genes, activation=False)
        self.dispersion_head = GraphConvolution(latent_size, num_genes, activation=False)
        self.dropout_head = GraphConvolution(latent_size, num_genes, activation=False)

        self.decoder = nn.Sequential(
            nn.Linear(latent_size, hidden_size),
            nn.LayerNorm(hidden_size),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_size, num_genes),
        )
        self.mask_predictor = nn.Linear(latent_size, num_genes)
        self.embedding_projector = nn.Sequential(
            nn.LayerNorm(latent_size),
            nn.Linear(latent_size, hidden_size),
            nn.GELU(),
            nn.Linear(hidden_size, hidden_size),
        )

    def encode(self, x: torch.Tensor, adj: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
        if x.ndim != 2 or x.shape[1] != self.num_genes:
            raise ValueError(f"x must be [batch, {self.num_genes}], got {tuple(x.shape)}")
        hidden = self.dropout(self.encoder(x, adj))
        return self.z_mean(hidden, adj), self.z_logvar(hidden, adj).clamp(-8.0, 8.0)

    def reparameterize(self, mean: torch.Tensor, logvar: torch.Tensor) -> torch.Tensor:
        if self.training:
            return mean + torch.randn_like(mean) * torch.exp(0.5 * logvar)
        return mean

    def forward(self, x: torch.Tensor, adj: torch.Tensor | None = None) -> dict[str, torch.Tensor]:
        if adj is None:
            adj = build_knn_adjacency(x, self.knn_k, self.edge_dropout if self.training else 0.0)
        z_mean, z_logvar = self.encode(x, adj)
        z = self.reparameterize(z_mean, z_logvar)
        mu = F.softplus(self.mean_head(z, adj)).clamp_min(1e-4)
        theta = F.softplus(self.dispersion_head(z, adj)).clamp_min(1e-4)
        pi_logits = self.dropout_head(z, adj)
        return {
            "embedding": self.embedding_projector(z_mean),
            "z": z,
            "z_mean": z_mean,
            "z_logvar": z_logvar,
            "zinb_mu": mu,
            "zinb_theta": theta,
            "zinb_pi_logits": pi_logits,
            "reconstruction": self.decoder(z),
            "mask_logits": self.mask_predictor(z),
            "adj": adj,
        }

    @torch.no_grad()
    def feature(self, x: torch.Tensor, batch_size: int = 512) -> torch.Tensor:
        outputs = []
        for start in range(0, x.shape[0], int(batch_size)):
            xb = x[start:start + int(batch_size)]
            adj = build_knn_adjacency(xb, self.knn_k, 0.0)
            z_mean, _ = self.encode(xb, adj)
            outputs.append(self.embedding_projector(z_mean))
        return torch.cat(outputs, dim=0)
