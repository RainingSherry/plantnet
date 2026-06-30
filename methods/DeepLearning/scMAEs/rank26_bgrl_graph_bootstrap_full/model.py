from __future__ import annotations

import copy

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


class DenseSAGELayer(nn.Module):
    """Dependency-free GraphSAGE-style layer over a dense mini-batch adjacency."""

    def __init__(self, in_features: int, out_features: int, dropout: float) -> None:
        super().__init__()
        if in_features <= 0 or out_features <= 0:
            raise ValueError("in_features and out_features must be positive")
        self.self_linear = nn.Linear(in_features, out_features, bias=False)
        self.neighbor_linear = nn.Linear(in_features, out_features, bias=False)
        self.norm = nn.LayerNorm(out_features)
        self.activation = nn.PReLU(1)
        self.dropout = nn.Dropout(dropout)

    def reset_parameters(self) -> None:
        self.self_linear.reset_parameters()
        self.neighbor_linear.reset_parameters()
        self.norm.reset_parameters()
        self.activation.weight.data.fill_(0.25)

    def forward(self, x: torch.Tensor, adjacency: torch.Tensor) -> torch.Tensor:
        if x.ndim != 2:
            raise ValueError(f"x must be [cells, features], got {tuple(x.shape)}")
        if adjacency.shape != (x.shape[0], x.shape[0]):
            raise ValueError(f"adjacency must be [cells, cells], got {tuple(adjacency.shape)}")
        norm_adj = normalize_adjacency(adjacency)
        neighbor = norm_adj @ x
        out = self.self_linear(x) + self.neighbor_linear(neighbor)
        return self.dropout(self.activation(self.norm(out)))


class BGRLDenseGraphEncoder(nn.Module):
    def __init__(self, num_genes: int, hidden_size: int, latent_size: int, dropout: float) -> None:
        super().__init__()
        if num_genes <= 0 or hidden_size <= 0 or latent_size <= 0:
            raise ValueError("num_genes, hidden_size, and latent_size must be positive")
        self.num_genes = int(num_genes)
        self.hidden_size = int(hidden_size)
        self.latent_size = int(latent_size)
        self.input_norm = nn.LayerNorm(num_genes)
        self.gnn1 = DenseSAGELayer(num_genes, hidden_size, dropout)
        self.gnn2 = DenseSAGELayer(hidden_size, hidden_size, dropout)
        self.gnn3 = DenseSAGELayer(hidden_size, latent_size, dropout)
        self.skip1 = nn.Linear(num_genes, hidden_size, bias=False)
        self.skip2 = nn.Linear(num_genes, hidden_size, bias=False)

    def reset_parameters(self) -> None:
        self.input_norm.reset_parameters()
        self.gnn1.reset_parameters()
        self.gnn2.reset_parameters()
        self.gnn3.reset_parameters()
        self.skip1.reset_parameters()
        self.skip2.reset_parameters()

    def forward(self, x: torch.Tensor, adjacency: torch.Tensor) -> torch.Tensor:
        if x.ndim != 2 or x.shape[1] != self.num_genes:
            raise ValueError(f"x must be [cells, {self.num_genes}], got {tuple(x.shape)}")
        x_norm = self.input_norm(x)
        h1 = self.gnn1(x_norm, adjacency)
        h2 = self.gnn2(h1 + self.skip1(x_norm), adjacency)
        return self.gnn3(h1 + h2 + self.skip2(x_norm), adjacency)


class MLPPredictor(nn.Module):
    def __init__(self, latent_size: int, predictor_hidden: int) -> None:
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(latent_size, predictor_hidden),
            nn.PReLU(1),
            nn.Linear(predictor_hidden, latent_size),
        )

    def reset_parameters(self) -> None:
        for module in self.modules():
            if isinstance(module, nn.Linear):
                module.reset_parameters()
            elif isinstance(module, nn.PReLU):
                module.weight.data.fill_(0.25)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.net(x)


class BGRLGraphBootstrapScMAE(nn.Module):
    """BGRL-style online/target graph bootstrap model for scRNA mini-batch graphs."""

    def __init__(
        self,
        num_genes: int,
        hidden_size: int = 128,
        latent_size: int = 64,
        predictor_hidden: int = 256,
        dropout: float = 0.1,
    ) -> None:
        super().__init__()
        self.num_genes = int(num_genes)
        self.latent_size = int(latent_size)
        self.online_encoder = BGRLDenseGraphEncoder(num_genes, hidden_size, latent_size, dropout)
        self.target_encoder = copy.deepcopy(self.online_encoder)
        self.target_encoder.reset_parameters()
        self.target_encoder.requires_grad_(False)
        self.predictor = MLPPredictor(latent_size, predictor_hidden)
        self.reconstruction_head = nn.Sequential(
            nn.Linear(latent_size, hidden_size),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_size, num_genes),
        )
        self.mask_head = nn.Linear(latent_size, num_genes)

    def trainable_parameters(self) -> list[nn.Parameter]:
        return list(self.online_encoder.parameters()) + list(self.predictor.parameters()) + list(self.reconstruction_head.parameters()) + list(self.mask_head.parameters())

    @torch.no_grad()
    def update_target_network(self, momentum: float) -> None:
        if not 0.0 <= float(momentum) <= 1.0:
            raise ValueError("momentum must be in [0, 1]")
        for online_param, target_param in zip(self.online_encoder.parameters(), self.target_encoder.parameters()):
            target_param.data.mul_(float(momentum)).add_(online_param.data, alpha=1.0 - float(momentum))

    def encode_online(self, x: torch.Tensor, adjacency: torch.Tensor) -> torch.Tensor:
        return self.online_encoder(x, adjacency)

    @torch.no_grad()
    def encode_target(self, x: torch.Tensor, adjacency: torch.Tensor) -> torch.Tensor:
        self.target_encoder.eval()
        return self.target_encoder(x, adjacency).detach()

    def forward(
        self,
        view1: torch.Tensor,
        adjacency1: torch.Tensor,
        view2: torch.Tensor,
        adjacency2: torch.Tensor,
    ) -> dict[str, torch.Tensor]:
        online_z1 = self.encode_online(view1, adjacency1)
        online_z2 = self.encode_online(view2, adjacency2)
        pred1 = self.predictor(online_z1)
        pred2 = self.predictor(online_z2)
        with torch.no_grad():
            target_z1 = self.encode_target(view1, adjacency1)
            target_z2 = self.encode_target(view2, adjacency2)
        reconstruction1 = self.reconstruction_head(online_z1)
        reconstruction2 = self.reconstruction_head(online_z2)
        mask_logits1 = self.mask_head(online_z1)
        mask_logits2 = self.mask_head(online_z2)
        return {
            "embedding": 0.5 * (online_z1 + online_z2),
            "online_z1": online_z1,
            "online_z2": online_z2,
            "prediction1": pred1,
            "prediction2": pred2,
            "target_z1": target_z1,
            "target_z2": target_z2,
            "reconstruction1": reconstruction1,
            "reconstruction2": reconstruction2,
            "mask_logits1": mask_logits1,
            "mask_logits2": mask_logits2,
        }

    @torch.no_grad()
    def feature(self, x: torch.Tensor, adjacency: torch.Tensor) -> torch.Tensor:
        return self.online_encoder(x, adjacency)
