from __future__ import annotations

import torch
from torch import nn
import torch.nn.functional as F


class ClusterCentricPrototypeAdapter(nn.Module):
    """Cluster-centric prototype propagation adapted from CCSM for cell latents."""

    def __init__(self, latent_dim: int, n_prototypes: int, temperature: float = 0.25, proto_weight: float = 0.15):
        super().__init__()
        self.n_prototypes = int(n_prototypes)
        self.temperature = float(temperature)
        self.proto_weight = float(proto_weight)
        self.prototypes = nn.Parameter(torch.randn(n_prototypes, latent_dim) * 0.02)
        self.proto_refine = nn.Sequential(
            nn.LayerNorm(latent_dim),
            nn.Linear(latent_dim, latent_dim),
            nn.GELU(),
            nn.Linear(latent_dim, latent_dim),
        )
        self.cell_gate = nn.Sequential(
            nn.LayerNorm(latent_dim * 3),
            nn.Linear(latent_dim * 3, latent_dim),
            nn.GELU(),
            nn.Linear(latent_dim, 1),
            nn.Sigmoid(),
        )

    @torch.no_grad()
    def set_prototypes(self, centers: torch.Tensor) -> None:
        if centers.shape != self.prototypes.shape:
            raise ValueError(f"Expected centers {tuple(self.prototypes.shape)}, got {tuple(centers.shape)}")
        self.prototypes.copy_(centers.to(self.prototypes.device, dtype=self.prototypes.dtype))

    def propagated_prototypes(self) -> tuple[torch.Tensor, torch.Tensor]:
        proto = self.prototypes
        sim = torch.matmul(F.normalize(proto, dim=1), F.normalize(proto, dim=1).t())
        sim = sim - torch.eye(self.n_prototypes, device=sim.device, dtype=sim.dtype) * 1e4
        graph = F.softmax(sim / max(self.temperature, 1e-4), dim=1)
        context = graph @ proto
        propagated = proto + 0.5 * self.proto_refine(context)
        return propagated, graph

    def forward(self, z: torch.Tensor) -> dict:
        propagated, proto_graph = self.propagated_prototypes()
        logits = torch.matmul(F.normalize(z, dim=1), F.normalize(propagated, dim=1).t()) / max(self.temperature, 1e-4)
        assign = F.softmax(logits, dim=1)
        context = assign @ propagated
        gate = self.cell_gate(torch.cat([z, context, context - z], dim=1))
        adapted = z + self.proto_weight * gate * (context - z)
        return {
            "adapted": adapted,
            "assignment": assign,
            "proto_context": context,
            "proto_graph": proto_graph,
            "proto_gate": gate,
            "prototypes": propagated,
        }


class PrototypeGraphScMAE(nn.Module):
    """Independent scMAE with cluster-centric prototype graph adapter."""

    def __init__(
        self,
        input_dim: int,
        n_prototypes: int,
        hidden_dim: int = 512,
        latent_dim: int = 32,
        dropout: float = 0.1,
        mask_prob: float = 0.4,
        proto_weight: float = 0.15,
        proto_temperature: float = 0.25,
    ):
        super().__init__()
        self.mask_prob = float(mask_prob)
        self.encoder = nn.Sequential(
            nn.Linear(input_dim, hidden_dim),
            nn.LayerNorm(hidden_dim),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim, hidden_dim // 2),
            nn.LayerNorm(hidden_dim // 2),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim // 2, latent_dim),
        )
        self.adapter = ClusterCentricPrototypeAdapter(latent_dim, n_prototypes, proto_temperature, proto_weight)
        self.mask_predictor = nn.Sequential(nn.Linear(latent_dim, hidden_dim // 2), nn.GELU(), nn.Linear(hidden_dim // 2, input_dim))
        self.decoder = nn.Sequential(
            nn.Linear(latent_dim, hidden_dim // 2),
            nn.LayerNorm(hidden_dim // 2),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim // 2, hidden_dim),
            nn.GELU(),
            nn.Linear(hidden_dim, input_dim),
        )

    def corrupt(self, x: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
        mask = (torch.rand_like(x) < self.mask_prob).float()
        return x * (1.0 - mask), mask

    def encode_with_adapter(self, x: torch.Tensor) -> dict:
        base = self.encoder(x)
        proto = self.adapter(base)
        proto["base_latent"] = base
        proto["latent"] = proto["adapted"]
        return proto

    def forward(self, x: torch.Tensor) -> dict:
        corrupted, mask = self.corrupt(x)
        masked = self.encode_with_adapter(corrupted)
        clean = self.encode_with_adapter(x)
        latent = masked["latent"]
        return {
            "latent": latent,
            "base_latent": masked["base_latent"],
            "clean_assignment": clean["assignment"],
            "assignment": masked["assignment"],
            "proto_context": masked["proto_context"],
            "proto_graph": masked["proto_graph"],
            "proto_gate": masked["proto_gate"],
            "prototypes": masked["prototypes"],
            "reconstruction": self.decoder(latent),
            "mask_logits": self.mask_predictor(latent),
            "mask": mask,
        }

    def feature(self, x: torch.Tensor) -> torch.Tensor:
        return self.encode_with_adapter(x)["latent"]
