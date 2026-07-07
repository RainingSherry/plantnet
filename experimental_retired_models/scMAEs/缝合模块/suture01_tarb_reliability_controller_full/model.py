from __future__ import annotations

import torch
from torch import nn
import torch.nn.functional as F


class LatentOperationController(nn.Module):
    """TARB-inspired operation selector rewritten for scRNA latent features.

    The original TARB chooses among image-restoration operations. Here the
    operations are small latent residual transforms. The full-gene scMAE encoder
    remains the anchor path, so setting module_weight=0 recovers the backbone.
    """

    def __init__(self, hidden_size: int, dropout: float = 0.05):
        super().__init__()
        self.hidden_size = int(hidden_size)
        self.ops = nn.ModuleList(
            [
                nn.Identity(),
                nn.Sequential(nn.LayerNorm(hidden_size), nn.Linear(hidden_size, hidden_size), nn.Mish(), nn.Dropout(dropout)),
                nn.Sequential(nn.LayerNorm(hidden_size), nn.Linear(hidden_size, hidden_size), nn.Tanh()),
                nn.Sequential(nn.LayerNorm(hidden_size), nn.Linear(hidden_size, hidden_size), nn.Sigmoid()),
                nn.Sequential(nn.LayerNorm(hidden_size), nn.Linear(hidden_size, hidden_size), nn.GELU(), nn.Linear(hidden_size, hidden_size)),
            ]
        )
        self.controller = nn.Sequential(
            nn.LayerNorm(hidden_size + 4),
            nn.Linear(hidden_size + 4, hidden_size),
            nn.Mish(),
            nn.Dropout(dropout),
            nn.Linear(hidden_size, len(self.ops)),
        )
        self.out_norm = nn.LayerNorm(hidden_size)

    def forward(self, z: torch.Tensor, reliability: torch.Tensor | None = None) -> tuple[torch.Tensor, torch.Tensor]:
        bsz = z.shape[0]
        if reliability is None:
            reliability = z.new_ones(bsz)
        reliability = reliability.clamp(0.0, 1.0).view(-1, 1)
        stats = torch.cat(
            [
                z.mean(dim=1, keepdim=True),
                z.std(dim=1, keepdim=True).clamp_min(1e-6),
                reliability,
                1.0 - reliability,
            ],
            dim=1,
        )
        weights = F.softmax(self.controller(torch.cat([z, stats], dim=1)), dim=1)
        outputs = []
        for op_id, op in enumerate(self.ops):
            val = op(z)
            if op_id == 2:
                # Variance-preserving projection: center per batch before adding
                # so it discourages collapse without changing mean too much.
                val = val - val.mean(dim=0, keepdim=True)
            elif op_id in {3, 4}:
                # Boundary-conservative/light-denoise operations are risky
                # smoothing-like terms, so reliability controls only these.
                val = val * reliability
            outputs.append(val * weights[:, op_id : op_id + 1])
        delta = torch.stack(outputs, dim=0).sum(dim=0)
        return self.out_norm(delta), weights


class SutureTARBScMAE(nn.Module):
    """Independent scMAE candidate with a TARB-style latent controller."""

    def __init__(
        self,
        num_genes: int,
        n_clusters: int,
        hidden_size: int = 128,
        dropout: float = 0.05,
        module_weight: float = 0.15,
    ):
        super().__init__()
        self.num_genes = int(num_genes)
        self.n_clusters = int(n_clusters)
        self.hidden_size = int(hidden_size)
        self.module_weight = float(module_weight)
        self.encoder = nn.Sequential(
            nn.Dropout(float(dropout)),
            nn.Linear(self.num_genes, 256),
            nn.LayerNorm(256),
            nn.Mish(),
            nn.Linear(256, hidden_size),
            nn.LayerNorm(hidden_size),
            nn.Mish(),
            nn.Linear(hidden_size, hidden_size),
        )
        self.controller = LatentOperationController(hidden_size, dropout=dropout)
        self.mask_predictor = nn.Linear(hidden_size, self.num_genes)
        self.decoder = nn.Linear(hidden_size + self.num_genes, self.num_genes)
        self.cluster_centers = nn.Parameter(torch.randn(self.n_clusters, hidden_size) * 0.02)

    def forward(self, x: torch.Tensor, reliability: torch.Tensor | None = None) -> dict[str, torch.Tensor]:
        base = self.encoder(x)
        delta, op_weights = self.controller(base, reliability)
        latent = base + self.module_weight * delta
        mask_logits = self.mask_predictor(latent)
        reconstruction = self.decoder(torch.cat([latent, mask_logits], dim=1))
        q = self.student_q(latent)
        return {
            "latent": latent,
            "base_latent": base,
            "operation_weights": op_weights,
            "mask_logits": mask_logits,
            "reconstruction": reconstruction,
            "cluster_q": q,
        }

    def student_q(self, latent: torch.Tensor) -> torch.Tensor:
        dist = torch.cdist(latent, self.cluster_centers).pow(2)
        q = 1.0 / (1.0 + dist)
        return q / q.sum(dim=1, keepdim=True).clamp_min(1e-8)

    @torch.no_grad()
    def feature(self, x: torch.Tensor) -> torch.Tensor:
        return self.forward(x, reliability=None)["latent"]

    def random_mask(self, x: torch.Tensor, mask_prob: float) -> tuple[torch.Tensor, torch.Tensor]:
        mask = (torch.rand_like(x) < float(mask_prob)).float()
        empty = mask.sum(dim=1) == 0
        if bool(empty.any()):
            cols = torch.randint(0, x.shape[1], (int(empty.sum()),), device=x.device)
            mask[empty, cols] = 1.0
        return x.masked_fill(mask.bool(), 0.0), mask

    @torch.no_grad()
    def initialize_centers(self, centers: torch.Tensor) -> None:
        if centers.shape != self.cluster_centers.shape:
            raise ValueError(f"center shape {tuple(centers.shape)} != {tuple(self.cluster_centers.shape)}")
        self.cluster_centers.copy_(centers)

    @staticmethod
    def target_distribution(q: torch.Tensor) -> torch.Tensor:
        weight = q.pow(2) / q.sum(dim=0, keepdim=True).clamp_min(1e-8)
        return weight / weight.sum(dim=1, keepdim=True).clamp_min(1e-8)
