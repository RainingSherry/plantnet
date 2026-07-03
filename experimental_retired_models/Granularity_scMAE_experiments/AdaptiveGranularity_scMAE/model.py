from __future__ import annotations

import torch
from torch import nn


class AdaptiveGranularityScMAE(nn.Module):
    """rank13 DEC scMAE backbone + optional rank29 SVD-anchor fusion.

    - expr encoder: Dropout->Linear(g,256)->LN->Mish->Linear(256,h)->LN->Mish->Linear(h,h)
    - optional anchor fusion (rank29): fuse(cat[expr, enc(anchor)]) + 0.25*enc(anchor)
      gives a stable geometric prior tied to the raw SVD manifold.
    - mask_predictor, decoder(concat[latent, mask_logits]) as in scMAE.
    - trainable DEC cluster_centers; student_q = Student-t soft assignment.
    """

    def __init__(self, num_genes: int, n_clusters: int, hidden_size: int = 128,
                 anchor_dim: int = 0, dropout: float = 0.05):
        super().__init__()
        self.num_genes = int(num_genes)
        self.n_clusters = int(n_clusters)
        self.hidden_size = int(hidden_size)
        self.anchor_dim = int(anchor_dim)
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
        if self.anchor_dim > 0:
            self.anchor_encoder = nn.Sequential(
                nn.Linear(self.anchor_dim, hidden_size),
                nn.LayerNorm(hidden_size),
                nn.Mish(),
            )
            self.fusion = nn.Sequential(
                nn.Linear(hidden_size * 2, hidden_size),
                nn.LayerNorm(hidden_size),
                nn.Mish(),
            )
        self.mask_predictor = nn.Linear(hidden_size, self.num_genes)
        self.decoder = nn.Linear(hidden_size + self.num_genes, self.num_genes)
        self.cluster_centers = nn.Parameter(torch.randn(self.n_clusters, hidden_size) * 0.02)

    def encode(self, x: torch.Tensor, anchor: torch.Tensor | None = None) -> torch.Tensor:
        expr = self.encoder(x)
        if self.anchor_dim > 0 and anchor is not None:
            anc = self.anchor_encoder(anchor)
            return self.fusion(torch.cat([expr, anc], dim=1)) + 0.25 * anc
        return expr

    def forward(self, x: torch.Tensor, anchor: torch.Tensor | None = None) -> dict[str, torch.Tensor]:
        latent = self.encode(x, anchor)
        mask_logits = self.mask_predictor(latent)
        reconstruction = self.decoder(torch.cat([latent, mask_logits], dim=1))
        q = self.student_q(latent)
        return {"latent": latent, "mask_logits": mask_logits, "reconstruction": reconstruction, "cluster_q": q}

    def student_q(self, latent: torch.Tensor) -> torch.Tensor:
        dist = torch.cdist(latent, self.cluster_centers).pow(2)
        q = 1.0 / (1.0 + dist)
        return q / q.sum(dim=1, keepdim=True).clamp_min(1e-8)

    @torch.no_grad()
    def feature(self, x: torch.Tensor, anchor: torch.Tensor | None = None) -> torch.Tensor:
        return self.encode(x, anchor)

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
    def sharpen(q: torch.Tensor) -> torch.Tensor:
        """DEC target distribution (hard-EM sharpening of responsibilities q)."""
        weight = q.pow(2) / q.sum(dim=0, keepdim=True).clamp_min(1e-8)
        return weight / weight.sum(dim=1, keepdim=True).clamp_min(1e-8)
