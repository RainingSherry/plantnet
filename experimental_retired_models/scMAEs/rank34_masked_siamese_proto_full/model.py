from __future__ import annotations

import math

import torch
from torch import nn
import torch.nn.functional as F


class SiameseBranch(nn.Module):
    def __init__(self, num_genes: int, hidden_size: int = 128, projection_dim: int = 128, dropout: float = 0.05):
        super().__init__()
        self.encoder = nn.Sequential(
            nn.Dropout(dropout),
            nn.Linear(num_genes, 256),
            nn.LayerNorm(256),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(256, hidden_size),
            nn.LayerNorm(hidden_size),
        )
        self.projector = nn.Sequential(
            nn.Linear(hidden_size, hidden_size * 2),
            nn.GELU(),
            nn.LayerNorm(hidden_size * 2),
            nn.Linear(hidden_size * 2, projection_dim),
        )

    def forward(self, x: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
        z = self.encoder(x)
        p = self.projector(z)
        return z, p


class MaskedSiameseProtoScMAE(nn.Module):
    """MSN-style masked-anchor/target-prototype scMAE."""

    def __init__(
        self,
        num_genes: int,
        hidden_size: int = 128,
        projection_dim: int = 128,
        num_prototypes: int = 128,
        dropout: float = 0.05,
    ):
        super().__init__()
        self.num_genes = int(num_genes)
        self.hidden_size = int(hidden_size)
        self.projection_dim = int(projection_dim)
        self.num_prototypes = int(num_prototypes)
        self.student = SiameseBranch(num_genes, hidden_size, projection_dim, dropout)
        self.teacher = SiameseBranch(num_genes, hidden_size, projection_dim, dropout)
        self.prototypes = nn.Parameter(torch.empty(num_prototypes, projection_dim))
        self.mask_predictor = nn.Linear(hidden_size, num_genes)
        self.decoder = nn.Sequential(
            nn.Linear(hidden_size + num_genes, 256),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(256, num_genes),
        )
        self.reset_parameters()

    def reset_parameters(self) -> None:
        for module in self.modules():
            if isinstance(module, nn.Linear):
                nn.init.xavier_uniform_(module.weight)
                if module.bias is not None:
                    nn.init.zeros_(module.bias)
            elif isinstance(module, nn.LayerNorm):
                nn.init.ones_(module.weight)
                nn.init.zeros_(module.bias)
        bound = 1.0 / math.sqrt(max(1, self.projection_dim))
        nn.init.uniform_(self.prototypes, -bound, bound)
        self.reset_teacher()

    @torch.no_grad()
    def reset_teacher(self) -> None:
        self.teacher.load_state_dict(self.student.state_dict(), strict=True)
        for param in self.teacher.parameters():
            param.requires_grad_(False)

    @torch.no_grad()
    def update_teacher(self, momentum: float) -> None:
        for student_param, teacher_param in zip(self.student.parameters(), self.teacher.parameters()):
            teacher_param.data.mul_(float(momentum)).add_(student_param.detach().data, alpha=1.0 - float(momentum))

    def forward_student(self, x: torch.Tensor) -> dict[str, torch.Tensor]:
        latent, projection = self.student(x)
        mask_logits = self.mask_predictor(latent)
        reconstruction = self.decoder(torch.cat([latent, mask_logits], dim=1))
        return {
            "latent": latent,
            "projection": projection,
            "mask_logits": mask_logits,
            "reconstruction": reconstruction,
        }

    @torch.no_grad()
    def forward_teacher(self, x: torch.Tensor) -> dict[str, torch.Tensor]:
        self.teacher.eval()
        latent, projection = self.teacher(x)
        return {"latent": latent, "projection": projection}

    @torch.no_grad()
    def feature(self, x: torch.Tensor, use_teacher: bool = False) -> torch.Tensor:
        branch = self.teacher if use_teacher else self.student
        return branch(x)[0]

    def prototype_probs(self, projection: torch.Tensor, temperature: float) -> torch.Tensor:
        query = F.normalize(projection, dim=-1)
        supports = F.normalize(self.prototypes, dim=-1)
        return F.softmax(query @ supports.t() / float(temperature), dim=-1)

    def mask_view(self, x: torch.Tensor, mask_prob: float) -> tuple[torch.Tensor, torch.Tensor]:
        mask = (torch.rand_like(x) < float(mask_prob)).float()
        empty = mask.sum(dim=1) == 0
        if bool(empty.any()):
            cols = torch.randint(0, x.shape[1], (int(empty.sum()),), device=x.device)
            mask[empty, cols] = 1.0
        return x.masked_fill(mask.bool(), 0.0), mask

    def weak_view(self, x: torch.Tensor, noise_std: float = 0.03, dropout_prob: float = 0.02) -> torch.Tensor:
        y = x
        if noise_std > 0:
            y = y + torch.randn_like(y) * float(noise_std)
        if dropout_prob > 0:
            keep = (torch.rand_like(y) >= float(dropout_prob)).float()
            y = y * keep
        return y
