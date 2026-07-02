from __future__ import annotations

import torch
from torch import nn


class MaskFeatGeneModuleScMAE(nn.Module):
    """
    MaskFeat-inspired scMAE with gene-module semantic feature prediction.

    The model keeps the scMAE mask prediction and masked expression decoder, then
    adds a semantic head that predicts module-level features for modules touched
    by the mask. Masked genes are replaced by a learnable per-gene sentinel.
    """

    def __init__(
        self,
        num_genes: int,
        num_modules: int = 32,
        module_feature_dim: int = 3,
        hidden_size: int = 128,
        dropout: float = 0.0,
    ):
        super().__init__()
        self.num_genes = int(num_genes)
        self.num_modules = int(num_modules)
        self.module_feature_dim = int(module_feature_dim)
        self.hidden_size = int(hidden_size)

        self.mask_value = nn.Parameter(torch.zeros(self.num_genes))
        self.encoder = nn.Sequential(
            nn.Dropout(p=float(dropout)),
            nn.Linear(self.num_genes, 256),
            nn.LayerNorm(256),
            nn.Mish(),
            nn.Linear(256, self.hidden_size),
            nn.LayerNorm(self.hidden_size),
            nn.Mish(),
            nn.Linear(self.hidden_size, self.hidden_size),
        )
        self.mask_predictor = nn.Linear(self.hidden_size, self.num_genes)
        self.expression_decoder = nn.Sequential(
            nn.Linear(self.hidden_size + self.num_genes, 256),
            nn.LayerNorm(256),
            nn.Mish(),
            nn.Linear(256, self.num_genes),
        )
        self.semantic_head = nn.Sequential(
            nn.Linear(self.hidden_size, self.hidden_size),
            nn.Mish(),
            nn.Linear(self.hidden_size, self.num_modules * self.module_feature_dim),
        )

    def sample_module_mask(
        self,
        x: torch.Tensor,
        module_ids: torch.Tensor,
        mask_prob: float,
        module_mask_prob: float,
    ) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
        random_mask = torch.bernoulli(float(mask_prob) * torch.ones_like(x)).bool()
        module_draw = torch.bernoulli(
            float(module_mask_prob) * torch.ones((x.shape[0], self.num_modules), device=x.device)
        ).bool()
        module_mask = module_draw[:, module_ids.long()]
        mask = (random_mask | module_mask).float()

        empty_rows = mask.sum(dim=1) == 0
        if bool(empty_rows.any()):
            cols = torch.randint(0, x.shape[1], (int(empty_rows.sum()),), device=x.device)
            mask[empty_rows, cols] = 1.0

        touched_modules = torch.zeros((x.shape[0], self.num_modules), dtype=torch.float32, device=x.device)
        touched_modules.scatter_reduce_(
            1,
            module_ids.long().unsqueeze(0).expand(x.shape[0], -1),
            mask,
            reduce="amax",
            include_self=False,
        )
        sentinel = self.mask_value.to(dtype=x.dtype, device=x.device).unsqueeze(0)
        corrupted = torch.where(mask.bool(), sentinel.expand_as(x), x)
        return corrupted, mask, touched_modules

    def forward(self, x: torch.Tensor) -> dict[str, torch.Tensor]:
        latent = self.encoder(x)
        mask_logits = self.mask_predictor(latent)
        reconstruction = self.expression_decoder(torch.cat([latent, mask_logits], dim=1))
        semantic = self.semantic_head(latent).view(x.shape[0], self.num_modules, self.module_feature_dim)
        return {
            "latent": latent,
            "mask_logits": mask_logits,
            "reconstruction": reconstruction,
            "semantic": semantic,
        }

    @torch.no_grad()
    def feature(self, x: torch.Tensor) -> torch.Tensor:
        return self.encoder(x)

