from __future__ import annotations

import torch
from torch import nn


class JOAOMaskPolicyScMAE(nn.Module):
    """Independent scMAE body with augmentation-aware projection heads."""

    def __init__(self, num_genes: int, hidden_size: int = 128, dropout: float = 0.0, n_policies: int = 4):
        super().__init__()
        self.num_genes = int(num_genes)
        self.hidden_size = int(hidden_size)
        self.n_policies = int(n_policies)
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
        self.decoder = nn.Linear(self.hidden_size + self.num_genes, self.num_genes)
        self.policy_heads = nn.ModuleList(
            [
                nn.Sequential(nn.Linear(self.hidden_size, self.hidden_size), nn.Mish(), nn.Linear(self.hidden_size, self.hidden_size))
                for _ in range(self.n_policies)
            ]
        )

    def forward(self, x: torch.Tensor, policy_id: int = 0) -> dict[str, torch.Tensor]:
        latent = self.encoder(x)
        mask_logits = self.mask_predictor(latent)
        reconstruction = self.decoder(torch.cat([latent, mask_logits], dim=1))
        projection = self.policy_heads[int(policy_id)](latent)
        return {
            "latent": latent,
            "mask_logits": mask_logits,
            "reconstruction": reconstruction,
            "projection": projection,
        }

    @torch.no_grad()
    def feature(self, x: torch.Tensor) -> torch.Tensor:
        return self.encoder(x)

    def corrupt(
        self,
        x: torch.Tensor,
        log_expr: torch.Tensor,
        module_ids: torch.Tensor,
        policy_id: int,
        mask_prob: float,
        module_mask_prob: float,
    ) -> tuple[torch.Tensor, torch.Tensor]:
        policy_id = int(policy_id)
        if policy_id == 0:
            should_swap = torch.bernoulli(float(mask_prob) * torch.ones_like(x)).bool()
            repl = x[torch.randperm(x.shape[0], device=x.device)] if x.shape[0] > 1 else x
            corrupted = torch.where(should_swap, repl, x)
            mask = (corrupted != x).float()
        elif policy_id == 1:
            mask = torch.bernoulli(float(mask_prob) * torch.ones_like(x)).float()
            corrupted = torch.where(mask.bool(), self.mask_value[None, :].expand_as(x), x)
        elif policy_id == 2:
            n_modules = int(module_ids.max().item()) + 1
            module_draw = torch.bernoulli(float(module_mask_prob) * torch.ones((x.shape[0], n_modules), device=x.device)).bool()
            mask = module_draw[:, module_ids.long()].float()
            corrupted = torch.where(mask.bool(), self.mask_value[None, :].expand_as(x), x)
        else:
            nonzero_weight = (log_expr > 1e-8).float() * 0.75 + 0.25
            mask = torch.bernoulli((float(mask_prob) * nonzero_weight).clamp(max=0.95)).float()
            corrupted = torch.where(mask.bool(), torch.zeros_like(x), x)
        empty = mask.sum(dim=1) == 0
        if bool(empty.any()):
            cols = torch.randint(0, x.shape[1], (int(empty.sum()),), device=x.device)
            mask[empty, cols] = 1.0
            corrupted[empty, cols] = self.mask_value[cols]
        return corrupted, mask

