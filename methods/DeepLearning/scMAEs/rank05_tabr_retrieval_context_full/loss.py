from __future__ import annotations

import torch
import torch.nn.functional as F
from torch import nn


class TabRContextLoss(nn.Module):
    def __init__(
        self,
        masked_data_weight: float = 0.75,
        mask_weight: float = 0.7,
        retrieval_weight: float = 0.1,
    ):
        super().__init__()
        self.masked_data_weight = float(masked_data_weight)
        self.mask_weight = float(mask_weight)
        self.retrieval_weight = float(retrieval_weight)

    def forward(self, outputs: dict[str, torch.Tensor], target: torch.Tensor, mask: torch.Tensor) -> tuple[torch.Tensor, dict[str, float]]:
        weights = mask * self.masked_data_weight + (1.0 - mask) * (1.0 - self.masked_data_weight)
        reconstruction_loss = torch.mul(weights, F.mse_loss(outputs["reconstruction"], target, reduction="none")).mean()
        mask_loss = F.binary_cross_entropy_with_logits(outputs["mask_logits"], mask.float())
        retrieval_loss = F.mse_loss(
            F.normalize(outputs["latent"], dim=1),
            F.normalize(outputs["context_latent"].detach(), dim=1),
        )
        loss = (1.0 - self.mask_weight) * reconstruction_loss + self.mask_weight * mask_loss + self.retrieval_weight * retrieval_loss
        parts = {
            "loss": float(loss.detach().cpu()),
            "reconstruction_loss": float(reconstruction_loss.detach().cpu()),
            "mask_loss": float(mask_loss.detach().cpu()),
            "retrieval_loss": float(retrieval_loss.detach().cpu()),
        }
        return loss, parts

