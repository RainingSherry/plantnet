from __future__ import annotations

import torch
from torch import nn
import torch.nn.functional as F


class MaskGITGeneTokenLoss(nn.Module):
    def __init__(
        self,
        masked_data_weight: float = 0.75,
        mask_weight: float = 0.55,
        token_weight: float = 0.16,
        replaced_weight: float = 0.08,
        confidence_weight: float = 0.01,
    ):
        super().__init__()
        self.masked_data_weight = float(masked_data_weight)
        self.mask_weight = float(mask_weight)
        self.token_weight = float(token_weight)
        self.replaced_weight = float(replaced_weight)
        self.confidence_weight = float(confidence_weight)

    def forward(
        self,
        out: dict[str, torch.Tensor],
        target: torch.Tensor,
        token_target: torch.Tensor,
        mask: torch.Tensor,
        replaced: torch.Tensor,
    ) -> tuple[torch.Tensor, dict[str, float]]:
        weights = mask * self.masked_data_weight + (1.0 - mask) * (1.0 - self.masked_data_weight)
        rec = (weights * F.smooth_l1_loss(out["reconstruction"], target, reduction="none")).mean()
        mask_loss = F.binary_cross_entropy_with_logits(out["mask_logits"], mask.float())
        token_ce = F.cross_entropy(out["token_logits"].permute(0, 2, 1), token_target, reduction="none")
        token_loss = (token_ce * mask).sum() / mask.sum().clamp_min(1.0)
        replaced_loss = F.binary_cross_entropy_with_logits(out["replaced_logits"], replaced.float())
        conf = out["token_confidence"]
        confidence_loss = torch.relu(0.55 - (conf * mask).sum() / mask.sum().clamp_min(1.0))
        scmae = (1.0 - self.mask_weight) * rec + self.mask_weight * mask_loss
        total = scmae + self.token_weight * token_loss + self.replaced_weight * replaced_loss + self.confidence_weight * confidence_loss
        return total, {
            "loss": float(total.detach().cpu()),
            "scmae_loss": float(scmae.detach().cpu()),
            "reconstruction_loss": float(rec.detach().cpu()),
            "mask_loss": float(mask_loss.detach().cpu()),
            "token_loss": float(token_loss.detach().cpu()),
            "replaced_loss": float(replaced_loss.detach().cpu()),
            "confidence_loss": float(confidence_loss.detach().cpu()),
            "masked_token_confidence": float(((conf * mask).sum() / mask.sum().clamp_min(1.0)).detach().cpu()),
            "replaced_fraction": float(replaced.detach().mean().cpu()),
        }
