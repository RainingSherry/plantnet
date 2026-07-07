from __future__ import annotations

import math

import torch
from torch import nn
import torch.nn.functional as F


class MaskedSiameseProtoLoss(nn.Module):
    def __init__(
        self,
        masked_data_weight: float = 0.75,
        mask_weight: float = 0.60,
        proto_weight: float = 0.12,
        memax_weight: float = 0.04,
        entropy_weight: float = 0.00,
        variance_weight: float = 0.01,
        temperature: float = 0.10,
        sharpen_temperature: float = 0.25,
    ):
        super().__init__()
        self.masked_data_weight = float(masked_data_weight)
        self.mask_weight = float(mask_weight)
        self.proto_weight = float(proto_weight)
        self.memax_weight = float(memax_weight)
        self.entropy_weight = float(entropy_weight)
        self.variance_weight = float(variance_weight)
        self.temperature = float(temperature)
        self.sharpen_temperature = float(sharpen_temperature)

    @staticmethod
    def sharpen(probs: torch.Tensor, temperature: float) -> torch.Tensor:
        p = probs.clamp_min(1e-8).pow(1.0 / float(temperature))
        return p / p.sum(dim=1, keepdim=True).clamp_min(1e-8)

    @staticmethod
    def variance_loss(z: torch.Tensor) -> torch.Tensor:
        std = torch.sqrt(z.var(dim=0) + 1e-4)
        return F.relu(0.5 - std).mean()

    def forward(
        self,
        model,
        student: dict[str, torch.Tensor],
        teacher: dict[str, torch.Tensor],
        target_expr: torch.Tensor,
        mask: torch.Tensor,
        proto_weight_scale: float = 1.0,
    ) -> tuple[torch.Tensor, dict[str, float]]:
        weights = mask * self.masked_data_weight + (1.0 - mask) * (1.0 - self.masked_data_weight)
        rec = (weights * F.smooth_l1_loss(student["reconstruction"], target_expr, reduction="none")).mean()
        mask_loss = F.binary_cross_entropy_with_logits(student["mask_logits"], mask.float())
        scmae = (1.0 - self.mask_weight) * rec + self.mask_weight * mask_loss

        anchor_probs = model.prototype_probs(student["projection"], self.temperature)
        with torch.no_grad():
            target_probs = model.prototype_probs(teacher["projection"], self.temperature)
            target_probs = self.sharpen(target_probs, self.sharpen_temperature)
        proto_ce = torch.mean(torch.sum(-target_probs * anchor_probs.clamp_min(1e-8).log(), dim=1))
        avg_probs = anchor_probs.mean(dim=0).clamp_min(1e-8)
        memax = torch.sum(avg_probs * avg_probs.log()) + math.log(float(avg_probs.numel()))
        entropy = torch.mean(torch.sum(-anchor_probs * anchor_probs.clamp_min(1e-8).log(), dim=1))
        var_loss = self.variance_loss(student["latent"])
        scale = float(proto_weight_scale)
        total = scmae + scale * self.proto_weight * proto_ce + self.memax_weight * memax + self.entropy_weight * entropy + self.variance_weight * var_loss
        with torch.no_grad():
            target_conf = target_probs.max(dim=1).values.mean()
            used_proto = float(torch.unique(target_probs.argmax(dim=1)).numel())
        return total, {
            "loss": float(total.detach().cpu()),
            "scmae_loss": float(scmae.detach().cpu()),
            "reconstruction_loss": float(rec.detach().cpu()),
            "mask_loss": float(mask_loss.detach().cpu()),
            "proto_ce_loss": float(proto_ce.detach().cpu()),
            "memax_loss": float(memax.detach().cpu()),
            "anchor_entropy": float(entropy.detach().cpu()),
            "variance_loss": float(var_loss.detach().cpu()),
            "target_confidence": float(target_conf.detach().cpu()),
            "target_used_prototypes": used_proto,
            "proto_weight_scale": scale,
        }
