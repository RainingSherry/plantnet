from __future__ import annotations

import torch
import torch.nn.functional as F
from torch import nn


def zinb_negative_log_likelihood(
    x: torch.Tensor,
    mean: torch.Tensor,
    theta: torch.Tensor,
    dropout_logits: torch.Tensor,
    eps: float = 1e-8,
) -> torch.Tensor:
    x = x.clamp_min(0.0)
    mean = mean.clamp_min(eps)
    theta = theta.clamp_min(eps)
    nb_case = (
        torch.lgamma(theta + x)
        - torch.lgamma(theta)
        - torch.lgamma(x + 1.0)
        + theta * (torch.log(theta + eps) - torch.log(theta + mean + eps))
        + x * (torch.log(mean + eps) - torch.log(theta + mean + eps))
    )
    log_pi = F.logsigmoid(dropout_logits)
    log_not_pi = F.logsigmoid(-dropout_logits)
    zero_nb = theta * (torch.log(theta + eps) - torch.log(theta + mean + eps))
    zero_case = torch.logsumexp(torch.stack([log_pi, log_not_pi + zero_nb], dim=0), dim=0)
    nonzero_case = log_not_pi + nb_case
    log_prob = torch.where(x < eps, zero_case, nonzero_case)
    return -log_prob.mean()


class ScVGAEZINBLoss(nn.Module):
    def __init__(
        self,
        masked_data_weight: float = 0.75,
        mask_weight: float = 0.7,
        zinb_weight: float = 0.01,
        kl_weight: float = 0.001,
    ):
        super().__init__()
        self.masked_data_weight = float(masked_data_weight)
        self.mask_weight = float(mask_weight)
        self.zinb_weight = float(zinb_weight)
        self.kl_weight = float(kl_weight)

    def forward(
        self,
        outputs: dict[str, torch.Tensor],
        log_expr_target: torch.Tensor,
        count_target: torch.Tensor,
        mask: torch.Tensor,
    ) -> tuple[torch.Tensor, dict[str, float]]:
        weights = mask * self.masked_data_weight + (1.0 - mask) * (1.0 - self.masked_data_weight)
        reconstruction_loss = torch.mul(
            weights,
            F.mse_loss(outputs["reconstruction"], log_expr_target, reduction="none"),
        ).mean()
        mask_loss = F.binary_cross_entropy_with_logits(outputs["mask_logits"], mask.float())
        zinb_loss = zinb_negative_log_likelihood(
            count_target,
            outputs["zinb_mean"],
            outputs["zinb_theta"],
            outputs["zinb_dropout_logits"],
        )
        kl_loss = -0.5 * torch.mean(1.0 + outputs["logvar"] - outputs["latent"].pow(2) - outputs["logvar"].exp())
        scmae_loss = (1.0 - self.mask_weight) * reconstruction_loss + self.mask_weight * mask_loss
        loss = scmae_loss + self.zinb_weight * zinb_loss + self.kl_weight * kl_loss
        return loss, {
            "loss": float(loss.detach().cpu()),
            "scmae_loss": float(scmae_loss.detach().cpu()),
            "reconstruction_loss": float(reconstruction_loss.detach().cpu()),
            "mask_loss": float(mask_loss.detach().cpu()),
            "zinb_loss": float(zinb_loss.detach().cpu()),
            "kl_loss": float(kl_loss.detach().cpu()),
        }

