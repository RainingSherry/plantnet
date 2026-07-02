from __future__ import annotations

import torch
import torch.nn.functional as F
from torch import nn


def nb_log_prob(counts: torch.Tensor, mean: torch.Tensor, theta: torch.Tensor, size_factor: torch.Tensor) -> torch.Tensor:
    mu = mean * size_factor.view(-1, 1)
    theta = theta.view(1, -1)
    return (
        torch.lgamma(counts + theta)
        - torch.lgamma(theta)
        - torch.lgamma(counts + 1.0)
        + theta * (torch.log(theta + 1e-8) - torch.log(theta + mu + 1e-8))
        + counts * (torch.log(mu + 1e-8) - torch.log(theta + mu + 1e-8))
    )


def zinb_nll(counts: torch.Tensor, mean: torch.Tensor, theta: torch.Tensor, pi_logits: torch.Tensor, size_factor: torch.Tensor) -> torch.Tensor:
    nb_lp = nb_log_prob(counts, mean, theta, size_factor)
    log_pi = F.logsigmoid(pi_logits)
    log_not_pi = F.logsigmoid(-pi_logits)
    zero_case = torch.logaddexp(log_pi, log_not_pi + nb_lp)
    nonzero_case = log_not_pi + nb_lp
    log_prob = torch.where(counts <= 1e-8, zero_case, nonzero_case)
    return -log_prob.mean()


def gaussian_kernel(x: torch.Tensor, y: torch.Tensor, bandwidths: tuple[float, ...] = (0.5, 1.0, 2.0, 4.0)) -> torch.Tensor:
    dist = torch.cdist(x, y, p=2).pow(2)
    out = x.new_zeros(dist.shape)
    for bw in bandwidths:
        out = out + torch.exp(-dist / (2.0 * bw * bw))
    return out / float(len(bandwidths))


def mmd_rbf(x: torch.Tensor, y: torch.Tensor) -> torch.Tensor:
    kxx = gaussian_kernel(x, x)
    kyy = gaussian_kernel(y, y)
    kxy = gaussian_kernel(x, y)
    return kxx.mean() + kyy.mean() - 2.0 * kxy.mean()


class ScInfoVAELoss(nn.Module):
    """scMAE objective plus InfoVAE MMD/MI and optional raw-count ZINB likelihood."""

    def __init__(
        self,
        masked_data_weight: float = 0.75,
        mask_loss_weight: float = 0.7,
        mmd_weight: float = 0.001,
        kl_weight: float = 0.0001,
        zinb_weight: float = 0.005,
        variance_weight: float = 0.0,
    ):
        super().__init__()
        self.masked_data_weight = float(masked_data_weight)
        self.mask_loss_weight = float(mask_loss_weight)
        self.mmd_weight = float(mmd_weight)
        self.kl_weight = float(kl_weight)
        self.zinb_weight = float(zinb_weight)
        self.variance_weight = float(variance_weight)

    def forward(
        self,
        outputs: dict[str, torch.Tensor],
        target_expr: torch.Tensor,
        mask: torch.Tensor,
        raw_counts: torch.Tensor | None,
        size_factor: torch.Tensor | None,
    ) -> tuple[torch.Tensor, dict[str, float], torch.Tensor]:
        rec_raw = F.mse_loss(outputs["reconstruction"], target_expr, reduction="none")
        weights = mask * self.masked_data_weight + (1.0 - mask) * (1.0 - self.masked_data_weight)
        rec_loss = (weights * rec_raw).mean()
        mask_loss = F.binary_cross_entropy_with_logits(outputs["mask_logits"], mask)
        kl = -0.5 * torch.mean(1.0 + outputs["logvar"] - outputs["mu"].pow(2) - outputs["logvar"].exp())
        prior = torch.randn_like(outputs["z_sample"])
        mmd = mmd_rbf(outputs["z_sample"], prior)
        zinb = outputs["embedding"].new_tensor(0.0)
        if raw_counts is not None and size_factor is not None and self.zinb_weight > 0.0:
            zinb = zinb_nll(raw_counts, outputs["nb_mean"], outputs["nb_theta"], outputs["nb_dropout_logits"], size_factor)
        std = torch.sqrt(outputs["embedding"].var(dim=0, unbiased=False) + 1e-4)
        variance_loss = F.relu(0.5 - std).mean()
        scmae_loss = (1.0 - self.mask_loss_weight) * rec_loss + self.mask_loss_weight * mask_loss
        total = scmae_loss + self.mmd_weight * mmd + self.kl_weight * kl + self.zinb_weight * zinb + self.variance_weight * variance_loss
        return total, {
            "loss": float(total.detach().cpu()),
            "scmae_loss": float(scmae_loss.detach().cpu()),
            "recon_loss": float(rec_loss.detach().cpu()),
            "mask_loss": float(mask_loss.detach().cpu()),
            "mmd_loss": float(mmd.detach().cpu()),
            "kl_loss": float(kl.detach().cpu()),
            "zinb_loss": float(zinb.detach().cpu()),
            "variance_loss": float(variance_loss.detach().cpu()),
        }, rec_raw.detach().mean(dim=0)
