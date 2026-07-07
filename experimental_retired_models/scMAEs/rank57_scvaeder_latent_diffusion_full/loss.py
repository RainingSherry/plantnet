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


class ScVAEDerLoss(nn.Module):
    """scMAE objective plus latent diffusion denoising and optional raw-count ZINB."""

    def __init__(
        self,
        masked_data_weight: float = 0.75,
        mask_loss_weight: float = 0.7,
        diffusion_weight: float = 0.02,
        kl_weight: float = 0.0001,
        zinb_weight: float = 0.005,
        variance_weight: float = 0.0,
    ):
        super().__init__()
        self.masked_data_weight = float(masked_data_weight)
        self.mask_loss_weight = float(mask_loss_weight)
        self.diffusion_weight = float(diffusion_weight)
        self.kl_weight = float(kl_weight)
        self.zinb_weight = float(zinb_weight)
        self.variance_weight = float(variance_weight)

    def forward(
        self,
        model,
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
        t = torch.rand(outputs["mu"].shape[0], device=outputs["mu"].device).clamp(0.02, 0.98)
        alpha_bar = torch.cos(0.5 * torch.pi * t).pow(2).clamp(1e-4, 0.9999)
        noise = torch.randn_like(outputs["mu"])
        z_noisy = torch.sqrt(alpha_bar).view(-1, 1) * outputs["mu"].detach() + torch.sqrt(1.0 - alpha_bar).view(-1, 1) * noise
        noise_pred = model.predict_noise(z_noisy, t)
        diffusion = F.mse_loss(noise_pred, noise)
        zinb = outputs["embedding"].new_tensor(0.0)
        if raw_counts is not None and size_factor is not None and self.zinb_weight > 0.0:
            zinb = zinb_nll(raw_counts, outputs["nb_mean"], outputs["nb_theta"], outputs["nb_dropout_logits"], size_factor)
        std = torch.sqrt(outputs["embedding"].var(dim=0, unbiased=False) + 1e-4)
        variance_loss = F.relu(0.5 - std).mean()
        scmae_loss = (1.0 - self.mask_loss_weight) * rec_loss + self.mask_loss_weight * mask_loss
        total = scmae_loss + self.diffusion_weight * diffusion + self.kl_weight * kl + self.zinb_weight * zinb + self.variance_weight * variance_loss
        return total, {
            "loss": float(total.detach().cpu()),
            "scmae_loss": float(scmae_loss.detach().cpu()),
            "recon_loss": float(rec_loss.detach().cpu()),
            "mask_loss": float(mask_loss.detach().cpu()),
            "diffusion_loss": float(diffusion.detach().cpu()),
            "kl_loss": float(kl.detach().cpu()),
            "zinb_loss": float(zinb.detach().cpu()),
            "variance_loss": float(variance_loss.detach().cpu()),
        }, rec_raw.detach().mean(dim=0)
