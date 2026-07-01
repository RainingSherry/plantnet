from __future__ import annotations

import torch
import torch.nn.functional as F
from torch import nn


def nb_nll(counts: torch.Tensor, mean: torch.Tensor, theta: torch.Tensor, size_factor: torch.Tensor) -> torch.Tensor:
    mu = mean * size_factor.view(-1, 1)
    theta = theta.view(1, -1)
    log_prob = (
        torch.lgamma(counts + theta)
        - torch.lgamma(theta)
        - torch.lgamma(counts + 1.0)
        + theta * (torch.log(theta + 1e-8) - torch.log(theta + mu + 1e-8))
        + counts * (torch.log(mu + 1e-8) - torch.log(theta + mu + 1e-8))
    )
    return -log_prob.mean()


class MultiFacetVAELoss(nn.Module):
    """scMAE loss plus multi-facet VAE uncertainty and optional NB raw-count likelihood."""

    def __init__(
        self,
        masked_data_weight: float = 0.75,
        mask_loss_weight: float = 0.7,
        kl_weight: float = 0.01,
        nb_weight: float = 0.02,
        entropy_weight: float = 0.002,
        variance_weight: float = 0.0,
    ):
        super().__init__()
        self.masked_data_weight = float(masked_data_weight)
        self.mask_loss_weight = float(mask_loss_weight)
        self.kl_weight = float(kl_weight)
        self.nb_weight = float(nb_weight)
        self.entropy_weight = float(entropy_weight)
        self.variance_weight = float(variance_weight)

    def forward(
        self,
        model,
        outputs: dict,
        target_expr: torch.Tensor,
        mask: torch.Tensor,
        raw_counts: torch.Tensor | None,
        size_factor: torch.Tensor | None,
    ) -> tuple[torch.Tensor, dict[str, float], torch.Tensor]:
        rec_raw = F.smooth_l1_loss(outputs["reconstruction"], target_expr, reduction="none")
        weights = mask * self.masked_data_weight + (1.0 - mask) * (1.0 - self.masked_data_weight)
        rec_loss = (weights * rec_raw).mean()
        mask_loss = F.binary_cross_entropy_with_logits(outputs["mask_logits"], mask)
        kl = outputs["embedding"].new_tensor(0.0)
        entropy = outputs["embedding"].new_tensor(0.0)
        for j, (mu, logvar, z, resp) in enumerate(zip(outputs["mus"], outputs["logvars"], outputs["zs"], outputs["responsibilities"])):
            log_q = -0.5 * (((z - mu).pow(2) / torch.exp(logvar)) + logvar + torch.log(torch.tensor(2.0 * torch.pi, device=z.device))).sum(dim=1)
            diff = z[:, None, :] - model.mog_mu[j][None, :, :]
            comp_logvar = model.mog_logvar[j][None, :, :].clamp(-8.0, 6.0)
            log_comp = -0.5 * (diff.pow(2) / torch.exp(comp_logvar) + comp_logvar + torch.log(torch.tensor(2.0 * torch.pi, device=z.device))).sum(dim=-1)
            log_p = torch.logsumexp(log_comp + F.log_softmax(model.mog_logits[j], dim=0).view(1, -1), dim=1)
            kl = kl + (log_q - log_p).mean()
            entropy = entropy + (-(resp.clamp_min(1e-8) * resp.clamp_min(1e-8).log()).sum(dim=1).mean())
        nb_loss = outputs["embedding"].new_tensor(0.0)
        if raw_counts is not None and size_factor is not None and self.nb_weight > 0:
            nb_loss = nb_nll(raw_counts, outputs["nb_mean"], outputs["nb_theta"], size_factor)
        std = torch.sqrt(outputs["embedding"].var(dim=0, unbiased=False) + 1e-4)
        variance_loss = F.relu(0.5 - std).mean()
        scmae_loss = (1.0 - self.mask_loss_weight) * rec_loss + self.mask_loss_weight * mask_loss
        total = scmae_loss + self.kl_weight * kl + self.nb_weight * nb_loss - self.entropy_weight * entropy + self.variance_weight * variance_loss
        return total, {
            "loss": float(total.detach().cpu()),
            "scmae_loss": float(scmae_loss.detach().cpu()),
            "recon_loss": float(rec_loss.detach().cpu()),
            "mask_loss": float(mask_loss.detach().cpu()),
            "kl_loss": float(kl.detach().cpu()),
            "nb_loss": float(nb_loss.detach().cpu()),
            "facet_entropy": float(entropy.detach().cpu()),
            "variance_loss": float(variance_loss.detach().cpu()),
        }, rec_raw.detach().mean(dim=0)
