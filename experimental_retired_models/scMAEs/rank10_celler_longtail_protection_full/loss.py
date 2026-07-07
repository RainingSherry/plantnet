from __future__ import annotations

import torch
from torch import nn
import torch.nn.functional as F


class CellerLongTailLoss(nn.Module):
    """scMAE loss plus unsupervised Gaussian Inflation and HDM terms."""

    def __init__(
        self,
        masked_data_weight: float = 0.75,
        mask_weight: float = 0.65,
        token_weight: float = 0.08,
        prototype_weight: float = 0.04,
        hdm_weight: float = 0.02,
        consistency_weight: float = 0.03,
        balance_weight: float = 0.02,
        rare_recon_boost: float = 1.0,
        tail_delta: float = 0.25,
        hard_topk: int = 8,
    ):
        super().__init__()
        self.masked_data_weight = float(masked_data_weight)
        self.mask_weight = float(mask_weight)
        self.token_weight = float(token_weight)
        self.prototype_weight = float(prototype_weight)
        self.hdm_weight = float(hdm_weight)
        self.consistency_weight = float(consistency_weight)
        self.balance_weight = float(balance_weight)
        self.rare_recon_boost = float(rare_recon_boost)
        self.tail_delta = float(tail_delta)
        self.hard_topk = int(hard_topk)

    def forward(
        self,
        out: dict[str, torch.Tensor],
        weak_out: dict[str, torch.Tensor],
        target: torch.Tensor,
        token_target: torch.Tensor,
        mask: torch.Tensor,
        rare_risk: torch.Tensor,
        boundary_risk: torch.Tensor,
        prototype_prior: torch.Tensor,
    ) -> tuple[torch.Tensor, dict[str, float]]:
        sample_risk = rare_risk.clamp(0, 1)
        recon_boost = (1.0 + self.rare_recon_boost * sample_risk).unsqueeze(1)
        gene_weight = mask * self.masked_data_weight + (1.0 - mask) * (1.0 - self.masked_data_weight)
        rec = (recon_boost * gene_weight * F.smooth_l1_loss(out["reconstruction"], target, reduction="none")).mean()
        mask_loss = F.binary_cross_entropy_with_logits(out["mask_logits"], mask.float())
        scmae = (1.0 - self.mask_weight) * rec + self.mask_weight * mask_loss

        token_ce = F.cross_entropy(out["token_logits"].permute(0, 2, 1), token_target, reduction="none")
        token_loss = (token_ce * mask * recon_boost.squeeze(1).unsqueeze(1)).sum() / mask.sum().clamp_min(1.0)

        logits = out["proto_logits"]
        probs = F.softmax(logits.detach(), dim=1)
        pseudo = probs.argmax(dim=1)
        confidence = probs.max(dim=1).values
        core_weight = (confidence.detach() * (1.0 - boundary_risk).clamp(0, 1) * (1.0 - 0.5 * rare_risk).clamp(0.1, 1.0)).detach()

        prior = prototype_prior.to(logits.device).clamp_min(1.0)
        inflation = torch.log(prior.max() / prior)
        inflated_logits = logits + self.tail_delta * inflation.unsqueeze(0)
        proto_ce = F.cross_entropy(inflated_logits, pseudo, reduction="none")
        proto_loss = (proto_ce * core_weight).sum() / core_weight.sum().clamp_min(1.0)

        hard_logits = logits.masked_fill(F.one_hot(pseudo, logits.shape[1]).bool(), -1e9)
        k = min(self.hard_topk, max(1, logits.shape[1] - 1))
        hard_vals = hard_logits.topk(k, dim=1).values
        target_vals = logits.gather(1, pseudo[:, None])
        hdm = F.softplus(hard_vals - target_vals + 0.2).mean(dim=1)
        hdm_loss = (hdm * (rare_risk + boundary_risk + 0.25).detach()).mean()

        z1 = F.normalize(out["latent"], dim=1)
        z2 = F.normalize(weak_out["latent"], dim=1)
        consistency = ((z1 - z2).pow(2).sum(dim=1) * (1.0 - boundary_risk).detach()).mean()

        mean_prob = F.softmax(logits, dim=1).mean(dim=0).clamp_min(1e-8)
        uniform = torch.full_like(mean_prob, 1.0 / mean_prob.numel())
        balance = F.kl_div(mean_prob.log(), uniform, reduction="sum")

        total = (
            scmae
            + self.token_weight * token_loss
            + self.prototype_weight * proto_loss
            + self.hdm_weight * hdm_loss
            + self.consistency_weight * consistency
            + self.balance_weight * balance
        )
        return total, {
            "loss": float(total.detach().cpu()),
            "scmae_loss": float(scmae.detach().cpu()),
            "reconstruction_loss": float(rec.detach().cpu()),
            "mask_loss": float(mask_loss.detach().cpu()),
            "token_loss": float(token_loss.detach().cpu()),
            "prototype_loss": float(proto_loss.detach().cpu()),
            "hdm_loss": float(hdm_loss.detach().cpu()),
            "consistency_loss": float(consistency.detach().cpu()),
            "balance_loss": float(balance.detach().cpu()),
            "core_weight": float(core_weight.mean().detach().cpu()),
        }

