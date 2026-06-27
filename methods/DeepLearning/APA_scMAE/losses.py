from __future__ import annotations

import torch
import torch.nn.functional as F


def student_losses(
    x: torch.Tensor,
    x_recon: torch.Tensor,
    pred_mask_logits: torch.Tensor,
    effective_mask: torch.Tensor,
    *,
    masked_data_weight: float,
    gamma: float,
) -> dict[str, torch.Tensor]:
    effective_mask = effective_mask.float()
    weights = effective_mask * float(masked_data_weight) + (1.0 - effective_mask) * (1.0 - float(masked_data_weight))
    rec = (weights * (x_recon - x).pow(2)).mean()
    mask = F.binary_cross_entropy_with_logits(pred_mask_logits, effective_mask, reduction="mean")
    loss = (1.0 - float(gamma)) * rec + float(gamma) * mask
    return {"loss_student": loss, "loss_rec": rec, "loss_mask": mask}


def mask_entropy(logits: torch.Tensor) -> torch.Tensor:
    probs = torch.softmax(logits, dim=1)
    entropy = -(probs * torch.log(probs + 1.0e-8)).sum(dim=1)
    return -entropy.mean()


def balance_loss(mask_soft: torch.Tensor) -> torch.Tensor:
    per_gene = mask_soft.mean(dim=0)
    return per_gene.var(unbiased=False)


def distortion_loss(mask_soft: torch.Tensor, delta: torch.Tensor) -> torch.Tensor:
    selected_delta = (mask_soft * delta).sum() / (mask_soft.sum() + 1.0e-8)
    all_delta = delta.mean().detach()
    return -selected_delta / (all_delta + 1.0e-8)


def coverage_loss(mask_soft: torch.Tensor, target_ratio: float) -> torch.Tensor:
    return (mask_soft.mean() - float(target_ratio)).pow(2)


def generator_losses(
    x: torch.Tensor,
    x_recon: torch.Tensor,
    effective_mask_st: torch.Tensor,
    logits: torch.Tensor,
    mask_soft: torch.Tensor,
    delta: torch.Tensor,
    *,
    mask_ratio: float,
    lambda_entropy: float,
    lambda_balance: float,
    lambda_distortion: float,
    lambda_coverage: float,
) -> dict[str, torch.Tensor]:
    err = (x_recon - x).pow(2)
    hard = -(effective_mask_st * err).sum() / (effective_mask_st.sum() + 1.0e-8)
    ent = mask_entropy(logits)
    bal = balance_loss(mask_soft)
    dist = distortion_loss(mask_soft, delta)
    cov = coverage_loss(mask_soft, mask_ratio)
    total = hard + float(lambda_entropy) * ent + float(lambda_balance) * bal + float(lambda_distortion) * dist + float(lambda_coverage) * cov
    return {
        "loss_generator": total,
        "loss_hard": hard,
        "loss_entropy": ent,
        "loss_balance": bal,
        "loss_distortion": dist,
        "loss_coverage": cov,
    }


def mask_diagnostics(mask_hard: torch.Tensor, effective_mask: torch.Tensor, logits: torch.Tensor) -> dict[str, float]:
    mask = mask_hard.detach().float()
    eff = effective_mask.detach().float()
    probs = torch.softmax(logits.detach(), dim=1)
    entropy = -(probs * torch.log(probs + 1.0e-8)).sum(dim=1).mean()
    per_gene = mask.mean(dim=0)
    sorted_vals = torch.sort(per_gene).values
    n = sorted_vals.numel()
    if n == 0 or float(sorted_vals.sum().cpu()) == 0.0:
        gini = torch.tensor(0.0, device=mask.device)
    else:
        idx = torch.arange(1, n + 1, device=mask.device, dtype=sorted_vals.dtype)
        gini = (2.0 * (idx * sorted_vals).sum() / (n * sorted_vals.sum())) - (n + 1.0) / n
    return {
        "mask_ratio": float(mask.mean().cpu()),
        "effective_mask_ratio": float(eff.mean().cpu()),
        "zero_to_zero_rate": float((mask * (1.0 - eff)).sum().cpu() / (mask.sum().cpu() + 1.0e-8)),
        "mask_entropy": float(entropy.cpu()),
        "mask_gini": float(gini.cpu()),
        "top_gene_concentration": float(per_gene.max().cpu()) if per_gene.numel() else 0.0,
    }
