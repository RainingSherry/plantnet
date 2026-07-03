from __future__ import annotations

import torch


def mask_gini(values: torch.Tensor) -> torch.Tensor:
    x = values.flatten().float()
    if x.numel() == 0:
        return torch.tensor(0.0, device=values.device)
    sorted_x, _ = torch.sort(x)
    n = sorted_x.numel()
    denom = sorted_x.sum().clamp_min(1.0e-8)
    index = torch.arange(1, n + 1, device=x.device, dtype=torch.float32)
    return (torch.sum((2 * index - n - 1) * sorted_x) / (n * denom)).clamp_min(0.0)


def generator_regularization(
    mask_soft: torch.Tensor,
    logits: torch.Tensor,
    eligibility: torch.Tensor,
    x: torch.Tensor,
    replacement: torch.Tensor,
    *,
    tau: float,
    distortion_min: float,
    distortion_max: float,
) -> dict[str, torch.Tensor]:
    eps = 1.0e-8
    per_gene = mask_soft.mean(dim=0)
    r = per_gene / (per_gene.sum() + eps)
    u = eligibility.float().sum(dim=0)
    u = u / (u.sum() + eps)
    coverage = torch.sum(r * (torch.log(r + eps) - torch.log(u + eps)))

    sigma = x.std(dim=0, unbiased=False).clamp_min(eps)
    distortion = (replacement.detach() - x).abs() / sigma.view(1, -1)
    mean_dist = (mask_soft * distortion).sum() / (mask_soft.sum() + eps)
    distortion_loss = torch.relu(torch.tensor(float(distortion_min), device=x.device) - mean_dist).pow(2)
    distortion_loss = distortion_loss + torch.relu(mean_dist - torch.tensor(float(distortion_max), device=x.device)).pow(2)

    masked_logits = logits.float().masked_fill(~eligibility.bool(), -1.0e9)
    p = torch.softmax(masked_logits / max(float(tau), eps), dim=1) * eligibility.float()
    p = p / (p.sum(dim=1, keepdim=True) + eps)
    n_eligible = eligibility.float().sum(dim=1).clamp_min(2.0)
    entropy = -(p * torch.log(p + eps)).sum(dim=1) / torch.log(n_eligible)
    entropy_loss = -entropy.mean()
    return {
        "coverage_loss": coverage,
        "distortion_loss": distortion_loss,
        "entropy_loss": entropy_loss,
        "per_gene_mask_rate": per_gene.detach(),
        "per_cell_mask_rate": mask_soft.mean(dim=1).detach(),
        "mask_entropy": entropy.mean().detach(),
        "mask_gini": mask_gini(per_gene).detach(),
        "distortion_mean": mean_dist.detach(),
    }

