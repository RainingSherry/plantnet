from __future__ import annotations

import torch


def sample_gumbel(shape, device) -> torch.Tensor:
    u = torch.rand(shape, device=device, dtype=torch.float32).clamp_(1.0e-6, 1.0 - 1.0e-6)
    return -torch.log(-torch.log(u))


def relaxed_topk_straight_through(
    logits: torch.Tensor,
    k_i: torch.Tensor,
    tau: float,
    eligibility: torch.Tensor,
    add_gumbel: bool = True,
) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
    logits32 = logits.float().masked_fill(~eligibility.bool(), -1.0e9)
    if add_gumbel:
        logits32 = logits32 + sample_gumbel(logits32.shape, logits32.device)
    b, g = logits32.shape
    k_i = k_i.to(device=logits32.device, dtype=torch.long)
    mask_hard = torch.zeros((b, g), dtype=torch.float32, device=logits32.device)
    for row in range(b):
        k = int(k_i[row].item())
        if k > 0:
            valid = int(eligibility[row].sum().item())
            k = min(k, valid)
            if k > 0:
                mask_hard[row, torch.topk(logits32[row], k=k).indices] = 1.0

    tau_value = max(float(tau), 1.0e-6)
    eligibility_f = eligibility.float()
    valid_counts = eligibility_f.sum(dim=1)
    target = torch.minimum(k_i.float(), valid_counts).view(-1, 1)

    valid_min = torch.where(eligibility, logits32, torch.full_like(logits32, float("inf"))).amin(dim=1, keepdim=True)
    valid_max = torch.where(eligibility, logits32, torch.full_like(logits32, float("-inf"))).amax(dim=1, keepdim=True)
    valid_min = torch.where(torch.isfinite(valid_min), valid_min, torch.zeros_like(valid_min))
    valid_max = torch.where(torch.isfinite(valid_max), valid_max, torch.zeros_like(valid_max))

    margin = 50.0 * tau_value
    low = valid_min - margin
    high = valid_max + margin
    for _ in range(40):
        mid = (low + high) * 0.5
        probs = torch.sigmoid((logits32 - mid) / tau_value) * eligibility_f
        too_many = probs.sum(dim=1, keepdim=True) > target
        low = torch.where(too_many, mid, low)
        high = torch.where(too_many, high, mid)

    threshold = (low + high) * 0.5
    mask_soft = torch.sigmoid((logits32 - threshold) / tau_value) * eligibility_f
    mask_soft = torch.where(target <= 0, torch.zeros_like(mask_soft), mask_soft)
    mask_soft = torch.where(target >= valid_counts.view(-1, 1), eligibility_f, mask_soft)
    mask_st = mask_hard + mask_soft - mask_soft.detach()
    return mask_hard.to(logits.dtype), mask_soft.to(logits.dtype), mask_st.to(logits.dtype)
