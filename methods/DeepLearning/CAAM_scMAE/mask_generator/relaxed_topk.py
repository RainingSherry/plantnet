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

    max_k = int(k_i.max().item()) if k_i.numel() else 0
    khot = torch.zeros_like(logits32)
    onehot_approx = []
    for step in range(max_k):
        active = (k_i > step).float().view(-1, 1)
        masked_logits = logits32 + torch.log(torch.clamp(1.0 - khot, min=1.0e-6))
        probs = torch.softmax(masked_logits / max(float(tau), 1.0e-6), dim=1)
        probs = probs * eligibility.float() * active
        khot = khot + probs
        onehot_approx.append(probs)
    mask_soft = torch.stack(onehot_approx, dim=0).sum(dim=0) if onehot_approx else torch.zeros_like(logits32)
    mask_soft = mask_soft * eligibility.float()
    mask_st = mask_hard + mask_soft - mask_soft.detach()
    return mask_hard.to(logits.dtype), mask_soft.to(logits.dtype), mask_st.to(logits.dtype)

