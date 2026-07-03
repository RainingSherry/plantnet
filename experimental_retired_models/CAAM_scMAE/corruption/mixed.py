from __future__ import annotations

import torch


class MixedCorruption:
    corruption_type = "mixed"

    def corrupt(self, x: torch.Tensor, mask: torch.Tensor, metadata: dict | None = None):
        r = torch.rand_like(x)
        mask_swap = mask.bool() & (r < 0.34)
        mask_zero = mask.bool() & (r >= 0.34) & (r < 0.67) & (x > 0)
        mask_noise = mask.bool() & ~(mask_swap | mask_zero)
        perm = torch.argsort(torch.rand(x.shape, device=x.device), dim=0)
        shuffled = torch.gather(x, dim=0, index=perm)
        sigma = 0.1 * x.std(dim=0, unbiased=False).clamp_min(1.0e-8)
        noisy = torch.clamp(x + torch.randn_like(x) * sigma.view(1, -1), min=0.0)
        x_tilde = x.clone()
        x_tilde = torch.where(mask_swap, shuffled, x_tilde)
        x_tilde = torch.where(mask_zero, torch.zeros_like(x_tilde), x_tilde)
        x_tilde = torch.where(mask_noise, noisy, x_tilde)
        info = {
            "selected_mask_rate": float(mask.float().mean().detach().cpu()),
            "effective_changed_rate": float((x_tilde != x).float().mean().detach().cpu()),
            "corruption_type": self.corruption_type,
            "mask_swap": mask_swap.float(),
            "mask_zero": mask_zero.float(),
            "mask_noise": mask_noise.float(),
        }
        return x_tilde, x_tilde, info

