from __future__ import annotations

import torch


def generator_loss(
    loss_rec_masked: torch.Tensor,
    loss_mask: torch.Tensor,
    regularizers: dict[str, torch.Tensor],
    *,
    beta_mask_loss: float,
    lambda_coverage: float,
    lambda_distortion: float,
    lambda_entropy: float,
) -> torch.Tensor:
    # Negative reconstruction/mask loss means generator maximizes student difficulty.
    return (
        -loss_rec_masked
        -float(beta_mask_loss) * loss_mask
        + float(lambda_coverage) * regularizers["coverage_loss"]
        + float(lambda_distortion) * regularizers["distortion_loss"]
        + float(lambda_entropy) * regularizers["entropy_loss"]
    )

