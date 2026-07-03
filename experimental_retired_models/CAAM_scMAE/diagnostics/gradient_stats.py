from __future__ import annotations

import torch.nn as nn

from methods.DeepLearning.CAAM_scMAE.models.common import grad_norm


def collect_gradient_stats(student: nn.Module, generator: nn.Module | None = None) -> dict:
    return {
        "student_grad_norm": grad_norm(student),
        "generator_grad_norm": grad_norm(generator) if generator is not None else 0.0,
        "encoder_grad_norm": grad_norm(student.encoder),
        "mask_head_grad_norm": grad_norm(student.mask_head),
        "decoder_grad_norm": grad_norm(student.decoder),
    }

