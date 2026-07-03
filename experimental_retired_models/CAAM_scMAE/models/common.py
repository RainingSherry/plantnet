from __future__ import annotations

import contextlib
from typing import Iterator

import torch
import torch.nn as nn


def trainable_parameter_count(module: nn.Module) -> int:
    return int(sum(p.numel() for p in module.parameters() if p.requires_grad))


def grad_norm(module: nn.Module) -> float:
    total = 0.0
    for p in module.parameters():
        if p.grad is not None:
            total += float(p.grad.detach().pow(2).sum().cpu())
    return float(total ** 0.5)


@contextlib.contextmanager
def freeze_module(module: nn.Module) -> Iterator[None]:
    states = [p.requires_grad for p in module.parameters()]
    was_training = module.training
    for p in module.parameters():
        p.requires_grad_(False)
    try:
        yield
    finally:
        for p, state in zip(module.parameters(), states):
            p.requires_grad_(state)
        module.train(was_training)

