from __future__ import annotations

import numpy as np
import torch
import torch.nn as nn


class GeneModuleTokenizer(nn.Module):
    def __init__(self, assignment: np.ndarray, token_dim: int) -> None:
        super().__init__()
        assignment = np.asarray(assignment, dtype=np.float32)
        denom = assignment.sum(axis=0, keepdims=True)
        denom[denom == 0.0] = 1.0
        normalized = assignment / denom
        self.register_buffer("assignment", torch.as_tensor(normalized, dtype=torch.float32), persistent=True)
        self.n_modules = int(normalized.shape[1])
        self.value_mlp = nn.Sequential(nn.Linear(1, token_dim), nn.Mish(), nn.Linear(token_dim, token_dim))
        self.module_embedding = nn.Parameter(torch.empty(self.n_modules, token_dim))
        nn.init.normal_(self.module_embedding, std=0.02)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        u = x @ self.assignment
        value = self.value_mlp(u.unsqueeze(-1))
        return value + self.module_embedding.view(1, self.n_modules, -1)

