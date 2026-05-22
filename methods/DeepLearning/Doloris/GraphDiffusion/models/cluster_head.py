import torch
import torch.nn as nn
import torch.nn.functional as F


class ClusterHead(nn.Module):
    def __init__(self, input_dim: int, n_clusters: int):
        super().__init__()
        self.input_dim = input_dim
        self.n_clusters = n_clusters
        self.prototypes = nn.Parameter(torch.randn(n_clusters, input_dim) * 0.02)

    def forward(self, z: torch.Tensor):
        distances = torch.cdist(z, self.prototypes, p=2) ** 2
        logits = -distances
        probs = F.softmax(logits, dim=-1)
        return logits, probs

    def target_distribution(self, probs: torch.Tensor) -> torch.Tensor:
        weight = probs ** 2 / probs.sum(dim=0, keepdim=True).clamp_min(1e-12)
        return weight / weight.sum(dim=1, keepdim=True).clamp_min(1e-12)
