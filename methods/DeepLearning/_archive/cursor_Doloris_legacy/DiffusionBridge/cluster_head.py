"""
Cluster head module: DEC-style clustering for diffusion bridge embeddings.

Implements a differentiable clustering head using the Student-t distribution
(similarly to DEC - Deep Embedded Clustering). The head learns cluster
prototypes and produces soft assignments that can be refined via self-training.

This enables end-to-end optimization of the bridge representation toward
cluster-separability.
"""
from __future__ import annotations

import torch
import torch.nn as nn
import torch.nn.functional as F


class ClusterHead(nn.Module):
    """
    DEC-style cluster head for learning cluster-separable representations.

    Uses a Student-t kernel (q distribution) and KL-divergence self-training
    (p distribution) to iteratively refine cluster boundaries.

    Architecture:
        encoder: MLP that transforms bridge embeddings to prototype space
        prototypes: learnable cluster centers (n_clusters, hidden_dim)
    """

    def __init__(
        self,
        input_dim: int,
        n_clusters: int,
        hidden_dim: int = 256,
        alpha: float = 1.0,
        init: str = "kmeans",
    ):
        """
        Args:
            input_dim: dimension of bridge target embeddings
            n_clusters: number of clusters
            hidden_dim: dimension of prototype space
            alpha: Student-t degrees of freedom parameter
            init: initialization strategy ('kmeans' or 'random')
        """
        super().__init__()
        self.alpha = alpha
        self.n_clusters = n_clusters

        self.encoder = nn.Sequential(
            nn.Linear(input_dim, hidden_dim),
            nn.LayerNorm(hidden_dim),
            nn.SiLU(),
            nn.Linear(hidden_dim, hidden_dim),
            nn.LayerNorm(hidden_dim),
            nn.SiLU(),
        )

        # Learnable cluster prototypes, initialized to small random values
        self.prototypes = nn.Parameter(torch.randn(n_clusters, hidden_dim) * 0.02)

    def forward(self, z: torch.Tensor) -> torch.Tensor:
        """
        Compute soft assignment probabilities (q distribution).

        Uses Student-t kernel (1 / (1 + d^2 / alpha)) to measure similarity
        between encoded embeddings and cluster prototypes.

        Args:
            z: bridge target embeddings (batch_size, input_dim)

        Returns:
            Soft assignment matrix q (batch_size, n_clusters)
        """
        h = self.encoder(z)
        # Squared Euclidean distance: \|h - mu_k\|^2
        dist = torch.cdist(h, self.prototypes, p=2).pow(2)
        # Student-t kernel: q_{ik} = (1 + d^2/alpha)^(-(alpha+1)/2)
        q = 1.0 / (1.0 + dist / self.alpha)
        q = q.pow((self.alpha + 1.0) / 2.0)
        q = q / q.sum(dim=1, keepdim=True).clamp_min(1e-8)
        return q

    @staticmethod
    def target_distribution(q: torch.Tensor) -> torch.Tensor:
        """
        Compute target distribution (p distribution) for KL-divergence self-training.

        Sharpens the soft assignments: high-confidence assignments get upweighted,
        low-confidence assignments get downweighted.

        Formula: p_{ik} = q_{ik}^2 / f_k / sum_j(q_{ij}^2 / f_j)
        where f_k = sum_j(q_{jk}) is the cluster soft frequency.

        Args:
            q: soft assignments from forward() (batch_size, n_clusters)

        Returns:
            Target distribution p (batch_size, n_clusters)
        """
        weight = q.pow(2) / q.sum(dim=0, keepdim=True).clamp_min(1e-8)
        return weight / weight.sum(dim=1, keepdim=True).clamp_min(1e-8)

    def hard_assign(self, z: torch.Tensor) -> torch.Tensor:
        """Get hard cluster assignments by argmax."""
        return self.forward(z).argmax(dim=1)

    def kl_divergence(self, q: torch.Tensor, p: torch.Tensor) -> torch.Tensor:
        """KL(q || p) for self-training."""
        return F.kl_div(
            q.clamp_min(1e-8).log(), p, reduction="batchmean"
        )

    def entropy(self, q: torch.Tensor) -> torch.Tensor:
        """Entropy of soft assignments: penalizes overconfident or diffuse distributions."""
        return -(q * q.clamp_min(1e-8).log()).sum(dim=1).mean()

    def initialize_from_kmeans(self, z: torch.Tensor, labels: torch.Tensor | None = None):
        """
        Initialize prototypes from KMeans clustering on embeddings.

        If labels are provided, prototypes are set to the mean of each cluster.
        Otherwise, run KMeans directly on the embeddings.

        Args:
            z: embeddings (batch_size, input_dim)
            labels: optional ground-truth cluster labels
        """
        if labels is not None:
            unique_labels = torch.unique(labels)
            with torch.no_grad():
                for cluster_id in unique_labels:
                    mask = labels == cluster_id
                    if mask.sum() > 0:
                        self.prototypes[cluster_id] = z[mask].mean(dim=0)
        else:
            from sklearn.cluster import KMeans

            z_np = z.detach().cpu().numpy()
            kmeans = KMeans(n_clusters=self.n_clusters, random_state=42, n_init=20)
            kmeans.fit(z_np)
            self.prototypes.data.copy_(torch.from_numpy(kmeans.cluster_centers_).float())

    def center_loss(self, z: torch.Tensor, q: torch.Tensor) -> torch.Tensor:
        """
        Compute center-based loss: encourages embeddings to move toward prototypes.

        Loss = sum_i(sum_k(q_{ik} * \|h(x_i) - mu_k\|^2))
        """
        h = self.encoder(z)
        dist = torch.cdist(h, self.prototypes, p=2).pow(2)
        return (dist * q).sum(dim=1).mean()
