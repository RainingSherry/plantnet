from typing import Sequence, Tuple

import torch
import torch.nn as nn
import torch.nn.functional as F

from models.cluster_head import ClusterHead
from models.gene_gat_encoder import GeneGATEncoder
from models.mask_diffusion_refiner import MaskDiffusionRefiner
from models.support_pooling import SupportPooling


class GraphDiffusionModel(nn.Module):
    def __init__(
        self,
        gene_input_dim: int,
        gene_hidden_dim: int,
        gene_output_dim: int,
        cell_dim: int,
        module_count: int,
        n_clusters: int,
        gat_layers: int = 2,
        dropout: float = 0.1,
    ):
        super().__init__()
        self.gene_encoder = GeneGATEncoder(
            input_dim=gene_input_dim,
            hidden_dim=gene_hidden_dim,
            output_dim=gene_output_dim,
            num_layers=gat_layers,
            dropout=dropout,
        )
        self.pooling = SupportPooling(
            gene_dim=gene_output_dim,
            cell_dim=cell_dim,
            module_count=module_count,
        )
        self.refiner = MaskDiffusionRefiner(cell_dim=cell_dim, steps=3, dropout=dropout)
        self.cluster_head = ClusterHead(input_dim=cell_dim, n_clusters=n_clusters)
        self.module_decoder = nn.Linear(cell_dim, module_count)

    def forward(
        self,
        gene_features: torch.Tensor,
        edge_index: torch.Tensor,
        edge_weight: torch.Tensor,
        support_indices: Sequence[torch.Tensor],
        support_weights: Sequence[torch.Tensor],
    ) -> Tuple[torch.Tensor, torch.Tensor, torch.Tensor, torch.Tensor, torch.Tensor, torch.Tensor]:
        gene_embeddings = self.gene_encoder(gene_features, edge_index=edge_index, edge_weight=edge_weight)
        cell_embeddings, pooled_modules = self.pooling(gene_embeddings, support_indices, support_weights)
        refined_embeddings = self.refiner(cell_embeddings)
        logits, probs = self.cluster_head(refined_embeddings)
        module_activation = self.module_decoder(refined_embeddings) + pooled_modules
        return gene_embeddings, cell_embeddings, refined_embeddings, logits, probs, pooled_modules, module_activation

    def compute_losses(
        self,
        refined_embeddings: torch.Tensor,
        logits: torch.Tensor,
        module_activation: torch.Tensor,
        pooled_modules: torch.Tensor,
        cluster_head: ClusterHead,
        consistency_weight: float = 0.2,
        module_weight: float = 0.1,
    ):
        probs = F.softmax(logits, dim=-1)
        target = cluster_head.target_distribution(probs).detach()
        cluster_loss = F.kl_div(probs.clamp_min(1e-8).log(), target, reduction="batchmean")
        module_loss = F.mse_loss(module_activation, pooled_modules.detach())
        variance_loss = -refined_embeddings.var(dim=0).mean()
        total = cluster_loss + module_weight * module_loss + consistency_weight * variance_loss
        return {
            "total": total,
            "cluster": cluster_loss,
            "module": module_loss,
            "variance": variance_loss,
        }
