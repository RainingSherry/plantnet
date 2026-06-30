from __future__ import annotations

import sys
from pathlib import Path

import torch
import torch.nn as nn

ROOT = Path(__file__).resolve().parents[5]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from methods.DeepLearning.scMAE.model import AutoEncoder

from .prototype import PrototypeHead


class MechanismScMAE(nn.Module):
    """Original scMAE plus optional lightweight prototype head."""

    def __init__(
        self,
        num_genes: int,
        hidden_size: int,
        dropout: float,
        masked_data_weight: float,
        mask_loss_weight: float,
        n_clusters: int,
        use_prototype_head: bool,
    ) -> None:
        super().__init__()
        self.backbone = AutoEncoder(
            num_genes=num_genes,
            hidden_size=hidden_size,
            dropout=dropout,
            masked_data_weight=masked_data_weight,
            mask_loss_weight=mask_loss_weight,
        )
        self.prototype_head = (
            PrototypeHead(n_clusters=n_clusters, dim=hidden_size)
            if use_prototype_head and n_clusters > 0
            else None
        )

    @property
    def has_prototypes(self) -> bool:
        return self.prototype_head is not None

    def scmae_loss(self, corrupted: torch.Tensor, target: torch.Tensor, mask: torch.Tensor):
        return self.backbone.loss_mask(corrupted, target, mask)

    def feature(self, x: torch.Tensor) -> torch.Tensor:
        return self.backbone.feature(x)

    def prototype_logits(self, z: torch.Tensor) -> torch.Tensor:
        if self.prototype_head is None:
            raise RuntimeError("Prototype head is disabled for this variant.")
        return self.prototype_head(z)

    def prototype_centers(self) -> torch.Tensor:
        if self.prototype_head is None:
            raise RuntimeError("Prototype head is disabled for this variant.")
        return self.prototype_head.centers

    def set_prototypes(self, centers) -> None:
        if self.prototype_head is None:
            raise RuntimeError("Prototype head is disabled for this variant.")
        self.prototype_head.set_centers(centers)

    def normalize_prototypes(self) -> None:
        if self.prototype_head is not None:
            self.prototype_head.normalize_centers()


def build_mechanism_scmae(num_genes: int, args, n_clusters: int) -> MechanismScMAE:
    return MechanismScMAE(
        num_genes=num_genes,
        hidden_size=args.hidden_size,
        dropout=args.dropout,
        masked_data_weight=args.masked_data_weight,
        mask_loss_weight=args.mask_loss_weight,
        n_clusters=n_clusters,
        use_prototype_head=bool(args.use_prototype or args.use_swav),
    )

