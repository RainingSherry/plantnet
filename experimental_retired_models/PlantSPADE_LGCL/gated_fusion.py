from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Tuple, Optional

import numpy as np
import scipy.sparse as sp
import torch
import torch.nn as nn
import torch.nn.functional as F

from .negative_sampling import NegativeSampler, NegativeSamplerConfig


class GatedFusion(nn.Module):
    """Lightweight adaptive fusion for frozen cell-level LGCL views."""

    def __init__(
        self,
        latent_dim: int,
        hidden_dim: int | None = None,
        dropout: float = 0.05,
        extra_dim: int = 0,
    ):
        super().__init__()
        self.latent_dim = int(latent_dim)
        self.extra_dim = int(extra_dim)
        hidden = int(hidden_dim) if hidden_dim and hidden_dim > 0 else max(8, self.latent_dim)
        self.gate_mlp = nn.Sequential(
            nn.Linear(3 * self.latent_dim + self.extra_dim, hidden),
            nn.GELU(),
            nn.Dropout(float(dropout)),
            nn.Linear(hidden, 3),
        )
        final = self.gate_mlp[-1]
        nn.init.zeros_(final.weight)
        nn.init.zeros_(final.bias)

    def forward(
        self,
        z_support: torch.Tensor,
        z_global: torch.Tensor,
        z_attention: torch.Tensor,
        batch_cells: torch.Tensor | None = None,
        extra: Optional[torch.Tensor] = None,
    ) -> Tuple[torch.Tensor, torch.Tensor]:
        if batch_cells is not None:
            support = z_support[batch_cells]
            global_view = z_global[batch_cells]
            attention = z_attention[batch_cells]
            extra_batch = extra[batch_cells] if extra is not None else None
        else:
            support = z_support
            global_view = z_global
            attention = z_attention
            extra_batch = extra

        stacked = torch.stack((support, global_view, attention), dim=1)
        gate_in = torch.cat((support, global_view, attention), dim=1)
        if extra_batch is not None:
            gate_in = torch.cat((gate_in, extra_batch), dim=1)
        gate_logits = self.gate_mlp(gate_in)
        gate = torch.softmax(gate_logits, dim=1)
        fused = torch.sum(gate.unsqueeze(-1) * stacked, dim=1)
        return F.normalize(fused, dim=1), gate


@dataclass
class GatedFusionConfig:
    epochs: int = 20
    batch_size: int = 2048
    pairs_per_epoch: int = 65536
    lr: float = 5e-3
    weight_decay: float = 1e-4
    hidden_dim: int = 0
    dropout: float = 0.05
    temperature: float = 0.2
    bpr_weight: float = 1.0
    contrastive_weight: float = 0.05
    consistency_weight: float = 0.05
    seed: int = 42
    negative_sampler: str = "random_zero"
    negative_neighbor_k: int = 15

    # Optional enhancements inspired by dynamic fusion literature
    use_cell_stats: bool = True
    gate_entropy_weight: float = 0.0
    gate_balance_weight: float = 0.0


class FrozenSupportPairSampler:
    def __init__(
        self,
        support: sp.csr_matrix,
        seed: int,
        negative_sampler: str,
        negative_neighbor_k: int,
    ):
        self.support = support.tocsr(copy=True)
        self.support.data = np.ones_like(self.support.data, dtype=np.float32)
        self.support.eliminate_zeros()
        self.support.sort_indices()
        coo = self.support.tocoo()
        self.edge_rows = coo.row.astype(np.int64, copy=False)
        self.edge_cols = coo.col.astype(np.int64, copy=False)
        if self.edge_rows.size == 0:
            raise ValueError("Support matrix has no non-zero cell-gene edges.")
        self.rng = np.random.default_rng(seed)
        self.negative_sampler = NegativeSampler(
            self.support,
            NegativeSamplerConfig(mode=negative_sampler, seed=seed, neighbor_k=negative_neighbor_k),
        )

    def sample(self, n_pairs: int) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        edge_idx = self.rng.integers(0, self.edge_rows.size, size=int(n_pairs), dtype=np.int64)
        cells = self.edge_rows[edge_idx]
        pos_genes = self.edge_cols[edge_idx]
        neg_genes = self.negative_sampler.sample(cells)
        return cells, pos_genes, neg_genes


def _as_float_tensor(array: np.ndarray, device: torch.device) -> torch.Tensor:
    return torch.as_tensor(np.asarray(array, dtype=np.float32), dtype=torch.float32, device=device)


def _symmetric_infonce(anchor: torch.Tensor, target: torch.Tensor, temperature: float) -> torch.Tensor:
    if anchor.shape[0] <= 1:
        return torch.zeros((), dtype=anchor.dtype, device=anchor.device)
    logits = F.normalize(anchor, dim=1) @ F.normalize(target, dim=1).T / float(temperature)
    labels = torch.arange(logits.shape[0], device=logits.device)
    return 0.5 * (F.cross_entropy(logits, labels) + F.cross_entropy(logits.T, labels))


def train_gated_fusion(
    z_support: np.ndarray,
    z_global: np.ndarray,
    z_attention: np.ndarray,
    gene_embedding: np.ndarray,
    support: sp.csr_matrix,
    device: torch.device,
    config: GatedFusionConfig,
) -> Tuple[np.ndarray, np.ndarray, Dict[str, list]]:
    """Train only the gate MLP against frozen support/global/attention views."""

    if z_support.shape != z_global.shape or z_support.shape != z_attention.shape:
        raise ValueError(
            "Gated fusion expects same-shaped cell embeddings: "
            f"support={z_support.shape}, global={z_global.shape}, attention={z_attention.shape}"
        )
    if z_support.ndim != 2:
        raise ValueError(f"Cell embeddings must be 2-D, got {z_support.shape}")

    latent_dim = int(z_support.shape[1])

    # Optional per-cell statistics as extra gating input (helps robustness across datasets)
    cell_stats_t: Optional[torch.Tensor] = None
    extra_dim = 0
    if config.use_cell_stats:
        deg = np.diff(support.tocsr().indptr).astype(np.float32)
        lib = np.asarray(support.sum(axis=1)).ravel().astype(np.float32)
        # log-scale then standardize
        stats = np.stack((np.log1p(deg), np.log1p(lib)), axis=1)
        stats = (stats - stats.mean(axis=0, keepdims=True)) / (stats.std(axis=0, keepdims=True) + 1e-6)
        cell_stats_t = _as_float_tensor(stats, device)
        extra_dim = int(cell_stats_t.shape[1])

    model = GatedFusion(
        latent_dim=latent_dim,
        hidden_dim=config.hidden_dim,
        dropout=config.dropout,
        extra_dim=extra_dim,
    ).to(device)
    z_support_t = _as_float_tensor(z_support, device)
    z_global_t = _as_float_tensor(z_global, device)
    z_attention_t = _as_float_tensor(z_attention, device)
    gene_t = F.normalize(_as_float_tensor(gene_embedding, device), dim=1)
    support_target_t = F.normalize(z_support_t, dim=1)
    global_target_t = F.normalize(z_global_t, dim=1)
    attention_target_t = F.normalize(z_attention_t, dim=1)
    mean_target_t = F.normalize(z_support_t + z_global_t + z_attention_t, dim=1)

    sampler = FrozenSupportPairSampler(
        support=support,
        seed=config.seed,
        negative_sampler=config.negative_sampler,
        negative_neighbor_k=config.negative_neighbor_k,
    )
    optimizer = torch.optim.AdamW(model.parameters(), lr=config.lr, weight_decay=config.weight_decay)
    rng = np.random.default_rng(config.seed + 313)
    n_cells = int(z_support.shape[0])
    pairs_per_epoch = max(1, int(config.pairs_per_epoch))
    batch_size = max(1, min(int(config.batch_size), n_cells))
    epochs = max(0, int(config.epochs))
    history: Dict[str, list] = {
        "loss": [],
        "bpr": [],
        "contrastive": [],
        "consistency": [],
        "gate_support_mean": [],
        "gate_global_mean": [],
        "gate_attention_mean": [],
    }

    for _ in range(epochs):
        model.train()
        cells_np, pos_np, neg_np = sampler.sample(pairs_per_epoch)
        cells = torch.as_tensor(cells_np, dtype=torch.long, device=device)
        pos = torch.as_tensor(pos_np, dtype=torch.long, device=device)
        neg = torch.as_tensor(neg_np, dtype=torch.long, device=device)
        fused_pairs, gate_pairs = model(z_support_t, z_global_t, z_attention_t, batch_cells=cells, extra=cell_stats_t)
        pos_scores = torch.sum(fused_pairs * gene_t[pos], dim=1)
        neg_scores = torch.sum(fused_pairs * gene_t[neg], dim=1)
        bpr = F.softplus(-(pos_scores - neg_scores)).mean()

        contrastive_cells = torch.as_tensor(
            rng.choice(n_cells, size=batch_size, replace=False),
            dtype=torch.long,
            device=device,
        )
        fused_batch, _ = model(
            z_support_t,
            z_global_t,
            z_attention_t,
            batch_cells=contrastive_cells,
            extra=cell_stats_t,
        )
        contrastive = (
            _symmetric_infonce(fused_batch, support_target_t[contrastive_cells], config.temperature)
            + _symmetric_infonce(fused_batch, global_target_t[contrastive_cells], config.temperature)
            + _symmetric_infonce(fused_batch, attention_target_t[contrastive_cells], config.temperature)
        ) / 3.0
        consistency = (
            1.0
            - F.cosine_similarity(fused_batch, mean_target_t[contrastive_cells], dim=1)
        ).mean()
        loss = (
            float(config.bpr_weight) * bpr
            + float(config.contrastive_weight) * contrastive
            + float(config.consistency_weight) * consistency
        )

        # Optional gate regularization to reduce view collapse
        if config.gate_entropy_weight and config.gate_entropy_weight > 0.0:
            gate_entropy = -(gate_pairs * torch.log(gate_pairs.clamp_min(1e-12))).sum(dim=1).mean()
            loss = loss - float(config.gate_entropy_weight) * gate_entropy
        if config.gate_balance_weight and config.gate_balance_weight > 0.0:
            gate_mean_t = gate_pairs.mean(dim=0)
            target = torch.full_like(gate_mean_t, 1.0 / 3.0)
            gate_balance = torch.sum((gate_mean_t - target) ** 2)
            loss = loss + float(config.gate_balance_weight) * gate_balance
        optimizer.zero_grad(set_to_none=True)
        loss.backward()
        optimizer.step()

        gate_mean = gate_pairs.detach().mean(dim=0).cpu().numpy()
        history["loss"].append(float(loss.detach().cpu()))
        history["bpr"].append(float(bpr.detach().cpu()))
        history["contrastive"].append(float(contrastive.detach().cpu()))
        history["consistency"].append(float(consistency.detach().cpu()))
        history["gate_support_mean"].append(float(gate_mean[0]))
        history["gate_global_mean"].append(float(gate_mean[1]))
        history["gate_attention_mean"].append(float(gate_mean[2]))

    model.eval()
    with torch.no_grad():
        fused, gate = model(z_support_t, z_global_t, z_attention_t, batch_cells=None, extra=cell_stats_t)
    return (
        fused.detach().cpu().numpy().astype(np.float32),
        gate.detach().cpu().numpy().astype(np.float32),
        history,
    )
