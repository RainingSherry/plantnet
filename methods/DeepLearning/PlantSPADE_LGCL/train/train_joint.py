from dataclasses import dataclass
from typing import Dict, Optional, Tuple

import numpy as np
import scipy.sparse as sp
import torch
import torch.nn as nn
import torch.nn.functional as F


def scipy_to_torch_sparse(matrix: sp.spmatrix, device: torch.device) -> torch.Tensor:
    coo = matrix.tocoo().astype(np.float32)
    indices = np.vstack((coo.row, coo.col)).astype(np.int64)
    return torch.sparse_coo_tensor(
        torch.from_numpy(indices),
        torch.from_numpy(coo.data),
        size=coo.shape,
        dtype=torch.float32,
        device=device,
    ).coalesce()


def normalized_bipartite_support(support: sp.csr_matrix) -> sp.csr_matrix:
    support = support.astype(np.float32).tocsr(copy=True)
    support.data = np.ones_like(support.data, dtype=np.float32)
    support.eliminate_zeros()
    support.sort_indices()
    coo = support.tocoo()
    row_deg = np.asarray(support.sum(axis=1)).ravel().astype(np.float32)
    col_deg = np.asarray(support.sum(axis=0)).ravel().astype(np.float32)
    denom = np.sqrt(row_deg[coo.row] * col_deg[coo.col])
    values = np.divide(1.0, denom, out=np.zeros_like(denom), where=denom > 0)
    norm = sp.coo_matrix((values.astype(np.float32), (coo.row, coo.col)), shape=support.shape)
    return norm.tocsr()


def sparse_dropout(matrix: torch.Tensor, p: float) -> torch.Tensor:
    if p <= 0.0:
        return matrix
    matrix = matrix.coalesce()
    values = F.dropout(matrix.values(), p=p, training=True)
    return torch.sparse_coo_tensor(matrix.indices(), values, matrix.shape, device=matrix.device).coalesce()


class PairSampler:
    def __init__(self, support: sp.csr_matrix, seed: int = 42):
        self.support = support.tocsr(copy=True)
        self.support.sort_indices()
        coo = self.support.tocoo()
        self.edge_rows = coo.row.astype(np.int64, copy=False)
        self.edge_cols = coo.col.astype(np.int64, copy=False)
        self.indptr = self.support.indptr.astype(np.int64, copy=False)
        self.indices = self.support.indices.astype(np.int64, copy=False)
        self.n_genes = self.support.shape[1]
        self.rng = np.random.default_rng(seed)
        if self.edge_rows.size == 0:
            raise ValueError("Support matrix has no non-zero cell-gene edges.")

    def _row_contains(self, row: int, col: int) -> bool:
        start = self.indptr[row]
        end = self.indptr[row + 1]
        row_indices = self.indices[start:end]
        pos = np.searchsorted(row_indices, col)
        return bool(pos < row_indices.size and row_indices[pos] == col)

    def _sample_negatives(self, cells: np.ndarray) -> np.ndarray:
        neg = self.rng.integers(0, self.n_genes, size=cells.shape[0], dtype=np.int64)
        for i, row in enumerate(cells):
            if self.indptr[row + 1] - self.indptr[row] >= self.n_genes:
                continue
            while self._row_contains(int(row), int(neg[i])):
                neg[i] = self.rng.integers(0, self.n_genes, dtype=np.int64)
        return neg

    def sample(self, n_pairs: int) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        edge_idx = self.rng.integers(0, self.edge_rows.size, size=n_pairs, dtype=np.int64)
        cells = self.edge_rows[edge_idx]
        pos_genes = self.edge_cols[edge_idx]
        neg_genes = self._sample_negatives(cells)
        return cells, pos_genes, neg_genes


class TopKSparseModuleLayer(nn.Module):
    def __init__(self, latent_dim: int, num_modules: int, top_k: int):
        super().__init__()
        self.num_modules = int(num_modules)
        self.top_k = int(top_k)
        self.module_queries = nn.Parameter(torch.empty(self.num_modules, latent_dim))
        nn.init.xavier_uniform_(self.module_queries)

    def forward(self, cell_emb: torch.Tensor, gene_emb: torch.Tensor) -> Dict[str, torch.Tensor]:
        k = min(self.top_k, gene_emb.shape[0])
        scores = self.module_queries @ gene_emb.T / np.sqrt(gene_emb.shape[1])
        top_scores, top_idx = torch.topk(scores, k=k, dim=1)
        weights = torch.softmax(top_scores, dim=1)
        module_emb = torch.sum(weights.unsqueeze(-1) * gene_emb[top_idx], dim=1)
        module_emb = F.normalize(module_emb, dim=1)
        activations = cell_emb @ module_emb.T
        return {
            "scores": scores,
            "top_scores": top_scores,
            "top_idx": top_idx,
            "top_weights": weights,
            "module_emb": module_emb,
            "activations": activations,
        }


class PlantSPADELGCL(nn.Module):
    def __init__(
        self,
        n_cells: int,
        n_genes: int,
        latent_dim: int,
        adj_norm: torch.Tensor,
        global_cell_embedding: np.ndarray,
        num_layers: int = 2,
        edge_dropout: float = 0.1,
        temperature: float = 0.2,
        num_modules: int = 16,
        module_top_k: int = 30,
    ):
        super().__init__()
        self.n_cells = int(n_cells)
        self.n_genes = int(n_genes)
        self.latent_dim = int(latent_dim)
        self.num_layers = int(num_layers)
        self.edge_dropout = float(edge_dropout)
        self.temperature = float(temperature)

        self.cell_embedding = nn.Embedding(self.n_cells, self.latent_dim)
        self.gene_embedding = nn.Embedding(self.n_genes, self.latent_dim)
        nn.init.xavier_uniform_(self.cell_embedding.weight)
        nn.init.xavier_uniform_(self.gene_embedding.weight)

        global_tensor = torch.as_tensor(global_cell_embedding, dtype=torch.float32)
        self.register_buffer("global_cell_embedding", global_tensor, persistent=False)
        self.global_projection = nn.Linear(global_tensor.shape[1], self.latent_dim, bias=False)
        self.register_buffer("adj_norm", adj_norm.coalesce(), persistent=False)

        self.module_layer: Optional[TopKSparseModuleLayer]
        if num_modules and num_modules > 0:
            self.module_layer = TopKSparseModuleLayer(self.latent_dim, num_modules, module_top_k)
        else:
            self.module_layer = None

    def propagate(self, edge_dropout: Optional[float] = None) -> Tuple[torch.Tensor, torch.Tensor]:
        adj = self.adj_norm
        if self.training:
            adj = sparse_dropout(adj, self.edge_dropout if edge_dropout is None else edge_dropout)
        adj_t = adj.transpose(0, 1).coalesce()

        cell_emb = self.cell_embedding.weight
        gene_emb = self.gene_embedding.weight
        cell_layers = [cell_emb]
        gene_layers = [gene_emb]
        for _ in range(self.num_layers):
            next_cell = torch.sparse.mm(adj, gene_layers[-1])
            next_gene = torch.sparse.mm(adj_t, cell_layers[-1])
            cell_layers.append(next_cell)
            gene_layers.append(next_gene)
        return torch.stack(cell_layers, dim=0).mean(dim=0), torch.stack(gene_layers, dim=0).mean(dim=0)

    def projected_global(self) -> torch.Tensor:
        return self.global_projection(self.global_cell_embedding)

    def forward(
        self,
        cells: torch.Tensor,
        pos_genes: torch.Tensor,
        neg_genes: torch.Tensor,
        contrastive_cells: Optional[torch.Tensor] = None,
        module_cells: Optional[torch.Tensor] = None,
    ) -> Dict[str, torch.Tensor]:
        local_cells, local_genes = self.propagate()
        user_emb = local_cells[cells]
        pos_emb = local_genes[pos_genes]
        neg_emb = local_genes[neg_genes]
        pos_scores = torch.sum(user_emb * pos_emb, dim=1)
        neg_scores = torch.sum(user_emb * neg_emb, dim=1)
        bpr = F.softplus(-(pos_scores - neg_scores)).mean()

        contrastive = torch.zeros((), dtype=torch.float32, device=cells.device)
        if contrastive_cells is not None and contrastive_cells.numel() > 1:
            local = F.normalize(local_cells[contrastive_cells], dim=1)
            global_view = F.normalize(self.projected_global()[contrastive_cells], dim=1)
            logits = local @ global_view.T / self.temperature
            labels = torch.arange(logits.shape[0], device=logits.device)
            contrastive = 0.5 * (F.cross_entropy(logits, labels) + F.cross_entropy(logits.T, labels))

        module_loss = torch.zeros((), dtype=torch.float32, device=cells.device)
        if self.module_layer is not None and module_cells is not None and module_cells.numel() > 0:
            module_out = self.module_layer(F.normalize(local_cells, dim=1), F.normalize(local_genes, dim=1))
            weights = torch.softmax(module_out["activations"][module_cells], dim=1)
            recon = weights @ module_out["module_emb"]
            module_loss = (1.0 - F.cosine_similarity(recon, F.normalize(local_cells[module_cells], dim=1), dim=1)).mean()

        return {
            "bpr": bpr,
            "contrastive": contrastive,
            "module": module_loss,
            "local_cells": local_cells,
            "local_genes": local_genes,
        }

    @torch.no_grad()
    def get_embeddings(self, normalize_output: bool = True) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        was_training = self.training
        self.eval()
        local_cells, local_genes = self.propagate(edge_dropout=0.0)
        global_cells = self.projected_global()
        if normalize_output:
            local_cells = F.normalize(local_cells, dim=1)
            local_genes = F.normalize(local_genes, dim=1)
            global_cells = F.normalize(global_cells, dim=1)
        if was_training:
            self.train()
        return (
            local_cells.detach().cpu().numpy().astype(np.float32),
            local_genes.detach().cpu().numpy().astype(np.float32),
            global_cells.detach().cpu().numpy().astype(np.float32),
        )

    @torch.no_grad()
    def module_top_genes(self, gene_names: np.ndarray, gene_embedding: Optional[np.ndarray] = None, top_k: Optional[int] = None):
        if self.module_layer is None:
            return []
        was_training = self.training
        self.eval()
        if gene_embedding is None:
            _, gene_emb, _ = self.get_embeddings(normalize_output=True)
            gene_tensor = torch.as_tensor(gene_emb, dtype=torch.float32, device=self.module_layer.module_queries.device)
        else:
            gene_tensor = torch.as_tensor(gene_embedding, dtype=torch.float32, device=self.module_layer.module_queries.device)
        out = self.module_layer(torch.empty((0, self.latent_dim), device=gene_tensor.device), F.normalize(gene_tensor, dim=1))
        idx = out["top_idx"].detach().cpu().numpy()
        scores = out["top_scores"].detach().cpu().numpy()
        weights = out["top_weights"].detach().cpu().numpy()
        k = idx.shape[1] if top_k is None else min(int(top_k), idx.shape[1])
        rows = []
        for module_id in range(idx.shape[0]):
            for rank in range(k):
                gene_id = int(idx[module_id, rank])
                rows.append(
                    {
                        "module": int(module_id),
                        "rank": int(rank + 1),
                        "gene": str(gene_names[gene_id]),
                        "gene_index": gene_id,
                        "score": float(scores[module_id, rank]),
                        "weight": float(weights[module_id, rank]),
                    }
                )
        if was_training:
            self.train()
        return rows


@dataclass
class LGCLTrainConfig:
    epochs: int = 80
    pairs_per_epoch: int = 262144
    contrastive_batch_size: int = 2048
    lr: float = 1e-3
    weight_decay: float = 1e-5
    contrastive_weight: float = 0.05
    module_weight: float = 0.001
    grad_clip: float = 5.0
    log_interval: int = 5
    seed: int = 42


def train_lgcl(
    model: PlantSPADELGCL,
    support: sp.csr_matrix,
    config: LGCLTrainConfig,
    device: torch.device,
) -> Dict[str, list]:
    sampler = PairSampler(support, seed=config.seed)
    optimizer = torch.optim.AdamW(model.parameters(), lr=config.lr, weight_decay=config.weight_decay)
    scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=max(1, config.epochs))
    rng = np.random.default_rng(config.seed + 1009)
    history = {"loss": [], "bpr": [], "contrastive": [], "module": []}

    for epoch in range(1, config.epochs + 1):
        model.train()
        cells_np, pos_np, neg_np = sampler.sample(config.pairs_per_epoch)
        cells = torch.as_tensor(cells_np, dtype=torch.long, device=device)
        pos = torch.as_tensor(pos_np, dtype=torch.long, device=device)
        neg = torch.as_tensor(neg_np, dtype=torch.long, device=device)

        csz = min(config.contrastive_batch_size, model.n_cells)
        contrastive_cells = torch.as_tensor(rng.choice(model.n_cells, size=csz, replace=False), dtype=torch.long, device=device)
        module_cells = contrastive_cells if model.module_layer is not None and config.module_weight > 0.0 else None

        optimizer.zero_grad(set_to_none=True)
        out = model(cells, pos, neg, contrastive_cells=contrastive_cells, module_cells=module_cells)
        loss = out["bpr"] + config.contrastive_weight * out["contrastive"] + config.module_weight * out["module"]
        loss.backward()
        if config.grad_clip and config.grad_clip > 0:
            torch.nn.utils.clip_grad_norm_(model.parameters(), config.grad_clip)
        optimizer.step()
        scheduler.step()

        metrics = {
            "loss": float(loss.detach().cpu()),
            "bpr": float(out["bpr"].detach().cpu()),
            "contrastive": float(out["contrastive"].detach().cpu()),
            "module": float(out["module"].detach().cpu()),
        }
        for key, value in metrics.items():
            history[key].append(value)
        if epoch == 1 or epoch == config.epochs or epoch % max(1, config.log_interval) == 0:
            print(
                f"Epoch {epoch:03d}/{config.epochs} "
                f"loss={metrics['loss']:.4f} bpr={metrics['bpr']:.4f} "
                f"cl={metrics['contrastive']:.4f} module={metrics['module']:.4f}"
            )
    return history
