from dataclasses import dataclass
from typing import Optional, Tuple

import numpy as np
import scipy.sparse as sp
import torch
import torch.nn as nn
import torch.nn.functional as F


@dataclass
class SparseAttentionWeights:
    cells: np.ndarray
    genes: np.ndarray
    weights: np.ndarray
    amplitudes: np.ndarray


class SupportGeneAttention(nn.Module):
    """Sparse gene-set cross-attention over each cell's observed support only.

    This is intentionally not a full Transformer block: every cell attends only to
    genes in its non-zero support set S_c, optionally truncated by expression
    amplitude. No dense [n_cells, n_genes] attention matrix is materialized.
    """

    def __init__(
        self,
        support,
        amplitude: sp.csr_matrix,
        gene_idf: torch.Tensor,
        top_k_genes: int = 128,
        beta: float = 0.1,
        gamma: float = 0.1,
        eta: float = 0.5,
        dropout: float = 0.1,
        trainable: bool = False,
    ):
        super().__init__()
        self.support = self._as_csr_support(support)
        self.amplitude = self._as_csr(amplitude)
        if self.support.shape != self.amplitude.shape:
            raise ValueError(f"support shape {self.support.shape} != amplitude shape {self.amplitude.shape}")
        self.top_k_genes = int(top_k_genes) if top_k_genes and top_k_genes > 0 else 0
        self.beta = float(beta)
        self.gamma = float(gamma)
        self.eta = float(eta)
        self.dropout = float(dropout)
        self.trainable = bool(trainable)
        if self.trainable:
            self.beta_param = nn.Parameter(torch.tensor(float(beta), dtype=torch.float32))
            self.gamma_param = nn.Parameter(torch.tensor(float(gamma), dtype=torch.float32))
            self.eta_param = nn.Parameter(torch.tensor(float(eta), dtype=torch.float32))
        else:
            self.register_parameter("beta_param", None)
            self.register_parameter("gamma_param", None)
            self.register_parameter("eta_param", None)
        self.register_buffer("gene_idf", gene_idf.detach().float().clone(), persistent=False)

        indptr, indices, values = self._build_truncated_edges()
        self.selected_indptr = indptr
        self.selected_indices = indices
        self.selected_amplitudes = values

    @staticmethod
    def _as_csr(matrix) -> sp.csr_matrix:
        if sp.issparse(matrix):
            out = matrix.tocsr().astype(np.float32)
        else:
            out = sp.csr_matrix(np.asarray(matrix, dtype=np.float32))
        out.data = np.nan_to_num(out.data, nan=0.0, posinf=0.0, neginf=0.0)
        out.data[out.data < 0.0] = 0.0
        out.eliminate_zeros()
        out.sort_indices()
        return out

    @classmethod
    def _as_csr_support(cls, support) -> sp.csr_matrix:
        if isinstance(support, tuple) and len(support) == 3:
            rows, cols, shape = support
            data = np.ones(len(rows), dtype=np.float32)
            support = sp.coo_matrix((data, (rows, cols)), shape=shape)
        out = cls._as_csr(support)
        out.data = np.ones_like(out.data, dtype=np.float32)
        out.eliminate_zeros()
        out.sort_indices()
        return out

    def _row_amplitudes(self, cell: int, genes: np.ndarray) -> np.ndarray:
        start = self.amplitude.indptr[cell]
        end = self.amplitude.indptr[cell + 1]
        amp_genes = self.amplitude.indices[start:end]
        amp_vals = self.amplitude.data[start:end]
        if amp_genes.size == 0:
            return np.zeros(genes.size, dtype=np.float32)
        if amp_genes.size == genes.size and np.array_equal(amp_genes, genes):
            return amp_vals.astype(np.float32, copy=False)
        positions = np.searchsorted(amp_genes, genes)
        values = np.zeros(genes.size, dtype=np.float32)
        valid = positions < amp_genes.size
        valid[valid] = amp_genes[positions[valid]] == genes[valid]
        if np.any(valid):
            values[valid] = amp_vals[positions[valid]]
        return values

    def _build_truncated_edges(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        n_cells = self.support.shape[0]
        indptr = np.zeros(n_cells + 1, dtype=np.int64)
        all_indices = []
        all_amplitudes = []
        for cell in range(n_cells):
            start = self.support.indptr[cell]
            end = self.support.indptr[cell + 1]
            genes = self.support.indices[start:end]
            if genes.size == 0:
                indptr[cell + 1] = indptr[cell]
                continue
            amplitudes = self._row_amplitudes(cell, genes)
            if self.top_k_genes > 0 and genes.size > self.top_k_genes:
                keep = np.argpartition(-amplitudes, self.top_k_genes - 1)[: self.top_k_genes]
                order = np.argsort(-amplitudes[keep], kind="stable")
                keep = keep[order]
                genes = genes[keep]
                amplitudes = amplitudes[keep]
            all_indices.append(genes.astype(np.int64, copy=False))
            all_amplitudes.append(amplitudes.astype(np.float32, copy=False))
            indptr[cell + 1] = indptr[cell] + genes.size
        if all_indices:
            indices = np.concatenate(all_indices).astype(np.int64, copy=False)
            amplitudes = np.concatenate(all_amplitudes).astype(np.float32, copy=False)
        else:
            indices = np.empty(0, dtype=np.int64)
            amplitudes = np.empty(0, dtype=np.float32)
        return indptr, indices, amplitudes

    def _coef(self, name: str, device: torch.device) -> torch.Tensor:
        param = getattr(self, f"{name}_param")
        if param is not None:
            return param.to(device=device)
        return torch.tensor(getattr(self, name), dtype=torch.float32, device=device)

    def _batch_edges(self, batch_cells: Optional[torch.Tensor], device: torch.device):
        n_cells = self.support.shape[0]
        if batch_cells is None:
            cells_np = np.arange(n_cells, dtype=np.int64)
        elif torch.is_tensor(batch_cells):
            cells_np = batch_cells.detach().cpu().numpy().astype(np.int64, copy=False)
        else:
            cells_np = np.asarray(batch_cells, dtype=np.int64)

        edge_positions = []
        edge_genes = []
        edge_amplitudes = []
        edge_cells = []
        for local_pos, cell in enumerate(cells_np):
            start = self.selected_indptr[cell]
            end = self.selected_indptr[cell + 1]
            if end <= start:
                continue
            count = end - start
            edge_positions.append(np.full(count, local_pos, dtype=np.int64))
            edge_cells.append(np.full(count, cell, dtype=np.int64))
            edge_genes.append(self.selected_indices[start:end])
            edge_amplitudes.append(self.selected_amplitudes[start:end])

        if edge_genes:
            local_pos_np = np.concatenate(edge_positions)
            cells_edge_np = np.concatenate(edge_cells)
            genes_np = np.concatenate(edge_genes)
            amps_np = np.concatenate(edge_amplitudes).astype(np.float32, copy=False)
        else:
            local_pos_np = np.empty(0, dtype=np.int64)
            cells_edge_np = np.empty(0, dtype=np.int64)
            genes_np = np.empty(0, dtype=np.int64)
            amps_np = np.empty(0, dtype=np.float32)

        return (
            torch.as_tensor(cells_np, dtype=torch.long, device=device),
            torch.as_tensor(local_pos_np, dtype=torch.long, device=device),
            torch.as_tensor(genes_np, dtype=torch.long, device=device),
            torch.as_tensor(amps_np, dtype=torch.float32, device=device),
            cells_edge_np,
            genes_np,
            amps_np,
        )

    def forward(
        self,
        cell_emb: torch.Tensor,
        gene_emb: torch.Tensor,
        batch_cells: Optional[torch.Tensor] = None,
        return_attention: bool = False,
    ):
        device = cell_emb.device
        batch_cells_t, edge_pos, edge_genes, edge_amp, edge_cells_np, edge_genes_np, edge_amp_np = self._batch_edges(
            batch_cells,
            device,
        )
        base = cell_emb[batch_cells_t]
        if edge_genes.numel() == 0:
            if return_attention:
                empty = SparseAttentionWeights(edge_cells_np, edge_genes_np, np.empty(0, dtype=np.float32), edge_amp_np)
                return base, empty
            return base

        dim = cell_emb.shape[1]
        query = base[edge_pos]
        key_value = gene_emb[edge_genes]
        logits = (query * key_value).sum(dim=1) / np.sqrt(float(dim))
        beta = self._coef("beta", device)
        gamma = self._coef("gamma", device)
        eta = self._coef("eta", device)
        if self.trainable or self.beta != 0.0:
            logits = logits + beta * torch.log1p(edge_amp)
        if self.trainable or self.gamma != 0.0:
            logits = logits + gamma * self.gene_idf.to(device=device)[edge_genes]

        n_batch = base.shape[0]
        max_logits = torch.full((n_batch,), -torch.inf, dtype=logits.dtype, device=device)
        max_logits.scatter_reduce_(0, edge_pos, logits, reduce="amax", include_self=True)
        exp_logits = torch.exp(logits - max_logits[edge_pos])
        denom = torch.zeros((n_batch,), dtype=logits.dtype, device=device)
        denom.scatter_add_(0, edge_pos, exp_logits)
        alpha = exp_logits / denom[edge_pos].clamp_min(1e-12)
        if self.training and self.dropout > 0.0:
            alpha = F.dropout(alpha, p=self.dropout, training=True)

        attended = torch.zeros_like(base)
        attended.index_add_(0, edge_pos, alpha.unsqueeze(1) * key_value)
        refined = base + eta * attended

        if return_attention:
            weights = SparseAttentionWeights(
                cells=edge_cells_np,
                genes=edge_genes_np,
                weights=alpha.detach().cpu().numpy().astype(np.float32),
                amplitudes=edge_amp_np.astype(np.float32, copy=False),
            )
            return refined, weights
        return refined


class TrainableSupportGeneAttentionRefiner(nn.Module):
    """Lightweight trainable sparse refiner over observed support sets.

    Only beta, gamma, and eta are optimized by default. Cell and gene
    embeddings are supplied by PlantSPADE-LGCL and remain external tensors.
    """

    def __init__(
        self,
        support,
        amplitude: sp.csr_matrix,
        gene_idf: torch.Tensor,
        top_k_genes: int = 128,
        beta: float = 0.1,
        gamma: float = 0.1,
        eta: float = 0.5,
        dropout: float = 0.1,
    ):
        super().__init__()
        self.attention = SupportGeneAttention(
            support=support,
            amplitude=amplitude,
            gene_idf=gene_idf,
            top_k_genes=top_k_genes,
            beta=beta,
            gamma=gamma,
            eta=eta,
            dropout=dropout,
            trainable=True,
        )

    def forward(
        self,
        cell_emb: torch.Tensor,
        gene_emb: torch.Tensor,
        batch_cells: Optional[torch.Tensor] = None,
        return_attention: bool = False,
    ):
        return self.attention(cell_emb, gene_emb, batch_cells=batch_cells, return_attention=return_attention)
