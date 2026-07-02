from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Tuple

import numpy as np
import scipy.sparse as sp
import torch
import torch.nn as nn
import torch.nn.functional as F

from ..support_gene_attention import SparseAttentionWeights


@dataclass
class LearnedTopKConfig:
    """Config for learned top-k routing over each cell's support set."""

    top_k_genes: int = 128
    # contribution weights (same semantics as SupportGeneAttention)
    beta_amplitude: float = 0.1
    gamma_idf: float = 0.1
    # similarity scaling: logits = sim / sqrt(d) * sim_scale
    sim_scale: float = 1.0
    # residual strength: refined = base + eta * attended
    eta: float = 0.5
    dropout: float = 0.1


class LearnedTopKSupportAttention(nn.Module):
    """Sparse gene-set attention with *learned* top-k selection.

    This differs from `SupportGeneAttention` in one key aspect:
    - top-k truncation is based on a *routing score* that uses
      cell-gene similarity + amplitude + IDF

    Motivation (from the CV "routing / sparse attention" family):
    - amplitude-only top-k tends to over-select housekeeping / high-expression genes
    - routing prefers genes that are both present *and* discriminative for the cell

    Computational constraints:
    - still only attends within each cell's observed support set S_c
    - no dense [n_cells, n_genes] matrix is materialized
    """

    def __init__(
        self,
        support,
        amplitude: sp.csr_matrix,
        gene_idf: torch.Tensor,
        config: LearnedTopKConfig,
    ):
        super().__init__()
        self.support = self._as_csr_support(support)
        self.amplitude = self._as_csr(amplitude)
        if self.support.shape != self.amplitude.shape:
            raise ValueError(f"support shape {self.support.shape} != amplitude shape {self.amplitude.shape}")

        self.config = config
        self.register_buffer("gene_idf", gene_idf.detach().float().clone(), persistent=False)

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

    @torch.no_grad()
    def _select_topk_edges(
        self,
        cell_emb: torch.Tensor,
        gene_emb: torch.Tensor,
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Compute a routing top-k per cell (CPU numpy indices).

        Returns CSR-like arrays (indptr, indices, amplitudes) over selected genes.
        """
        if cell_emb.ndim != 2 or gene_emb.ndim != 2:
            raise ValueError("Expected 2-D embeddings")
        if cell_emb.shape[1] != gene_emb.shape[1]:
            raise ValueError(f"Dim mismatch: cell_emb={cell_emb.shape}, gene_emb={gene_emb.shape}")

        n_cells = self.support.shape[0]
        top_k = int(self.config.top_k_genes) if self.config.top_k_genes and self.config.top_k_genes > 0 else 0

        # move to CPU for per-cell variable-length routing
        cell_cpu = cell_emb.detach().cpu()
        gene_cpu = gene_emb.detach().cpu()
        idf_cpu = self.gene_idf.detach().cpu()

        dim = int(cell_cpu.shape[1])
        inv_sqrt_d = 1.0 / np.sqrt(float(dim))

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

            amps = self._row_amplitudes(cell, genes)

            # routing score: similarity + beta * log1p(amp) + gamma * idf
            q = cell_cpu[cell]  # (d,)
            k = gene_cpu[genes]  # (deg, d)
            sim = (k * q.unsqueeze(0)).sum(dim=1).numpy() * float(inv_sqrt_d) * float(self.config.sim_scale)
            score = sim
            if self.config.beta_amplitude != 0.0:
                score = score + float(self.config.beta_amplitude) * np.log1p(amps)
            if self.config.gamma_idf != 0.0:
                score = score + float(self.config.gamma_idf) * idf_cpu[genes].numpy()

            if top_k > 0 and genes.size > top_k:
                keep = np.argpartition(-score, top_k - 1)[:top_k]
                order = np.argsort(-score[keep], kind="stable")
                keep = keep[order]
                genes = genes[keep]
                amps = amps[keep]

            all_indices.append(genes.astype(np.int64, copy=False))
            all_amplitudes.append(amps.astype(np.float32, copy=False))
            indptr[cell + 1] = indptr[cell] + genes.size

        if all_indices:
            indices = np.concatenate(all_indices).astype(np.int64, copy=False)
            amplitudes = np.concatenate(all_amplitudes).astype(np.float32, copy=False)
        else:
            indices = np.empty(0, dtype=np.int64)
            amplitudes = np.empty(0, dtype=np.float32)

        return indptr, indices, amplitudes

    def forward(
        self,
        cell_emb: torch.Tensor,
        gene_emb: torch.Tensor,
        batch_cells: Optional[torch.Tensor] = None,
        return_attention: bool = False,
    ):
        # Build selected sparse edges once per call.
        # Note: for now we recompute routing for the *full* dataset; this keeps code simple.
        # If needed we can optimize to batch routing later.
        indptr, indices, amps = self._select_topk_edges(cell_emb, gene_emb)

        device = cell_emb.device
        n_cells = self.support.shape[0]
        if batch_cells is None:
            cells_np = np.arange(n_cells, dtype=np.int64)
        else:
            cells_np = batch_cells.detach().cpu().numpy().astype(np.int64, copy=False)

        edge_pos = []
        edge_cells = []
        edge_genes = []
        edge_amp = []
        for local_pos, cell in enumerate(cells_np):
            s = indptr[cell]
            e = indptr[cell + 1]
            if e <= s:
                continue
            count = e - s
            edge_pos.append(np.full(count, local_pos, dtype=np.int64))
            edge_cells.append(np.full(count, cell, dtype=np.int64))
            edge_genes.append(indices[s:e])
            edge_amp.append(amps[s:e])

        base = cell_emb[torch.as_tensor(cells_np, dtype=torch.long, device=device)]
        if not edge_genes:
            if return_attention:
                empty = SparseAttentionWeights(
                    cells=np.empty(0, dtype=np.int64),
                    genes=np.empty(0, dtype=np.int64),
                    weights=np.empty(0, dtype=np.float32),
                    amplitudes=np.empty(0, dtype=np.float32),
                )
                return base, empty
            return base

        edge_pos_t = torch.as_tensor(np.concatenate(edge_pos), dtype=torch.long, device=device)
        edge_cells_np = np.concatenate(edge_cells)
        edge_genes_np = np.concatenate(edge_genes)
        edge_amp_np = np.concatenate(edge_amp).astype(np.float32, copy=False)

        edge_genes_t = torch.as_tensor(edge_genes_np, dtype=torch.long, device=device)
        edge_amp_t = torch.as_tensor(edge_amp_np, dtype=torch.float32, device=device)

        query = base[edge_pos_t]
        key_value = gene_emb[edge_genes_t]
        dim = cell_emb.shape[1]
        logits = (query * key_value).sum(dim=1) / np.sqrt(float(dim))

        # NOTE: truncation already used amp/idf/sim; here we keep logits as pure sim
        n_batch = base.shape[0]
        max_logits = torch.full((n_batch,), -torch.inf, dtype=logits.dtype, device=device)
        max_logits.scatter_reduce_(0, edge_pos_t, logits, reduce="amax", include_self=True)
        exp_logits = torch.exp(logits - max_logits[edge_pos_t])
        denom = torch.zeros((n_batch,), dtype=logits.dtype, device=device)
        denom.scatter_add_(0, edge_pos_t, exp_logits)
        alpha = exp_logits / denom[edge_pos_t].clamp_min(1e-12)
        if self.training and self.config.dropout > 0.0:
            alpha = F.dropout(alpha, p=float(self.config.dropout), training=True)

        attended = torch.zeros_like(base)
        attended.index_add_(0, edge_pos_t, alpha.unsqueeze(1) * key_value)
        refined = base + float(self.config.eta) * attended

        if return_attention:
            weights = SparseAttentionWeights(
                cells=edge_cells_np,
                genes=edge_genes_np,
                weights=alpha.detach().cpu().numpy().astype(np.float32),
                amplitudes=edge_amp_np,
            )
            return refined, weights
        return refined
