from __future__ import annotations

from dataclasses import dataclass
from typing import List, Optional

import numpy as np
import scipy.sparse as sp
from sklearn.decomposition import TruncatedSVD
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import normalize


@dataclass
class NegativeSamplerConfig:
    mode: str = "random_zero"
    seed: int = 42
    neighbor_k: int = 15
    max_conflict_candidates: int = 512


class NegativeSampler:
    """Sample zero-expression genes for BPR without materializing dense zeros."""

    VALID_MODES = {"random_zero", "idf_weighted_zero", "neighbor_conflict_zero"}

    def __init__(
        self,
        support: sp.csr_matrix,
        config: Optional[NegativeSamplerConfig] = None,
    ):
        self.config = config or NegativeSamplerConfig()
        if self.config.mode not in self.VALID_MODES:
            raise ValueError(f"Unknown negative sampler {self.config.mode!r}; expected one of {sorted(self.VALID_MODES)}")

        self.support = support.astype(np.float32).tocsr(copy=True)
        self.support.data = np.ones_like(self.support.data, dtype=np.float32)
        self.support.eliminate_zeros()
        self.support.sort_indices()
        self.indptr = self.support.indptr.astype(np.int64, copy=False)
        self.indices = self.support.indices.astype(np.int64, copy=False)
        self.n_cells, self.n_genes = self.support.shape
        self.rng = np.random.default_rng(self.config.seed)

        df = np.diff(self.support.tocsc().indptr).astype(np.float64)
        # This mode intentionally emphasizes globally common genes that are
        # absent in the current cell. They are high-information zeros despite
        # having low classical IDF.
        self.global_weights = np.log1p(df)
        if not np.any(self.global_weights > 0):
            self.global_weights = np.ones(self.n_genes, dtype=np.float64)
        self.global_weights = self.global_weights / self.global_weights.sum()

        self.conflict_candidates: Optional[List[np.ndarray]] = None
        self.conflict_weights: Optional[List[np.ndarray]] = None
        if self.config.mode == "neighbor_conflict_zero":
            self._build_neighbor_conflicts()

    def _row_contains(self, row: int, col: int) -> bool:
        start = self.indptr[row]
        end = self.indptr[row + 1]
        row_indices = self.indices[start:end]
        pos = np.searchsorted(row_indices, col)
        return bool(pos < row_indices.size and row_indices[pos] == col)

    def _sample_random_zero_one(self, row: int) -> int:
        if self.indptr[row + 1] - self.indptr[row] >= self.n_genes:
            return int(self.rng.integers(0, self.n_genes))
        neg = int(self.rng.integers(0, self.n_genes))
        while self._row_contains(row, neg):
            neg = int(self.rng.integers(0, self.n_genes))
        return neg

    def _sample_weighted_zero_one(self, row: int) -> int:
        if self.indptr[row + 1] - self.indptr[row] >= self.n_genes:
            return int(self.rng.integers(0, self.n_genes))
        for _ in range(32):
            neg = int(self.rng.choice(self.n_genes, p=self.global_weights))
            if not self._row_contains(row, neg):
                return neg
        return self._sample_random_zero_one(row)

    def _build_neighbor_conflicts(self) -> None:
        n_components = max(2, min(32, self.n_cells - 1, self.n_genes - 1))
        if n_components < 2 or self.n_cells < 3:
            self.conflict_candidates = [np.empty(0, dtype=np.int64) for _ in range(self.n_cells)]
            self.conflict_weights = [np.empty(0, dtype=np.float64) for _ in range(self.n_cells)]
            return

        svd = TruncatedSVD(n_components=n_components, random_state=self.config.seed)
        emb = svd.fit_transform(self.support).astype(np.float32)
        emb = normalize(emb, norm="l2", axis=1, copy=False)
        n_neighbors = min(max(2, self.config.neighbor_k + 1), self.n_cells)
        nn = NearestNeighbors(n_neighbors=n_neighbors, metric="cosine")
        nn.fit(emb)
        _, neighbor_idx = nn.kneighbors(emb, return_distance=True)

        csc = self.support.tocsc()
        self.conflict_candidates = []
        self.conflict_weights = []
        for cell in range(self.n_cells):
            neighbors = neighbor_idx[cell, 1:]
            if neighbors.size == 0:
                self.conflict_candidates.append(np.empty(0, dtype=np.int64))
                self.conflict_weights.append(np.empty(0, dtype=np.float64))
                continue
            gene_counts = np.asarray(self.support[neighbors].sum(axis=0)).ravel()
            start, end = self.indptr[cell], self.indptr[cell + 1]
            gene_counts[self.indices[start:end]] = 0.0
            nz = np.flatnonzero(gene_counts > 0)
            if nz.size > self.config.max_conflict_candidates:
                keep = np.argpartition(-gene_counts[nz], self.config.max_conflict_candidates - 1)[
                    : self.config.max_conflict_candidates
                ]
                nz = nz[keep]
            weights = gene_counts[nz].astype(np.float64)
            # A light global-frequency tie-break keeps ubiquitous housekeeping
            # genes from dominating when neighbor counts are tied.
            df = np.diff(csc.indptr).astype(np.float64)
            weights = weights * np.log1p(df[nz]).clip(min=1.0)
            weights = weights / weights.sum() if np.any(weights > 0) else np.ones(nz.size) / max(1, nz.size)
            self.conflict_candidates.append(nz.astype(np.int64, copy=False))
            self.conflict_weights.append(weights)

    def _sample_neighbor_conflict_one(self, row: int) -> int:
        assert self.conflict_candidates is not None
        assert self.conflict_weights is not None
        candidates = self.conflict_candidates[row]
        if candidates.size == 0:
            return self._sample_weighted_zero_one(row)
        weights = self.conflict_weights[row]
        return int(self.rng.choice(candidates, p=weights))

    def sample(self, cells: np.ndarray) -> np.ndarray:
        cells = np.asarray(cells, dtype=np.int64)
        out = np.empty(cells.shape[0], dtype=np.int64)
        if self.config.mode == "random_zero":
            for i, row in enumerate(cells):
                out[i] = self._sample_random_zero_one(int(row))
        elif self.config.mode == "idf_weighted_zero":
            for i, row in enumerate(cells):
                out[i] = self._sample_weighted_zero_one(int(row))
        elif self.config.mode == "neighbor_conflict_zero":
            for i, row in enumerate(cells):
                out[i] = self._sample_neighbor_conflict_one(int(row))
        return out
