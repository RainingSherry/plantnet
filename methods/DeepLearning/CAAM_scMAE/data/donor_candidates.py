from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import torch


class DonorCandidateProvider:
    """Label-free matched donor candidate provider.

    The provider is intentionally label-free. It may use batch code, library-size bins,
    and zero-ratio bins, but it must never use true cell-type labels or cluster labels.
    """

    def __init__(
        self,
        x: np.ndarray,
        batch_code: np.ndarray,
        library_size: np.ndarray,
        zero_ratio: np.ndarray,
        *,
        candidate_pool_size: int,
        library_size_bins: int,
        zero_ratio_bins: int,
        atol: float,
        rtol: float,
        seed: int,
    ) -> None:
        self.x_np = np.asarray(x, dtype=np.float32)
        self.n_cells, self.n_genes = self.x_np.shape
        if self.n_cells <= 1:
            raise ValueError(
                "CAAM matched donor corruption requires at least two cells; "
                "self-donor fallback is forbidden by the BDD."
            )
        self.atol = float(atol)
        self.rtol = float(rtol)
        self.rng = np.random.default_rng(seed)
        self.batch_code = np.asarray(batch_code, dtype=np.int64)
        self.library_bin = self._bin(library_size, library_size_bins)
        self.zero_bin = self._bin(zero_ratio, zero_ratio_bins)
        self.candidate_pool_size = int(candidate_pool_size)
        if self.candidate_pool_size <= 0:
            raise ValueError("candidate_pool_size must be positive.")
        self.candidates, self.stats = self._build_candidates()

    @staticmethod
    def _bin(values: np.ndarray, n_bins: int) -> np.ndarray:
        values = np.asarray(values, dtype=np.float32)
        if n_bins <= 1 or np.allclose(values, values[0]):
            return np.zeros(values.shape[0], dtype=np.int64)
        quantiles = np.quantile(values, np.linspace(0.0, 1.0, int(n_bins) + 1)[1:-1])
        return np.searchsorted(quantiles, values, side="right").astype(np.int64)

    def _sample_from_pool(self, pool: np.ndarray, size: int) -> np.ndarray:
        if pool.size == 0:
            raise ValueError("Cannot sample from empty donor pool.")
        replace = pool.size < size
        return self.rng.choice(pool, size=size, replace=replace).astype(np.int64)

    def _pool_for(self, i: int) -> tuple[np.ndarray, str]:
        all_idx = np.arange(self.n_cells, dtype=np.int64)
        not_self = all_idx[all_idx != i]
        same_batch = not_self[self.batch_code[not_self] == self.batch_code[i]]
        matched = same_batch[
            (self.library_bin[same_batch] == self.library_bin[i])
            & (self.zero_bin[same_batch] == self.zero_bin[i])
        ]
        if matched.size:
            return matched, "matched"
        nearest = same_batch[
            (np.abs(self.library_bin[same_batch] - self.library_bin[i]) <= 1)
            & (np.abs(self.zero_bin[same_batch] - self.zero_bin[i]) <= 1)
        ]
        if nearest.size:
            return nearest, "batch"
        global_same = not_self[
            (self.library_bin[not_self] == self.library_bin[i])
            & (self.zero_bin[not_self] == self.zero_bin[i])
        ]
        if global_same.size:
            return global_same, "global"
        return not_self, "global"

    def _build_candidates(self) -> tuple[np.ndarray, dict]:
        candidates = np.zeros((self.n_cells, self.candidate_pool_size), dtype=np.int64)
        levels: dict[str, int] = {"matched": 0, "batch": 0, "global": 0}
        for i in range(self.n_cells):
            pool, level = self._pool_for(i)
            if np.any(pool == i):
                raise RuntimeError("Internal donor-pool error: self index appears in donor pool.")
            levels[level] = levels.get(level, 0) + 1
            candidates[i] = self._sample_from_pool(pool, self.candidate_pool_size)
        if np.any(candidates == np.arange(self.n_cells, dtype=np.int64)[:, None]):
            raise RuntimeError("Internal donor-candidate error: self donor was sampled.")
        return candidates, {"fallback_levels": levels, "single_cell": False, "self_donor_forbidden": True}

    def save(self, save_dir: Path) -> None:
        save_dir.mkdir(parents=True, exist_ok=True)
        np.save(save_dir / "donor_candidate_indices.npy", self.candidates)
        with open(save_dir / "donor_candidate_stats.json", "w", encoding="utf-8") as handle:
            json.dump(self.stats, handle, indent=2)

    def sample_batch(
        self,
        batch_indices: torch.Tensor,
        x_full: torch.Tensor,
        device: torch.device,
    ) -> dict[str, torch.Tensor]:
        idx_cpu = batch_indices.detach().cpu().numpy().astype(np.int64)
        b = int(idx_cpu.shape[0])
        g = self.n_genes
        slot = torch.randint(
            low=0,
            high=self.candidate_pool_size,
            size=(b, g),
            device=device,
        )
        cand = torch.as_tensor(self.candidates[idx_cpu], dtype=torch.long, device=device)
        donor_indices = torch.gather(cand, dim=1, index=slot)
        if torch.any(donor_indices == batch_indices.to(device).long().view(-1, 1)):
            raise RuntimeError("Matched donor corruption selected self donor, which is forbidden.")
        gene_ids = torch.arange(g, device=device).view(1, g).expand(b, g)
        replacement = x_full[donor_indices, gene_ids]
        original = x_full[batch_indices.to(device).long()]
        eligibility = ~torch.isclose(replacement, original, atol=self.atol, rtol=self.rtol)
        return {
            "replacement": replacement,
            "eligibility": eligibility,
            "donor_indices": donor_indices,
        }
