from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from sklearn.cluster import MiniBatchKMeans
from sklearn.decomposition import PCA
from sklearn.metrics import pairwise_distances_argmin_min


def select_context_indices(
    x: np.ndarray,
    context_size: int,
    context_pca_dim: int,
    context_seed: int,
    save_dir: Path | None = None,
) -> np.ndarray:
    n = int(x.shape[0])
    size = max(1, min(int(context_size), n))
    if n <= size:
        indices = np.arange(n, dtype=np.int64)
    else:
        dim = max(1, min(int(context_pca_dim), x.shape[1], n - 1))
        emb = PCA(n_components=dim, random_state=context_seed).fit_transform(np.asarray(x, dtype=np.float32))
        centers = MiniBatchKMeans(n_clusters=size, random_state=context_seed, batch_size=min(1024, n), n_init=10).fit(emb).cluster_centers_
        nearest, _ = pairwise_distances_argmin_min(centers, emb)
        selected: list[int] = []
        seen: set[int] = set()
        for idx in nearest.astype(int).tolist():
            if idx not in seen:
                selected.append(idx)
                seen.add(idx)
        rng = np.random.default_rng(context_seed)
        if len(selected) < size:
            first = selected[0] if selected else int(rng.integers(0, n))
            if not selected:
                selected.append(first)
                seen.add(first)
            distances = np.linalg.norm(emb - emb[first], axis=1)
            while len(selected) < size:
                candidates = np.setdiff1d(np.arange(n), np.fromiter(seen, dtype=np.int64), assume_unique=False)
                if candidates.size == 0:
                    break
                next_idx = int(candidates[np.argmax(distances[candidates])])
                selected.append(next_idx)
                seen.add(next_idx)
                distances = np.minimum(distances, np.linalg.norm(emb - emb[next_idx], axis=1))
        indices = np.asarray(selected[:size], dtype=np.int64)
    if save_dir is not None:
        save_dir.mkdir(parents=True, exist_ok=True)
        np.save(save_dir / "context_indices.npy", indices)
        with open(save_dir / "context_selection.json", "w", encoding="utf-8") as handle:
            json.dump(
                {
                    "context_size": int(size),
                    "context_pca_dim": int(context_pca_dim),
                    "context_seed": int(context_seed),
                    "label_free": True,
                    "n_selected": int(indices.size),
                },
                handle,
                indent=2,
            )
    return indices

