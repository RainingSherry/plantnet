from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import torch


@dataclass
class MaskState:
    mode: str
    base_rate: float
    gene_probs: np.ndarray
    stats: dict
    modules: list[np.ndarray] | None = None


def _minmax(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=np.float32)
    finite = np.isfinite(values)
    out = np.zeros_like(values, dtype=np.float32)
    if not np.any(finite):
        return out
    lo = float(np.min(values[finite]))
    hi = float(np.max(values[finite]))
    if hi <= lo:
        out[finite] = 0.5
    else:
        out[finite] = (values[finite] - lo) / (hi - lo)
    return out


def gene_statistics(data: np.ndarray) -> dict:
    x = np.asarray(data, dtype=np.float32)
    mean = np.mean(x, axis=0)
    var = np.var(x, axis=0)
    abs_mean = np.mean(np.abs(x), axis=0)
    dropout = np.mean(np.isclose(x, 0.0), axis=0)
    dispersion = var / (np.abs(mean) + 1e-6)
    marker_risk = np.maximum(_minmax(abs_mean), _minmax(var))
    return {
        "mean": mean,
        "var": var,
        "abs_mean": abs_mean,
        "dropout": dropout,
        "dispersion": dispersion,
        "marker_risk": marker_risk,
    }


def _build_modules(data: np.ndarray, module_size: int) -> list[np.ndarray]:
    n_genes = int(data.shape[1])
    module_size = max(2, int(module_size))
    if n_genes <= module_size:
        return [np.arange(n_genes, dtype=np.int64)]
    sample = np.asarray(data[: min(data.shape[0], 2048)], dtype=np.float32)
    sample = sample - sample.mean(axis=0, keepdims=True)
    try:
        _, _, vt = np.linalg.svd(sample, full_matrices=False)
        order = np.argsort(vt[0])
    except np.linalg.LinAlgError:
        order = np.arange(n_genes)
    return [order[start : start + module_size].astype(np.int64) for start in range(0, n_genes, module_size)]


def build_mask_state(
    data: np.ndarray,
    mode: str,
    base_rate: float,
    p_min: float,
    p_max: float,
    strength: float,
    module_size: int,
) -> MaskState:
    stats = gene_statistics(data)
    base = float(base_rate)
    mode = str(mode)
    if mode == "random":
        score = np.ones(data.shape[1], dtype=np.float32)
    elif mode == "variance_adaptive":
        score = 1.0 + float(strength) * (_minmax(stats["var"]) - 0.5)
    elif mode == "dropout_adaptive":
        balanced_dropout = 1.0 - np.abs(_minmax(stats["dropout"]) - 0.5) * 2.0
        score = 1.0 + float(strength) * (balanced_dropout - 0.5)
    elif mode == "marker_safe":
        score = 1.0 + float(strength) * (0.5 - stats["marker_risk"])
    elif mode == "module_block":
        module_var = _minmax(stats["var"])
        score = 1.0 + float(strength) * (module_var - 0.5)
    else:
        raise ValueError(f"Unknown mask mode: {mode}")
    probs = np.clip(base * score.astype(np.float32), float(p_min), float(p_max)).astype(np.float32)
    modules = _build_modules(data, module_size) if mode == "module_block" else None
    payload = {
        "mode": mode,
        "base_rate": base,
        "configured_mean_rate": float(np.mean(probs)),
        "p_min": float(np.min(probs)),
        "p_max": float(np.max(probs)),
        "n_modules": 0 if modules is None else int(len(modules)),
        "mean_variance": float(np.mean(stats["var"])),
        "mean_dropout": float(np.mean(stats["dropout"])),
        "mean_dispersion": float(np.mean(stats["dispersion"])),
        "mean_marker_risk": float(np.mean(stats["marker_risk"])),
    }
    return MaskState(mode=mode, base_rate=base, gene_probs=probs, stats=payload, modules=modules)


def apply_replacement_noise(
    x: torch.Tensor,
    mask_state: MaskState | None,
    mask_ratio: float,
) -> tuple[torch.Tensor, torch.Tensor]:
    if mask_state is None:
        probs = float(mask_ratio) * torch.ones_like(x)
    elif mask_state.mode == "module_block" and mask_state.modules:
        probs = torch.zeros_like(x)
        gene_probs = torch.as_tensor(mask_state.gene_probs, dtype=x.dtype, device=x.device)
        for module in mask_state.modules:
            module_idx = torch.as_tensor(module, dtype=torch.long, device=x.device)
            module_prob = torch.mean(gene_probs[module_idx]).clamp(0.0, 1.0)
            cell_draw = torch.bernoulli(module_prob * torch.ones((x.shape[0], 1), device=x.device, dtype=x.dtype))
            probs[:, module_idx] = cell_draw
        should_swap = probs
        replacement = x if x.shape[0] <= 1 else x[torch.randperm(x.shape[0], device=x.device)]
        corrupted = torch.where(should_swap.bool(), replacement, x)
        return corrupted, (corrupted != x).float()
    else:
        gene_probs = torch.as_tensor(mask_state.gene_probs, dtype=x.dtype, device=x.device)
        probs = gene_probs.view(1, -1).expand_as(x)
    should_swap = torch.bernoulli(probs.clamp(0.0, 1.0))
    replacement = x if x.shape[0] <= 1 else x[torch.randperm(x.shape[0], device=x.device)]
    corrupted = torch.where(should_swap.bool(), replacement, x)
    mask = (corrupted != x).float()
    return corrupted, mask

