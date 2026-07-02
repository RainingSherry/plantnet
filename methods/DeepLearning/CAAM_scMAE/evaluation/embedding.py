from __future__ import annotations

import numpy as np
import torch
from torch.utils.data import DataLoader


@torch.no_grad()
def extract_embeddings(student, dataset, device: torch.device, batch_size: int, context_x: torch.Tensor | None = None, context_indices: torch.Tensor | None = None):
    student.eval()
    if context_x is not None and context_indices is not None:
        student.refresh_context_cache(context_x.to(device), context_indices.to(device))
    loader = DataLoader(dataset, batch_size=batch_size, shuffle=False, drop_last=False, num_workers=0)
    n = len(dataset)
    latent_dim = int(student.feature(next(iter(loader))["x"].to(device), indices=next(iter(loader))["index"].to(device)).shape[1])
    out = np.zeros((n, latent_dim), dtype=np.float32)
    for batch in loader:
        idx = batch["index"].long()
        z = student.feature(batch["x"].to(device), indices=idx.to(device))
        out[idx.numpy()] = z.detach().cpu().numpy().astype(np.float32)
    return np.nan_to_num(out, nan=0.0, posinf=0.0, neginf=0.0)

