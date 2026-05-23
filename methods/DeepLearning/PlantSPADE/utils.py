import json
import os
import random
from typing import Any

import numpy as np
import torch


LABEL_CANDIDATES = [
    "cell_type",
    "Celltype",
    "celltype",
    "cell_label",
    "label",
    "CellType",
    "cell_type1",
    "celltype_after",
    "seurat_clusters",
    "integrated_snn_res.0.3",
]


def set_seed(seed: int) -> None:
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed(seed)
        torch.cuda.manual_seed_all(seed)
        torch.backends.cudnn.deterministic = True
        torch.backends.cudnn.benchmark = False


def get_device(gpu: int = 0, no_cuda: bool = False) -> torch.device:
    use_cuda = (not no_cuda) and torch.cuda.is_available()
    return torch.device(f"cuda:{gpu}" if use_cuda else "cpu")


def ensure_dir(path: str) -> str:
    os.makedirs(path, exist_ok=True)
    return path


def save_json(data: dict[str, Any], path: str) -> None:
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(data, handle, indent=2, ensure_ascii=False)


def to_numpy(array_like: Any) -> np.ndarray:
    if isinstance(array_like, np.ndarray):
        return array_like
    if torch.is_tensor(array_like):
        return array_like.detach().cpu().numpy()
    return np.asarray(array_like)
