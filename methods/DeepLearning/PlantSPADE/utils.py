import json
import os
import random

import numpy as np
import torch


LABEL_CANDIDATES = [
    "cell_type",
    "Celltype",
    "celltype",
    "cell_label",
    "label",
    "labels",
    "Cluster",
    "cluster",
    "clusters",
    "Seurat_clusters",
    "cell_cluster",
    "paul15_clusters",
]


def ensure_dir(path: str) -> str:
    os.makedirs(path, exist_ok=True)
    return path


def save_json(payload: dict, path: str) -> None:
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, default=str)


def set_seed(seed: int = 42) -> None:
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False


def get_device(gpu: int = 0, no_cuda: bool = False) -> torch.device:
    if no_cuda or not torch.cuda.is_available():
        return torch.device("cpu")
    return torch.device(f"cuda:{gpu}")

