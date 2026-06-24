from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import torch
from torch.utils.data import Dataset


@dataclass
class EvaluationBundle:
    labels: np.ndarray
    label_names: np.ndarray
    label_key: str | None


class CAAMExpressionDataset(Dataset):
    """Training dataset. It intentionally never returns labels."""

    def __init__(
        self,
        x: np.ndarray,
        batch_code: np.ndarray | None = None,
        library_size: np.ndarray | None = None,
        zero_ratio: np.ndarray | None = None,
    ) -> None:
        self.x = torch.as_tensor(x, dtype=torch.float32)
        n = int(self.x.shape[0])
        self.batch_code = torch.as_tensor(
            np.zeros(n, dtype=np.int64) if batch_code is None else batch_code,
            dtype=torch.long,
        )
        self.library_size = torch.as_tensor(
            np.asarray(self.x.sum(dim=1), dtype=np.float32) if library_size is None else library_size,
            dtype=torch.float32,
        )
        self.zero_ratio = torch.as_tensor(
            np.asarray((self.x == 0).float().mean(dim=1), dtype=np.float32) if zero_ratio is None else zero_ratio,
            dtype=torch.float32,
        )

    def __len__(self) -> int:
        return int(self.x.shape[0])

    def __getitem__(self, idx: int) -> dict[str, torch.Tensor | int]:
        return {
            "index": int(idx),
            "x": self.x[idx],
            "batch_code": self.batch_code[idx],
            "library_size": self.library_size[idx],
            "zero_ratio": self.zero_ratio[idx],
        }


def assert_no_training_labels(batch: dict) -> None:
    forbidden = {"label", "labels", "cell_type", "true_cluster", "n_clusters"}
    present = forbidden.intersection(batch.keys())
    if present:
        raise AssertionError(f"Training batch contains forbidden label fields: {sorted(present)}")

