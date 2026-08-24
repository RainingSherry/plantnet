from __future__ import annotations

import numpy as np
import pandas as pd


def remap_desc_labels(labels) -> np.ndarray:
    """Map arbitrary DESC cluster ids to dense integer labels without category NaNs."""
    mapped = pd.factorize(np.asarray(labels), sort=True)[0]
    return mapped.astype(np.int64, copy=False)
