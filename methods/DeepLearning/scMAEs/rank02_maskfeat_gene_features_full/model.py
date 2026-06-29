from __future__ import annotations

import torch
import torch.nn as nn

class GenePatchMaskFeatMAE(nn.Module):
    def __init__(self):
        super().__init__()
