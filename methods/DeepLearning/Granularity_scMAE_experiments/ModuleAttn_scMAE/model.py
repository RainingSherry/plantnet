from __future__ import annotations

import torch
from torch import nn


class ModuleAttnScMAE(nn.Module):
    """scMAE backbone whose encoder is replaced by gene-module tokenization +
    module-level self-attention, then DEC cluster centers on the latent.

    Only the ENCODER differs from AdaptiveSwitch (which uses a plain MLP). DEC
    head + decoder + student_q are identical, so any change is attributable to
    the gene-module structure.

    Encoder flow:
      x (n, G) --grouped masked proj--> M module tokens (n, M, d)
                --+ module posemb, LayerNorm--> self-attention (L layers)
                --flatten (n, M*d) --> Linear --> latent (n, hidden)

    `module_of` : LongTensor (G,) giving each gene's module id in [0, M).
    Grouped projection = a (G x (M*d)) weight masked so gene g only feeds the
    d slots of its own module -> respects co-expression grouping, solves the
    "genes have no natural order" problem that breaks naive 1D conv.
    """

    def __init__(self, num_genes: int, n_clusters: int, module_of: torch.Tensor,
                 n_modules: int, token_dim: int = 16, hidden_size: int = 128,
                 n_heads: int = 4, n_layers: int = 2, dropout: float = 0.05,
                 use_attn: bool = True):
        super().__init__()
        self.num_genes = int(num_genes)
        self.n_clusters = int(n_clusters)
        self.M = int(n_modules)
        self.d = int(token_dim)
        self.hidden_size = int(hidden_size)
        self.use_attn = bool(use_attn)

        # grouped masked projection: gene g -> only its module's d slots
        self.in_proj = nn.Linear(self.num_genes, self.M * self.d)
        mask = torch.zeros(self.M * self.d, self.num_genes)
        mo = module_of.long()
        for g in range(self.num_genes):
            m = int(mo[g])
            mask[m * self.d:(m + 1) * self.d, g] = 1.0
        self.register_buffer("proj_mask", mask)  # (M*d, G)

        self.module_posemb = nn.Parameter(torch.randn(1, self.M, self.d) * 0.02)
        self.token_norm = nn.LayerNorm(self.d)
        self.input_dropout = nn.Dropout(dropout)

        if self.use_attn:
            layer = nn.TransformerEncoderLayer(
                d_model=self.d, nhead=n_heads, dim_feedforward=self.d * 4,
                dropout=dropout, batch_first=True, activation="gelu")
            self.attn = nn.TransformerEncoder(layer, num_layers=n_layers)
        else:
            self.attn = None

        self.to_latent = nn.Sequential(
            nn.Linear(self.M * self.d, hidden_size),
            nn.LayerNorm(hidden_size),
            nn.Mish(),
            nn.Linear(hidden_size, hidden_size),
        )
        self.mask_predictor = nn.Linear(hidden_size, self.num_genes)
        self.decoder = nn.Linear(hidden_size + self.num_genes, self.num_genes)
        self.cluster_centers = nn.Parameter(torch.randn(self.n_clusters, hidden_size) * 0.02)

    def encode(self, x: torch.Tensor) -> torch.Tensor:
        x = self.input_dropout(x)
        masked_w = self.in_proj.weight * self.proj_mask  # enforce grouping
        tokens = torch.nn.functional.linear(x, masked_w, self.in_proj.bias)  # (n, M*d)
        tokens = tokens.view(-1, self.M, self.d) + self.module_posemb
        tokens = self.token_norm(tokens)
        if self.attn is not None:
            tokens = self.attn(tokens)
        return self.to_latent(tokens.reshape(tokens.shape[0], -1))

    def forward(self, x: torch.Tensor) -> dict[str, torch.Tensor]:
        latent = self.encode(x)
        mask_logits = self.mask_predictor(latent)
        reconstruction = self.decoder(torch.cat([latent, mask_logits], dim=1))
        q = self.student_q(latent)
        return {"latent": latent, "mask_logits": mask_logits, "reconstruction": reconstruction, "cluster_q": q}

    def student_q(self, latent: torch.Tensor) -> torch.Tensor:
        dist = torch.cdist(latent, self.cluster_centers).pow(2)
        q = 1.0 / (1.0 + dist)
        return q / q.sum(dim=1, keepdim=True).clamp_min(1e-8)

    @torch.no_grad()
    def feature(self, x: torch.Tensor) -> torch.Tensor:
        return self.encode(x)

    def random_mask(self, x: torch.Tensor, mask_prob: float) -> tuple[torch.Tensor, torch.Tensor]:
        mask = (torch.rand_like(x) < float(mask_prob)).float()
        empty = mask.sum(dim=1) == 0
        if bool(empty.any()):
            cols = torch.randint(0, x.shape[1], (int(empty.sum()),), device=x.device)
            mask[empty, cols] = 1.0
        return x.masked_fill(mask.bool(), 0.0), mask

    @torch.no_grad()
    def initialize_centers(self, centers: torch.Tensor) -> None:
        self.cluster_centers.copy_(centers)

    @staticmethod
    def sharpen(q: torch.Tensor) -> torch.Tensor:
        weight = q.pow(2) / q.sum(dim=0, keepdim=True).clamp_min(1e-8)
        return weight / weight.sum(dim=1, keepdim=True).clamp_min(1e-8)
