from __future__ import annotations

import numpy as np
import torch
from torch import nn


class PlainMLPEncoder(nn.Module):
    """The rank13 / GatedNeighborMix MLP encoder, verbatim. Strongest proven
    backbone; used as the base for the marker-aware masking experiment so the
    only new signal is the per-gene reconstruction weighting / biased masking."""

    def __init__(self, n_genes: int, latent_dim: int = 128, hidden: int = 256, dropout: float = 0.05) -> None:
        super().__init__()
        self.net = nn.Sequential(
            nn.Dropout(dropout),
            nn.Linear(n_genes, hidden),
            nn.LayerNorm(hidden),
            nn.Mish(),
            nn.Linear(hidden, latent_dim),
            nn.LayerNorm(latent_dim),
            nn.Mish(),
            nn.Linear(latent_dim, latent_dim),
        )

    def forward(self, x: torch.Tensor, r: torch.Tensor | None = None) -> torch.Tensor:
        return self.net(x)


class CellAxisOnMLPEncoder(nn.Module):
    """scMAE-level MLP base + reliability-gated prototype cell-axis on TOP of it.

    The gene-module tokenizer is dropped entirely (it was empirically dead
    weight: more modules made clustering worse). The only new mechanism is the
    cell axis -- cells refine their embedding by cross-attending to K learnable
    prototypes -- and it now runs on the full-gene MLP features, not on lossy
    pooled tokens. A zero-initialized scalar gate `g` means training starts as
    pure scMAE MLP and the cell axis is only mixed in if it actually helps, so
    the encoder can never underperform the MLP base by construction.

    z = proj([h, r_i * g * PrototypeAttn(h, P)]),  h = MLP(x)
    """

    def __init__(
        self,
        n_genes: int,
        latent_dim: int = 128,
        mlp_hidden: int = 256,
        n_prototypes: int = 16,
        proto_dim: int = 128,
        dropout: float = 0.1,
    ) -> None:
        super().__init__()
        self.mlp_base = nn.Sequential(
            nn.Dropout(dropout),
            nn.Linear(n_genes, mlp_hidden),
            nn.LayerNorm(mlp_hidden),
            nn.Mish(),
            nn.Linear(mlp_hidden, proto_dim),
            nn.LayerNorm(proto_dim),
            nn.Mish(),
        )
        self.prototypes = nn.Parameter(torch.randn(n_prototypes, proto_dim) * 0.02)
        self.q = nn.Linear(proto_dim, proto_dim)
        self.kv = nn.Linear(proto_dim, proto_dim * 2)
        self.scale = proto_dim ** -0.5
        # zero-init scalar gate: starts as pure MLP, learns to add the cell axis
        self.axis_gate = nn.Parameter(torch.zeros(1))
        self.proj = nn.Sequential(nn.Linear(proto_dim * 2, latent_dim), nn.LayerNorm(latent_dim))

    def forward(self, x: torch.Tensor, r: torch.Tensor | None = None) -> torch.Tensor:
        h = self.mlp_base(x)                                 # [B, d] scMAE-level features
        qh = self.q(h)
        k, v = self.kv(self.prototypes).chunk(2, dim=-1)     # [K, d], [K, d]
        attn = torch.softmax(qh @ k.t() * self.scale, dim=-1)  # [B, K]
        ctx = attn @ v                                       # [B, d] cell-axis signal
        ctx = self.axis_gate * ctx                           # global learnable gate (0-init)
        if r is not None:                                    # per-cell reliability gate
            ctx = r.view(-1, 1) * ctx
        return self.proj(torch.cat([h, ctx], dim=-1))        # [B, latent_dim]


class DualAxisEncoder(nn.Module):
    """Gene-axis (always on) + prototype cell-axis (reliability-gated) encoder.

    - Gene axis: expression is pooled into M gene modules (shared gene-row
      weights via a hard [G, M] assignment), tokenized, then run through a few
      self-attention layers. This is the per-cell representation base and is
      NEVER gated -- it plays the role of the self-anchored scMAE backbone.
    - Cell axis: the gene-axis cell summary cross-attends to K learnable
      prototypes (K << N, so it scales). This is the "cells inform each other"
      signal. It is multiplied by the per-cell reliability r_i, so rare/boundary
      cells (r->0) fall back to the pure gene-axis path (= pure scMAE behavior),
      exactly where neighbor smoothing is unsafe.

    Output: latent [B, latent_dim], a drop-in replacement for the MLP encoder in
    GatedNeighborMixScMAE, so every downstream head is unchanged.
    """

    def __init__(
        self,
        n_genes: int,
        assignment: np.ndarray,
        token_dim: int = 48,
        latent_dim: int = 128,
        n_prototypes: int = 16,
        heads: int = 4,
        gene_layers: int = 2,
        dropout: float = 0.1,
        use_mlp_base: bool = True,
        mlp_hidden: int = 256,
    ) -> None:
        super().__init__()
        A = np.asarray(assignment, dtype=np.float32)
        denom = A.sum(axis=0, keepdims=True)
        denom[denom == 0.0] = 1.0
        self.register_buffer("assignment", torch.as_tensor(A / denom, dtype=torch.float32))
        self.n_modules = int(A.shape[1])
        self.value_mlp = nn.Sequential(
            nn.Linear(1, token_dim), nn.Mish(), nn.Linear(token_dim, token_dim)
        )
        self.module_emb = nn.Parameter(torch.randn(self.n_modules, token_dim) * 0.02)
        self.gene_axis = nn.ModuleList(
            [
                nn.TransformerEncoderLayer(
                    d_model=token_dim,
                    nhead=heads,
                    dim_feedforward=token_dim * 4,
                    dropout=dropout,
                    activation="gelu",
                    batch_first=True,
                )
                for _ in range(gene_layers)
            ]
        )
        self.prototypes = nn.Parameter(torch.randn(n_prototypes, token_dim) * 0.02)
        self.q = nn.Linear(token_dim, token_dim)
        self.kv = nn.Linear(token_dim, token_dim * 2)
        self.scale = token_dim ** -0.5
        # Optional full-gene MLP base path (fusion). Guarantees scMAE-level
        # per-gene fidelity; the dual-axis (gene-module + cell-axis) signal is
        # concatenated as an augmentation so fine markers are never pooled away.
        self.use_mlp_base = bool(use_mlp_base)
        if self.use_mlp_base:
            self.mlp_base = nn.Sequential(
                nn.Dropout(dropout),
                nn.Linear(n_genes, mlp_hidden),
                nn.LayerNorm(mlp_hidden),
                nn.Mish(),
                nn.Linear(mlp_hidden, token_dim),
                nn.LayerNorm(token_dim),
                nn.Mish(),
            )
            proj_in = token_dim * 3
        else:
            proj_in = token_dim * 2
        self.proj = nn.Sequential(nn.Linear(proj_in, latent_dim), nn.LayerNorm(latent_dim))

    def forward(self, x: torch.Tensor, r: torch.Tensor | None = None) -> torch.Tensor:
        u = x @ self.assignment                              # [B, M] module expression
        tok = self.value_mlp(u.unsqueeze(-1)) + self.module_emb.unsqueeze(0)  # [B, M, d]
        for blk in self.gene_axis:                           # gene-axis self-attention
            tok = blk(tok)
        cell = tok.mean(dim=1)                               # [B, d] gene-axis-only base
        qh = self.q(cell)                                    # [B, d]
        k, v = self.kv(self.prototypes).chunk(2, dim=-1)     # [K, d], [K, d]
        attn = torch.softmax(qh @ k.t() * self.scale, dim=-1)  # [B, K]
        ctx = attn @ v                                       # [B, d] cell-axis signal
        if r is not None:                                    # reliability gate
            ctx = r.view(-1, 1) * ctx
        parts = [cell, ctx]
        if self.use_mlp_base:
            parts.append(self.mlp_base(x))                   # [B, d] full-gene base
        return self.proj(torch.cat(parts, dim=-1))           # [B, latent_dim]


class DualAxisGatedScMAE(nn.Module):
    """rank13 DEC scMAE backbone with the MLP encoder swapped for DualAxisEncoder.

    mask_predictor / decoder / cluster_centers / student_q / random_mask /
    target_distribution are identical to GatedNeighborMixScMAE so the proven DEC
    behavior and the gated NeighborMix training loop stay unchanged. Only the
    encoder path (and the optional per-cell reliability r) is new.
    """

    def __init__(
        self,
        num_genes: int,
        n_clusters: int,
        assignment: np.ndarray,
        hidden_size: int = 128,
        dropout: float = 0.05,
        token_dim: int = 48,
        n_prototypes: int = 16,
        heads: int = 4,
        gene_layers: int = 2,
        attn_dropout: float = 0.1,
        use_mlp_base: bool = True,
        mlp_hidden: int = 256,
        encoder_type: str = "dualaxis",
        gene_weight: np.ndarray | None = None,
    ) -> None:
        super().__init__()
        if encoder_type not in {"dualaxis", "cellaxis_on_mlp", "mlp"}:
            raise ValueError(f"Unknown encoder_type={encoder_type!r}")
        self.num_genes = int(num_genes)
        self.n_clusters = int(n_clusters)
        self.hidden_size = int(hidden_size)
        self.encoder_type = encoder_type
        # label-free per-gene marker weight (bimodality), for biased masking.
        if gene_weight is None:
            gw = torch.ones(self.num_genes, dtype=torch.float32)
        else:
            gw = torch.as_tensor(np.asarray(gene_weight, dtype=np.float32))
        self.register_buffer("gene_weight", gw)
        if encoder_type == "cellaxis_on_mlp":
            self.encoder = CellAxisOnMLPEncoder(
                n_genes=num_genes,
                latent_dim=hidden_size,
                mlp_hidden=mlp_hidden,
                n_prototypes=n_prototypes,
                proto_dim=hidden_size,
                dropout=attn_dropout,
            )
        elif encoder_type == "mlp":
            self.encoder = PlainMLPEncoder(
                n_genes=num_genes, latent_dim=hidden_size, hidden=mlp_hidden, dropout=dropout
            )
        else:
            self.encoder = DualAxisEncoder(
                n_genes=num_genes,
                assignment=assignment,
                token_dim=token_dim,
                latent_dim=hidden_size,
                n_prototypes=n_prototypes,
                heads=heads,
                gene_layers=gene_layers,
                dropout=attn_dropout,
                use_mlp_base=use_mlp_base,
                mlp_hidden=mlp_hidden,
            )
        self.mask_predictor = nn.Linear(hidden_size, self.num_genes)
        self.decoder = nn.Linear(hidden_size + self.num_genes, self.num_genes)
        self.cluster_centers = nn.Parameter(torch.randn(self.n_clusters, hidden_size) * 0.02)

    def forward(self, x: torch.Tensor, r: torch.Tensor | None = None) -> dict[str, torch.Tensor]:
        latent = self.encoder(x, r)
        mask_logits = self.mask_predictor(latent)
        reconstruction = self.decoder(torch.cat([latent, mask_logits], dim=1))
        q = self.student_q(latent)
        return {"latent": latent, "mask_logits": mask_logits, "reconstruction": reconstruction, "cluster_q": q}

    def student_q(self, latent: torch.Tensor) -> torch.Tensor:
        dist = torch.cdist(latent, self.cluster_centers).pow(2)
        q = 1.0 / (1.0 + dist)
        return q / q.sum(dim=1, keepdim=True).clamp_min(1e-8)

    @torch.no_grad()
    def feature(self, x: torch.Tensor, r: torch.Tensor | None = None) -> torch.Tensor:
        return self.encoder(x, r)

    def random_mask(self, x: torch.Tensor, mask_prob: float, marker_bias: float = 0.0) -> tuple[torch.Tensor, torch.Tensor]:
        """Zero-masking (rank13). If marker_bias>0, per-gene mask probability is
        tilted toward high-marker (bimodal) genes: p_g = mask_prob * (1 + bias*(w_g_norm - mean)),
        so marker genes are masked more often and the encoder must reconstruct
        them from context (anti-smoothing). marker_bias=0 recovers uniform masking."""
        if marker_bias > 0.0:
            w = self.gene_weight.to(x.device)
            wn = w / w.mean().clamp_min(1e-8)                 # mean 1
            p = (float(mask_prob) * (1.0 + float(marker_bias) * (wn - 1.0))).clamp(0.01, 0.99)
            mask = (torch.rand_like(x) < p.view(1, -1)).float()
        else:
            mask = (torch.rand_like(x) < float(mask_prob)).float()
        empty = mask.sum(dim=1) == 0
        if bool(empty.any()):
            cols = torch.randint(0, x.shape[1], (int(empty.sum()),), device=x.device)
            mask[empty, cols] = 1.0
        return x.masked_fill(mask.bool(), 0.0), mask

    @torch.no_grad()
    def initialize_centers(self, centers: torch.Tensor) -> None:
        if centers.shape != self.cluster_centers.shape:
            raise ValueError(f"center shape {tuple(centers.shape)} != {tuple(self.cluster_centers.shape)}")
        self.cluster_centers.copy_(centers)

    @staticmethod
    def target_distribution(q: torch.Tensor) -> torch.Tensor:
        weight = q.pow(2) / q.sum(dim=0, keepdim=True).clamp_min(1e-8)
        return weight / weight.sum(dim=1, keepdim=True).clamp_min(1e-8)
