from __future__ import annotations

import torch
import torch.nn as nn
import torch.nn.functional as F


class SharedGeneContext(nn.Module):
    def __init__(self, n_genes: int, token_dim: int) -> None:
        super().__init__()
        self.gene_embedding = nn.Embedding(n_genes, token_dim)
        self.stat_encoder = nn.Sequential(
            nn.Linear(4, token_dim),
            nn.LayerNorm(token_dim),
            nn.Mish(),
            nn.Linear(token_dim, token_dim),
        )

    def forward(self, gene_stats: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
        gene_ids = torch.arange(gene_stats.shape[0], device=gene_stats.device)
        return self.gene_embedding(gene_ids), self.stat_encoder(gene_stats.float())


class PrototypeAttention(nn.Module):
    def __init__(self, query_dim: int, proto_dim: int, token_dim: int) -> None:
        super().__init__()
        self.query = nn.Linear(query_dim, token_dim)
        self.key = nn.Linear(proto_dim, token_dim)
        self.value = nn.Linear(proto_dim, token_dim)
        self.scale = token_dim ** -0.5

    def forward(self, query: torch.Tensor, prototypes: torch.Tensor) -> torch.Tensor:
        q = self.query(query).unsqueeze(1)
        k = self.key(prototypes).unsqueeze(0)
        v = self.value(prototypes).unsqueeze(0)
        weights = torch.softmax(torch.matmul(q, k.transpose(1, 2)) * self.scale, dim=-1)
        return torch.matmul(weights, v).squeeze(1)


def straight_through_topk(logits: torch.Tensor, mask_ratio: float, temperature: float) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
    n_genes = int(logits.shape[1])
    k = max(1, min(n_genes, int(float(mask_ratio) * n_genes)))
    topk = torch.topk(logits, k=k, dim=1).indices
    hard = torch.zeros_like(logits).scatter_(1, topk, 1.0)
    soft = torch.softmax(logits / max(float(temperature), 1.0e-6), dim=1) * float(k)
    soft = soft.clamp(0.0, 1.0)
    st = hard + soft - soft.detach()
    return hard, soft, st


class APAGenerator(nn.Module):
    def __init__(self, *, n_genes: int, token_dim: int, proto_dim: int, attention_heads: int, dropout: float) -> None:
        super().__init__()
        if token_dim % attention_heads != 0:
            raise ValueError("token_dim must be divisible by attention_heads")
        self.pair_encoder = nn.Sequential(
            nn.Linear(4, token_dim),
            nn.LayerNorm(token_dim),
            nn.Mish(),
            nn.Linear(token_dim, token_dim),
        )
        self.fusion = nn.Sequential(
            nn.Linear(token_dim * 3, token_dim),
            nn.LayerNorm(token_dim),
            nn.Mish(),
        )
        self.gene_attention = nn.MultiheadAttention(
            embed_dim=token_dim,
            num_heads=attention_heads,
            dropout=dropout,
            batch_first=True,
        )
        self.proto_attention = PrototypeAttention(token_dim, proto_dim, token_dim)
        self.mask_fusion = nn.Sequential(
            nn.Linear(token_dim * 3 + 2, token_dim),
            nn.LayerNorm(token_dim),
            nn.Mish(),
            nn.Linear(token_dim, 1),
        )
        self.n_genes = int(n_genes)

    def forward(
        self,
        x: torch.Tensor,
        replacement: torch.Tensor,
        effective: torch.Tensor,
        gene_vec: torch.Tensor,
        stat_vec: torch.Tensor,
        prototypes: torch.Tensor,
        *,
        mask_ratio: float,
        temperature: float,
    ) -> dict[str, torch.Tensor]:
        delta = (replacement - x).abs()
        pair_input = torch.stack([x, replacement, delta, effective.float()], dim=-1)
        pair_vec = self.pair_encoder(pair_input)
        gene_context = torch.cat([gene_vec, stat_vec], dim=-1).unsqueeze(0).expand(x.shape[0], -1, -1)
        h_gen = self.fusion(torch.cat([pair_vec, gene_context], dim=-1))
        h_gene, _ = self.gene_attention(h_gen, h_gen, h_gen, need_weights=False)
        cell_summary = h_gene.mean(dim=1)
        proto_context = self.proto_attention(cell_summary, prototypes).unsqueeze(1).expand_as(h_gene)
        h_mask = torch.cat([h_gen, h_gene, proto_context, delta.unsqueeze(-1), effective.float().unsqueeze(-1)], dim=-1)
        logits = self.mask_fusion(h_mask).squeeze(-1)
        hard, soft, st = straight_through_topk(logits, mask_ratio, temperature)
        return {"logits": logits, "mask_hard": hard, "mask_soft": soft, "mask_st": st, "delta": delta}


class APAStudent(nn.Module):
    def __init__(self, *, n_genes: int, token_dim: int, cell_dim: int, proto_dim: int, attention_heads: int, dropout: float) -> None:
        super().__init__()
        if token_dim % attention_heads != 0:
            raise ValueError("token_dim must be divisible by attention_heads")
        self.value_encoder = nn.Sequential(
            nn.Linear(1, token_dim),
            nn.LayerNorm(token_dim),
            nn.Mish(),
            nn.Linear(token_dim, token_dim),
        )
        self.fusion = nn.Sequential(
            nn.Linear(token_dim * 3, token_dim),
            nn.LayerNorm(token_dim),
            nn.Mish(),
        )
        self.cell_token = nn.Parameter(torch.zeros(1, 1, token_dim))
        nn.init.normal_(self.cell_token, std=0.02)
        self.gene_attention = nn.MultiheadAttention(
            embed_dim=token_dim,
            num_heads=attention_heads,
            dropout=dropout,
            batch_first=True,
        )
        self.proto_attention = PrototypeAttention(token_dim, proto_dim, token_dim)
        self.pool = nn.Sequential(
            nn.Linear(token_dim * 2, cell_dim),
            nn.LayerNorm(cell_dim),
            nn.Mish(),
            nn.Linear(cell_dim, cell_dim),
        )
        self.mask_head = nn.Sequential(
            nn.Linear(token_dim + cell_dim, token_dim),
            nn.Mish(),
            nn.Linear(token_dim, 1),
        )
        self.decoder = nn.Sequential(
            nn.Linear(token_dim + cell_dim + 1 + token_dim, token_dim),
            nn.Mish(),
            nn.Linear(token_dim, 1),
        )
        self.n_genes = int(n_genes)
        self.cell_dim = int(cell_dim)

    def encode_tokens(
        self,
        x_tilde: torch.Tensor,
        gene_vec: torch.Tensor,
        stat_vec: torch.Tensor,
        prototypes: torch.Tensor,
    ) -> tuple[torch.Tensor, torch.Tensor]:
        value_vec = self.value_encoder(x_tilde.unsqueeze(-1))
        gene_context = torch.cat([gene_vec, stat_vec], dim=-1).unsqueeze(0).expand(x_tilde.shape[0], -1, -1)
        h0 = self.fusion(torch.cat([value_vec, gene_context], dim=-1))
        cell_token = self.cell_token.expand(x_tilde.shape[0], -1, -1)
        encoded, _ = self.gene_attention(
            torch.cat([cell_token, h0], dim=1),
            torch.cat([cell_token, h0], dim=1),
            torch.cat([cell_token, h0], dim=1),
            need_weights=False,
        )
        c_raw = encoded[:, 0, :]
        h_gene = encoded[:, 1:, :]
        p_context = self.proto_attention(c_raw, prototypes)
        z = self.pool(torch.cat([c_raw, p_context], dim=-1))
        return z, h_gene

    def forward(
        self,
        x_tilde: torch.Tensor,
        gene_vec: torch.Tensor,
        stat_vec: torch.Tensor,
        prototypes: torch.Tensor,
    ) -> dict[str, torch.Tensor]:
        z, h_gene = self.encode_tokens(x_tilde, gene_vec, stat_vec, prototypes)
        z_gene = z.unsqueeze(1).expand(-1, x_tilde.shape[1], -1)
        pred_mask_logits = self.mask_head(torch.cat([h_gene, z_gene], dim=-1)).squeeze(-1)
        gene_vec_batch = gene_vec.unsqueeze(0).expand(x_tilde.shape[0], -1, -1)
        recon_input = torch.cat([h_gene, z_gene, pred_mask_logits.unsqueeze(-1), gene_vec_batch], dim=-1)
        reconstruction = self.decoder(recon_input).squeeze(-1)
        return {"z": z, "gene_tokens": h_gene, "mask_logits": pred_mask_logits, "x_recon": reconstruction}

    def feature(
        self,
        x_clean: torch.Tensor,
        gene_vec: torch.Tensor,
        stat_vec: torch.Tensor,
        prototypes: torch.Tensor,
    ) -> torch.Tensor:
        z, _ = self.encode_tokens(x_clean, gene_vec, stat_vec, prototypes)
        return z


class APAModel(nn.Module):
    def __init__(
        self,
        *,
        n_genes: int,
        token_dim: int,
        cell_dim: int,
        proto_dim: int,
        attention_heads: int,
        dropout: float,
    ) -> None:
        super().__init__()
        self.shared = SharedGeneContext(n_genes, token_dim)
        self.generator = APAGenerator(
            n_genes=n_genes,
            token_dim=token_dim,
            proto_dim=proto_dim,
            attention_heads=attention_heads,
            dropout=dropout,
        )
        self.student = APAStudent(
            n_genes=n_genes,
            token_dim=token_dim,
            cell_dim=cell_dim,
            proto_dim=proto_dim,
            attention_heads=attention_heads,
            dropout=dropout,
        )

    def shared_context(self, gene_stats: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
        return self.shared(gene_stats)


def grad_norm(module: nn.Module | None) -> float:
    if module is None:
        return 0.0
    total = 0.0
    for param in module.parameters():
        if param.grad is not None:
            total += float(param.grad.detach().pow(2).sum().cpu())
    return float(total ** 0.5)


def freeze(module: nn.Module, frozen: bool) -> None:
    for param in module.parameters():
        param.requires_grad_(not frozen)
