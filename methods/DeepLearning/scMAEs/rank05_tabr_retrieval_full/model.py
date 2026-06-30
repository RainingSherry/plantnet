from __future__ import annotations

import torch
import torch.nn as nn
import torch.nn.functional as F


class ResidualMLPBlock(nn.Module):
    def __init__(self, hidden_size: int, expansion: float, dropout: float) -> None:
        super().__init__()
        block_size = int(hidden_size * float(expansion))
        self.net = nn.Sequential(
            nn.LayerNorm(hidden_size),
            nn.Linear(hidden_size, block_size),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(block_size, hidden_size),
            nn.Dropout(dropout),
        )

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return x + self.net(x)


class TabRRetrievalScMAE(nn.Module):
    """Retrieval-augmented scMAE following the TabR context mechanism.

    The query is a masked expression vector. Candidate cells are encoded to keys,
    nearest contexts are retrieved in key space, and context values are aggregated
    with a TabR-style residual value term:

        value = candidate_value(candidate_x) + T(query_key - context_key)
        query_state = query_state + softmax(sim(query_key, context_key)) @ value
    """

    def __init__(
        self,
        num_genes: int,
        hidden_size: int = 128,
        encoder_blocks: int = 2,
        predictor_blocks: int = 2,
        expansion: float = 2.0,
        context_dropout: float = 0.1,
        dropout: float = 0.1,
    ) -> None:
        super().__init__()
        if num_genes <= 0:
            raise ValueError("num_genes must be positive")
        self.num_genes = int(num_genes)
        self.hidden_size = int(hidden_size)
        self.input = nn.Sequential(
            nn.Dropout(dropout),
            nn.Linear(num_genes, hidden_size),
            nn.LayerNorm(hidden_size),
            nn.Mish(),
        )
        self.encoder_blocks = nn.ModuleList(
            [ResidualMLPBlock(hidden_size, expansion, dropout) for _ in range(encoder_blocks)]
        )
        self.key = nn.Linear(hidden_size, hidden_size)
        self.value_encoder = nn.Sequential(
            nn.Linear(num_genes, hidden_size),
            nn.LayerNorm(hidden_size),
            nn.GELU(),
            nn.Linear(hidden_size, hidden_size),
        )
        self.context_transform = nn.Sequential(
            nn.Linear(hidden_size, int(hidden_size * expansion)),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(int(hidden_size * expansion), hidden_size, bias=False),
        )
        self.context_dropout = nn.Dropout(context_dropout)
        self.predictor_blocks = nn.ModuleList(
            [ResidualMLPBlock(hidden_size, expansion, dropout) for _ in range(predictor_blocks)]
        )
        self.mask_predictor = nn.Linear(hidden_size, num_genes)
        self.decoder = nn.Sequential(
            nn.LayerNorm(hidden_size),
            nn.GELU(),
            nn.Linear(hidden_size, num_genes),
        )

    def encode_key(self, x: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
        if x.ndim != 2 or x.shape[1] != self.num_genes:
            raise ValueError(f"x must be [batch, {self.num_genes}], got {tuple(x.shape)}")
        state = self.input(x)
        for block in self.encoder_blocks:
            state = block(state)
        key = self.key(state)
        return state, key

    @torch.no_grad()
    def search_context(
        self,
        query_key: torch.Tensor,
        candidate_key: torch.Tensor,
        context_size: int,
        query_indices: torch.Tensor | None = None,
    ) -> torch.Tensor:
        if query_key.ndim != 2 or candidate_key.ndim != 2 or query_key.shape[1] != candidate_key.shape[1]:
            raise ValueError("query_key and candidate_key must be [n, hidden] with matching hidden size")
        k_eff = min(int(context_size) + (1 if query_indices is not None else 0), candidate_key.shape[0])
        distances = torch.cdist(query_key.float(), candidate_key.float(), p=2.0)
        context_idx = distances.topk(k_eff, largest=False, dim=1).indices
        if query_indices is not None:
            if query_indices.ndim != 1 or query_indices.shape[0] != query_key.shape[0]:
                raise ValueError("query_indices must be [batch]")
            self_mask = context_idx == query_indices.view(-1, 1)
            distances_sel = distances.gather(1, context_idx).masked_fill(self_mask, float("inf"))
            order = distances_sel.argsort(dim=1)
            context_idx = context_idx.gather(1, order[:, : int(context_size)])
        elif context_idx.shape[1] > int(context_size):
            context_idx = context_idx[:, : int(context_size)]
        return context_idx

    def forward(
        self,
        x: torch.Tensor,
        candidate_x: torch.Tensor,
        context_size: int,
        query_indices: torch.Tensor | None = None,
    ) -> dict[str, torch.Tensor]:
        query_state, query_key = self.encode_key(x)
        with torch.no_grad():
            _, candidate_key = self.encode_key(candidate_x)
            context_idx = self.search_context(query_key.detach(), candidate_key.detach(), context_size, query_indices)
            context_key = candidate_key[context_idx].detach()
        context_x = candidate_x[context_idx]
        context_value = self.value_encoder(context_x)
        values = context_value + self.context_transform(query_key[:, None, :] - context_key)
        similarities = (
            -query_key.square().sum(dim=1, keepdim=True)
            + 2.0 * (query_key[:, None, :] @ context_key.transpose(1, 2)).squeeze(1)
            - context_key.square().sum(dim=2)
        )
        probs = self.context_dropout(F.softmax(similarities, dim=1))
        state = query_state + (probs[:, None, :] @ values).squeeze(1)
        for block in self.predictor_blocks:
            state = block(state)
        return {
            "embedding": state,
            "reconstruction": self.decoder(state),
            "mask_logits": self.mask_predictor(state),
            "context_idx": context_idx,
            "context_probs": probs,
        }

    @torch.no_grad()
    def feature(self, x: torch.Tensor, candidate_x: torch.Tensor, context_size: int, batch_size: int = 512) -> torch.Tensor:
        outputs = []
        for start in range(0, x.shape[0], int(batch_size)):
            xb = x[start:start + int(batch_size)]
            idx = torch.arange(start, start + xb.shape[0], device=x.device)
            outputs.append(self.forward(xb, candidate_x, context_size, idx)["embedding"])
        return torch.cat(outputs, dim=0)
