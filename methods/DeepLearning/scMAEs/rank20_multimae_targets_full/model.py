from __future__ import annotations

import itertools
from collections import OrderedDict
from typing import Mapping

import torch
import torch.nn as nn
from torch.distributions.dirichlet import Dirichlet


TASK_ORDER = ("expr", "rank", "stat")


class TransformerBlock(nn.Module):
    def __init__(self, dim: int, num_heads: int, dropout: float, mlp_ratio: int = 4) -> None:
        super().__init__()
        if dim % num_heads != 0:
            raise ValueError("dim must be divisible by num_heads")
        self.norm1 = nn.LayerNorm(dim)
        self.attn = nn.MultiheadAttention(dim, num_heads, dropout=dropout, batch_first=True)
        self.drop = nn.Dropout(dropout)
        self.norm2 = nn.LayerNorm(dim)
        self.mlp = nn.Sequential(
            nn.Linear(dim, dim * mlp_ratio),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(dim * mlp_ratio, dim),
            nn.Dropout(dropout),
        )

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        q = self.norm1(x)
        attn, _ = self.attn(q, q, q, need_weights=False)
        x = x + self.drop(attn)
        return x + self.mlp(self.norm2(x))


class CrossAttentionBlock(nn.Module):
    def __init__(self, dim: int, num_heads: int, dropout: float, mlp_ratio: int = 4) -> None:
        super().__init__()
        if dim % num_heads != 0:
            raise ValueError("dim must be divisible by num_heads")
        self.query_norm = nn.LayerNorm(dim)
        self.context_norm = nn.LayerNorm(dim)
        self.attn = nn.MultiheadAttention(dim, num_heads, dropout=dropout, batch_first=True)
        self.drop = nn.Dropout(dropout)
        self.out_norm = nn.LayerNorm(dim)
        self.mlp = nn.Sequential(
            nn.Linear(dim, dim * mlp_ratio),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(dim * mlp_ratio, dim),
            nn.Dropout(dropout),
        )

    def forward(self, queries: torch.Tensor, context: torch.Tensor) -> torch.Tensor:
        attn, _ = self.attn(
            self.query_norm(queries),
            self.context_norm(context),
            self.context_norm(context),
            need_weights=False,
        )
        queries = queries + self.drop(attn)
        return queries + self.mlp(self.out_norm(queries))


class TaskInputAdapter(nn.Module):
    def __init__(self, input_dim: int, hidden_size: int, dropout: float) -> None:
        super().__init__()
        if input_dim <= 0:
            raise ValueError("input_dim must be positive")
        self.proj = nn.Linear(input_dim, hidden_size)
        self.norm = nn.LayerNorm(hidden_size)
        self.drop = nn.Dropout(dropout)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        if x.ndim != 3:
            raise ValueError(f"task input must be [batch, patches, features], got {tuple(x.shape)}")
        return self.drop(self.norm(self.proj(x)))


class TaskOutputDecoder(nn.Module):
    """MultiMAE-style task decoder: unshuffle, add mask tokens, then cross-attend."""

    def __init__(
        self,
        task: str,
        output_dim: int,
        num_patches: int,
        hidden_size: int,
        decoder_size: int,
        decoder_depth: int,
        decoder_heads: int,
        dropout: float,
        context_tasks: tuple[str, ...] = TASK_ORDER,
    ) -> None:
        super().__init__()
        if task not in context_tasks:
            raise ValueError(f"unknown task {task!r}")
        if decoder_size % decoder_heads != 0:
            raise ValueError("decoder_size must be divisible by decoder_heads")
        self.task = task
        self.num_patches = int(num_patches)
        self.output_dim = int(output_dim)
        self.proj_context = nn.Linear(hidden_size, decoder_size)
        self.mask_token = nn.Parameter(torch.zeros(1, 1, decoder_size))
        self.position = nn.Parameter(torch.randn(1, self.num_patches, decoder_size) * 0.02)
        self.task_embeddings = nn.ParameterDict(
            {name: nn.Parameter(torch.randn(1, 1, decoder_size) * 0.02) for name in context_tasks}
        )
        self.cross_attn = CrossAttentionBlock(decoder_size, decoder_heads, dropout)
        self.decoder_blocks = nn.ModuleList(
            [TransformerBlock(decoder_size, decoder_heads, dropout) for _ in range(decoder_depth)]
        )
        self.out_norm = nn.LayerNorm(decoder_size)
        self.out_proj = nn.Linear(decoder_size, self.output_dim)

    def _context_embeddings(self, input_info: Mapping[str, object], batch_size: int, device: torch.device) -> torch.Tensor:
        pieces = []
        task_infos = input_info["tasks"]
        if not isinstance(task_infos, Mapping):
            raise ValueError("input_info['tasks'] must be a mapping")
        for name, info in task_infos.items():
            num_tokens = int(info["num_tokens"])
            if num_tokens != self.num_patches:
                raise ValueError(f"task {name} has {num_tokens} tokens; expected {self.num_patches}")
            task_emb = self.task_embeddings[name].expand(batch_size, num_tokens, -1)
            pieces.append(task_emb + self.position.to(device=device))
        return torch.cat(pieces, dim=1)

    def forward(
        self,
        encoder_tokens: torch.Tensor,
        input_info: Mapping[str, object],
        ids_keep: torch.Tensor,
        ids_restore: torch.Tensor,
    ) -> torch.Tensor:
        if encoder_tokens.ndim != 3:
            raise ValueError(f"encoder_tokens must be [batch, tokens, hidden], got {tuple(encoder_tokens.shape)}")
        batch_size = encoder_tokens.shape[0]
        context_tokens = self.proj_context(encoder_tokens)
        num_global = int(input_info["num_global_tokens"])
        context_without_global = context_tokens[:, :-num_global] if num_global else context_tokens
        num_task_tokens = int(input_info["num_task_tokens"])
        missing = num_task_tokens - context_without_global.shape[1]
        if missing < 0:
            raise ValueError("encoder contains more non-global tokens than the complete task sequence")
        mask_tokens = self.mask_token.expand(batch_size, missing, -1).to(dtype=context_tokens.dtype)
        context_with_mask = torch.cat([context_without_global, mask_tokens], dim=1)
        gather_restore = ids_restore.unsqueeze(-1).expand(-1, -1, context_with_mask.shape[-1])
        context_with_mask = torch.gather(context_with_mask, dim=1, index=gather_restore)
        context_with_mask = context_with_mask + self._context_embeddings(input_info, batch_size, context_tokens.device)

        task_info = input_info["tasks"][self.task]
        queries = context_with_mask[:, int(task_info["start_idx"]) : int(task_info["end_idx"])]
        gather_keep = ids_keep.unsqueeze(-1).expand(-1, -1, context_with_mask.shape[-1])
        decoder_context = torch.gather(context_with_mask, dim=1, index=gather_keep)
        if num_global:
            decoder_context = torch.cat([decoder_context, context_tokens[:, -num_global:]], dim=1)

        decoded = self.cross_attn(queries, decoder_context)
        for block in self.decoder_blocks:
            decoded = block(decoded)
        return self.out_proj(self.out_norm(decoded))


class MultiMAETargetsScMAE(nn.Module):
    """MultiMAE adapted to RNA expression, rank, and patch-statistic targets."""

    def __init__(
        self,
        num_genes: int,
        patch_size: int = 20,
        hidden_size: int = 128,
        depth: int = 3,
        num_heads: int = 4,
        decoder_size: int = 128,
        decoder_depth: int = 1,
        decoder_heads: int = 4,
        num_global_tokens: int = 1,
        dropout: float = 0.1,
    ) -> None:
        super().__init__()
        if num_genes <= 0 or patch_size <= 0:
            raise ValueError("num_genes and patch_size must be positive")
        if hidden_size % num_heads != 0:
            raise ValueError("hidden_size must be divisible by num_heads")
        self.num_genes = int(num_genes)
        self.patch_size = int(patch_size)
        self.num_patches = (self.num_genes + self.patch_size - 1) // self.patch_size
        self.num_global_tokens = int(num_global_tokens)
        task_dims = OrderedDict((("expr", self.patch_size), ("rank", self.patch_size), ("stat", 3)))
        self.input_adapters = nn.ModuleDict(
            {task: TaskInputAdapter(dim, hidden_size, dropout) for task, dim in task_dims.items()}
        )
        self.task_embeddings = nn.ParameterDict(
            {task: nn.Parameter(torch.randn(1, 1, hidden_size) * 0.02) for task in TASK_ORDER}
        )
        self.position = nn.Parameter(torch.randn(1, self.num_patches, hidden_size) * 0.02)
        self.global_tokens = nn.Parameter(torch.randn(1, self.num_global_tokens, hidden_size) * 0.02)
        self.encoder = nn.ModuleList([TransformerBlock(hidden_size, num_heads, dropout) for _ in range(depth)])
        self.encoder_norm = nn.LayerNorm(hidden_size)
        self.output_adapters = nn.ModuleDict(
            {
                task: TaskOutputDecoder(
                    task=task,
                    output_dim=dim,
                    num_patches=self.num_patches,
                    hidden_size=hidden_size,
                    decoder_size=decoder_size,
                    decoder_depth=decoder_depth,
                    decoder_heads=decoder_heads,
                    dropout=dropout,
                )
                for task, dim in task_dims.items()
            }
        )

    def sample_alphas(self, batch_size: int, n_tasks: int, alpha: float, device: torch.device, eps: float = 1e-5) -> torch.Tensor:
        valid_choices = torch.tensor(list(itertools.product([0, 1], repeat=n_tasks))[1:], dtype=torch.float32, device=device)
        chosen = torch.randint(0, valid_choices.shape[0], (batch_size,), device=device)
        return valid_choices.index_select(0, chosen) * float(alpha) + eps

    def generate_input_info(self, input_task_tokens: Mapping[str, torch.Tensor]) -> OrderedDict:
        info: OrderedDict[str, object] = OrderedDict()
        info["tasks"] = OrderedDict()
        start = 0
        for task, tokens in input_task_tokens.items():
            num_tokens = int(tokens.shape[1])
            info["tasks"][task] = {"num_tokens": num_tokens, "start_idx": start, "end_idx": start + num_tokens}
            start += num_tokens
        info["num_task_tokens"] = start
        info["num_global_tokens"] = self.num_global_tokens
        return info

    def generate_random_masks(
        self,
        input_tokens: Mapping[str, torch.Tensor],
        num_encoded_tokens: int,
        alpha: float,
        sample_tasks_uniformly: bool,
    ) -> tuple[dict[str, torch.Tensor], torch.Tensor, torch.Tensor]:
        batch_size = next(iter(input_tokens.values())).shape[0]
        device = next(iter(input_tokens.values())).device
        n_tasks = len(input_tokens)
        if not 0 < int(num_encoded_tokens) <= sum(tokens.shape[1] for tokens in input_tokens.values()):
            raise ValueError("num_encoded_tokens must be in the complete task-token range")
        if sample_tasks_uniformly:
            alphas = self.sample_alphas(batch_size, n_tasks, alpha, device)
            task_dist = Dirichlet(alphas).sample()
        else:
            task_dist = Dirichlet(torch.full((n_tasks,), float(alpha), device=device)).sample((batch_size,))
        samples_per_task = (task_dist * int(num_encoded_tokens)).round().long()

        task_masks = []
        token_counts = [tokens.shape[1] for tokens in input_tokens.values()]
        for task_idx, num_tokens in enumerate(token_counts):
            noise = torch.rand(batch_size, num_tokens, device=device)
            shuffled = torch.argsort(noise, dim=1)
            arange = torch.arange(num_tokens, device=device).unsqueeze(0).expand(batch_size, -1)
            ranked = torch.gather(arange, dim=1, index=shuffled)
            task_masks.append(torch.where(ranked < samples_per_task[:, task_idx].unsqueeze(1), 0, 1))

        mask_all = torch.cat(task_masks, dim=1)
        ids_shuffle = torch.argsort(mask_all + torch.rand_like(mask_all.float()), dim=1)
        ids_restore = torch.argsort(ids_shuffle, dim=1)
        ids_keep = ids_shuffle[:, : int(num_encoded_tokens)]
        adjusted_mask = torch.ones_like(mask_all)
        adjusted_mask[:, : int(num_encoded_tokens)] = 0
        adjusted_mask = torch.gather(adjusted_mask, dim=1, index=ids_restore)
        split_masks = torch.split(adjusted_mask, token_counts, dim=1)
        return {task: mask for task, mask in zip(input_tokens.keys(), split_masks)}, ids_keep, ids_restore

    def encode_task_inputs(self, task_inputs: Mapping[str, torch.Tensor]) -> OrderedDict:
        tokens: OrderedDict[str, torch.Tensor] = OrderedDict()
        for task in TASK_ORDER:
            if task not in task_inputs:
                raise ValueError(f"missing task input {task!r}")
            x = task_inputs[task]
            if x.ndim != 3 or x.shape[1] != self.num_patches:
                raise ValueError(f"{task} input must be [batch, {self.num_patches}, features], got {tuple(x.shape)}")
            tokens[task] = self.input_adapters[task](x)
            tokens[task] = tokens[task] + self.position.to(dtype=tokens[task].dtype, device=tokens[task].device)
            tokens[task] = tokens[task] + self.task_embeddings[task].to(dtype=tokens[task].dtype, device=tokens[task].device)
        return tokens

    def forward(
        self,
        task_inputs: Mapping[str, torch.Tensor],
        num_encoded_tokens: int,
        alpha: float = 1.0,
        sample_tasks_uniformly: bool = False,
    ) -> dict[str, torch.Tensor | dict[str, torch.Tensor]]:
        input_task_tokens = self.encode_task_inputs(task_inputs)
        input_info = self.generate_input_info(input_task_tokens)
        task_masks, ids_keep, ids_restore = self.generate_random_masks(
            input_task_tokens, num_encoded_tokens, alpha, sample_tasks_uniformly
        )
        all_tokens = torch.cat(list(input_task_tokens.values()), dim=1)
        gathered = torch.gather(all_tokens, dim=1, index=ids_keep.unsqueeze(-1).expand(-1, -1, all_tokens.shape[-1]))
        global_tokens = self.global_tokens.expand(gathered.shape[0], -1, -1).to(dtype=gathered.dtype, device=gathered.device)
        encoded = torch.cat([gathered, global_tokens], dim=1)
        for block in self.encoder:
            encoded = block(encoded)
        encoded = self.encoder_norm(encoded)
        preds = {
            task: decoder(encoded, input_info, ids_keep, ids_restore)
            for task, decoder in self.output_adapters.items()
        }
        embedding = encoded[:, -self.num_global_tokens :].mean(dim=1) if self.num_global_tokens else encoded.mean(dim=1)
        preds["embedding"] = embedding
        preds["task_masks"] = task_masks
        return preds

    @torch.no_grad()
    def feature_from_tasks(self, task_inputs: Mapping[str, torch.Tensor]) -> torch.Tensor:
        input_task_tokens = self.encode_task_inputs(task_inputs)
        all_tokens = torch.cat(list(input_task_tokens.values()), dim=1)
        global_tokens = self.global_tokens.expand(all_tokens.shape[0], -1, -1).to(dtype=all_tokens.dtype, device=all_tokens.device)
        encoded = torch.cat([all_tokens, global_tokens], dim=1)
        for block in self.encoder:
            encoded = block(encoded)
        encoded = self.encoder_norm(encoded)
        return encoded[:, -self.num_global_tokens :].mean(dim=1) if self.num_global_tokens else encoded.mean(dim=1)
