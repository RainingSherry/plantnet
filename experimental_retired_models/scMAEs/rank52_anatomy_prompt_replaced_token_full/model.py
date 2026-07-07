from __future__ import annotations

import torch
from torch import nn


class AnatomyPromptScMAE(nn.Module):
    """scMAE with gene-module anatomy prompts and replaced-token heads."""

    def __init__(
        self,
        input_dim: int,
        hidden_size: int = 128,
        decoder_hidden: int = 128,
        n_modules: int = 32,
        n_token_bins: int = 8,
        dropout: float = 0.0,
    ):
        super().__init__()
        self.input_dim = int(input_dim)
        self.hidden_size = int(hidden_size)
        self.n_modules = int(n_modules)
        self.n_token_bins = int(n_token_bins)
        self.encoder_backbone = nn.Sequential(
            nn.Dropout(dropout),
            nn.Linear(input_dim, 256),
            nn.LayerNorm(256),
            nn.GELU(),
            nn.Linear(256, hidden_size),
            nn.LayerNorm(hidden_size),
            nn.GELU(),
        )
        self.module_prompt = nn.Embedding(n_modules, hidden_size)
        self.module_value_projector = nn.Sequential(
            nn.Linear(n_modules, hidden_size),
            nn.LayerNorm(hidden_size),
            nn.GELU(),
        )
        self.prompt_gate = nn.Sequential(nn.Linear(hidden_size * 2, hidden_size), nn.Sigmoid())
        self.prompt_norm = nn.LayerNorm(hidden_size)
        self.mask_predictor = nn.Linear(hidden_size, input_dim)
        self.decoder = nn.Sequential(
            nn.Linear(hidden_size + input_dim, decoder_hidden),
            nn.GELU(),
            nn.Linear(decoder_hidden, input_dim),
        )
        self.replaced_head = nn.Linear(hidden_size, input_dim)
        self.token_head = nn.Linear(hidden_size, input_dim * n_token_bins)

    def module_summary(self, x: torch.Tensor, module_ids: torch.Tensor) -> torch.Tensor:
        one_hot = torch.nn.functional.one_hot(module_ids.long(), num_classes=self.n_modules).float().to(x.device)
        denom = one_hot.sum(dim=0).clamp_min(1.0)
        return x.matmul(one_hot) / denom.view(1, -1)

    def encode(self, x: torch.Tensor, module_ids: torch.Tensor) -> torch.Tensor:
        base = self.encoder_backbone(x)
        module_values = self.module_value_projector(self.module_summary(x, module_ids))
        prompt_bank = self.module_prompt(module_ids.long().to(x.device)).mean(dim=0, keepdim=True).expand_as(base)
        gate = self.prompt_gate(torch.cat([base, module_values], dim=1))
        return self.prompt_norm(base + gate * (module_values + prompt_bank))

    def forward(self, x: torch.Tensor, module_ids: torch.Tensor) -> dict[str, torch.Tensor]:
        z = self.encode(x, module_ids)
        mask_logits = self.mask_predictor(z)
        recon = self.decoder(torch.cat([z, mask_logits], dim=1))
        token_logits = self.token_head(z).view(x.shape[0], self.input_dim, self.n_token_bins)
        return {
            "embedding": z,
            "reconstruction": recon,
            "mask_logits": mask_logits,
            "replaced_logits": self.replaced_head(z),
            "token_logits": token_logits,
        }
