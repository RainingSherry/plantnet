from __future__ import annotations

import numpy as np
import torch
import torch.nn as nn

from .axial_encoder import AxialEncoder
from .decoder import Decoder
from .mask_head import MaskHead
from .mlp_encoder import MLPEncoder


class CAAMStudent(nn.Module):
    def __init__(self, encoder: nn.Module, decoder: Decoder, mask_head: MaskHead, encoder_type: str) -> None:
        super().__init__()
        self.encoder = encoder
        self.decoder = decoder
        self.mask_head = mask_head
        self.encoder_type = encoder_type

    def forward(self, x_tilde: torch.Tensor, mask: torch.Tensor | None = None, indices: torch.Tensor | None = None) -> dict:
        enc = self.encoder(x_tilde, query_indices=indices) if self.encoder_type == "axial" else self.encoder(x_tilde)
        z = enc["z"]
        mask_logits = self.mask_head(z)
        x_hat = self.decoder(z, mask_logits, oracle_mask=mask)
        return {**enc, "mask_logits": mask_logits, "x_hat": x_hat}

    def feature(self, x: torch.Tensor, indices: torch.Tensor | None = None) -> torch.Tensor:
        enc = self.encoder(x, query_indices=indices) if self.encoder_type == "axial" else self.encoder(x)
        return enc["z"]

    def refresh_context_cache(self, context_x: torch.Tensor, context_indices: torch.Tensor) -> None:
        if self.encoder_type == "axial":
            self.encoder.refresh_context_cache(context_x, context_indices)

    def context_cache_checksum(self) -> float:
        if self.encoder_type == "axial":
            return self.encoder.context_cache_checksum()
        return 0.0


def build_student(
    *,
    n_genes: int,
    config: dict,
    assignment: np.ndarray | None = None,
) -> CAAMStudent:
    latent_dim = int(config["model"]["latent_dim"])
    encoder_type = config["model"]["encoder_type"]
    if encoder_type == "mlp":
        encoder = MLPEncoder(
            n_genes=n_genes,
            latent_dim=latent_dim,
            hidden_dim=int(config["model"].get("mlp_hidden_dim", 256)),
            dropout=float(config["model"].get("dropout", 0.0)),
        )
    elif encoder_type == "axial":
        if assignment is None:
            raise ValueError("Axial encoder requires a gene module assignment matrix.")
        encoder = AxialEncoder(
            assignment=assignment,
            token_dim=int(config["axial"]["token_dim"]),
            latent_dim=latent_dim,
            gene_attention_heads=int(config["axial"]["gene_attention_heads"]),
            gene_attention_layers=int(config["axial"]["gene_attention_layers"]),
            attention_dropout=float(config["axial"]["attention_dropout"]),
        )
    else:
        raise ValueError(f"Unknown encoder_type={encoder_type!r}")
    mask_head = MaskHead(latent_dim, n_genes)
    decoder = Decoder(latent_dim, n_genes, conditioning=config["model"]["decoder_mask_conditioning"])
    return CAAMStudent(encoder=encoder, decoder=decoder, mask_head=mask_head, encoder_type=encoder_type)

