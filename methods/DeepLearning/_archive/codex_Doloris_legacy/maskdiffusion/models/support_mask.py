import torch
import torch.nn as nn


class SupportMaskNet(nn.Module):
    """Predict observed expression support probabilities for each cell-gene pair."""

    def __init__(self, num_genes: int, hidden_dims=None, dropout: float = 0.1):
        super().__init__()
        if hidden_dims is None:
            hidden_dims = [512, 256, 128]

        layers = []
        in_dim = num_genes
        for hidden_dim in hidden_dims:
            layers.extend(
                [
                    nn.Dropout(dropout),
                    nn.Linear(in_dim, hidden_dim),
                    nn.LayerNorm(hidden_dim),
                    nn.Mish(inplace=True),
                ]
            )
            in_dim = hidden_dim
        self.encoder = nn.Sequential(*layers)
        self.head = nn.Linear(in_dim, num_genes)

    def forward(self, x: torch.Tensor) -> dict:
        h = self.encoder(x)
        logits = self.head(h)
        return {
            "gene_activation_logits": logits,
            "gene_activation_prob": torch.sigmoid(logits),
            "encoder_output": h,
        }

