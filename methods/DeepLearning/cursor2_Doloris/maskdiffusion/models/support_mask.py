import torch
import torch.nn as nn
import torch.nn.functional as F


class SupportMaskNet(nn.Module):
    """Predicts per-gene activation probability for each cell.

    This module learns the *expression support structure*: for each cell c
    and gene g, it predicts P(X_{cg} > 0) based on the expression vector.

    The predicted probabilities are used as per-gene importance weights in the
    reconstruction loss, following the DOLORIS sparsity-masking philosophy:
      - Genes with high mask probability + observed expression -> full loss weight
      - Genes with low mask probability or not expressed -> reduced weight

    Architecture:
        X (n_genes)
          -> Linear(n_genes, hidden_dims[0]) + LayerNorm + GELU + Dropout
          -> [hidden_dims layers]
          -> Linear(hidden_dims[-1], n_genes)   (gene activation logits)
          -> Sigmoid -> gene_activation_prob
    """

    def __init__(
        self,
        n_genes: int,
        hidden_dims: list = None,
        dropout: float = 0.1,
    ):
        super().__init__()
        if hidden_dims is None:
            hidden_dims = [512, 256, 128]

        layers = []
        in_dim = n_genes
        for h_dim in hidden_dims:
            layers.extend([
                nn.Linear(in_dim, h_dim),
                nn.LayerNorm(h_dim),
                nn.GELU(),
                nn.Dropout(dropout),
            ])
            in_dim = h_dim

        self.backbone = nn.Sequential(*layers)
        self.gene_activation_head = nn.Linear(in_dim, n_genes)

    def forward(self, x: torch.Tensor) -> dict:
        """
        Args:
            x: Input expression vector, shape (batch, n_genes).
        Returns:
            dict with keys:
                - gene_activation_prob: (batch, n_genes) sigmoid probabilities
                - gene_activation_logits: (batch, n_genes) raw logits
        """
        h = self.backbone(x)
        logits = self.gene_activation_head(h)
        probs = torch.sigmoid(logits)
        return {
            "gene_activation_prob": probs,
            "gene_activation_logits": logits,
        }

    def get_gene_activation_loss(
        self,
        x: torch.Tensor,
        support: torch.Tensor = None,
        threshold: float = 0.5,
    ) -> torch.Tensor:
        """Binary cross-entropy loss for support prediction.

        Args:
            x: Input expression (batch, n_genes).
            support: Binary support target (batch, n_genes), values in {0, 1}.
                     If None, inferred as (x > 0).float().
            threshold: Ignored (kept for API compatibility).
        Returns:
            Scalar BCE loss (mean over all positions).
        """
        if support is None:
            support = (x > 0).float()
        output = self.forward(x)
        # BCE with logits = -log(sigmoid) for target=1 - log(1-sigmoid) for target=0
        loss = F.binary_cross_entropy_with_logits(
            output["gene_activation_logits"],
            support,
            reduction="mean",
        )
        return loss

    @torch.no_grad()
    def predict_support(self, x: torch.Tensor, threshold: float = 0.5) -> torch.Tensor:
        """Hard support prediction from expression input.

        Args:
            x: Input expression (batch, n_genes).
            threshold: Threshold for binarizing probabilities.
        Returns:
            Binary support tensor (batch, n_genes), values in {0, 1}.
        """
        probs = self.forward(x)["gene_activation_prob"]
        return (probs >= threshold).float()
