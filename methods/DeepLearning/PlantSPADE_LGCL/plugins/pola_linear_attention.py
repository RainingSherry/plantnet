"""
Polarity-aware Contrastive Loss (adapted from PolaFormer, ICLR 2025)
Paper: PolaFormer: Polarity-aware Linear Attention for Vision Transformers
Reference: https://arxiv.org/abs/2501.15061

The core idea: decompose Q/K into positive (ReLU) and negative (ReLU(-)) branches,
modeling both gene co-expression and anti-correlation patterns in cell embeddings.
"""
import torch
import torch.nn as nn
import torch.nn.functional as F


class PolaContrastiveLoss(nn.Module):
    """Polarity-aware contrastive loss for cell embedding alignment.

    Applies polarity decomposition to the standard cosine similarity, separating
    positive (ReLU(q)) and negative (ReLU(-q)) signal branches to capture
    both gene co-expression and anti-correlation patterns.
    """

    def __init__(self, dim: int, num_heads: int = 4, alpha: float = 4.0):
        super().__init__()
        assert dim % num_heads == 0
        self.dim = dim
        self.num_heads = num_heads
        self.head_dim = dim // num_heads
        self.alpha = alpha
        # Learnable power per head
        self.power = nn.Parameter(torch.zeros(num_heads))

    def forward(self, local: torch.Tensor, global_view: torch.Tensor, temperature: float) -> torch.Tensor:
        B = local.shape[0]
        if B <= 2:
            return torch.zeros((), dtype=local.dtype, device=local.device)

        # Multi-head reshape: (B, dim) -> (B, num_heads, head_dim)
        q = local.reshape(B, self.num_heads, self.head_dim)
        k = global_view.reshape(B, self.num_heads, self.head_dim)

        # Cosine similarity per head: (B, num_heads)
        q_norm = F.normalize(q, dim=-1)
        k_norm = F.normalize(k, dim=-1)
        cos_sim = (q_norm * k_norm).sum(dim=-1)  # (B, num_heads)

        # Learnable power function (PolaFormer)
        power = 1.0 + self.alpha * torch.sigmoid(self.power)  # (num_heads,)
        power = power.unsqueeze(0)  # (1, num_heads)

        # Positive polarity: ReLU(q) vs ReLU(k)
        cos_pos = F.relu(cos_sim)
        # Negative polarity: ReLU(-q) vs ReLU(-k)
        cos_neg = F.relu(-cos_sim)

        # Apply learnable power to separate branches
        polarity_score = (cos_pos ** power).mean(dim=-1) - (cos_neg ** power).mean(dim=-1)  # (B,)

        # N-way contrastive loss: use polarity-weighted cosine similarity
        # Standard cosine similarity for logits
        base_sim = F.normalize(local, dim=1) @ F.normalize(global_view, dim=1).T / temperature  # (B, B)
        # Modulate by polarity score (outer product weighting)
        pol_weight = polarity_score.unsqueeze(0) * polarity_score.unsqueeze(1)  # (B, B)
        logits = base_sim + 0.1 * pol_weight  # modulate with polarity

        labels = torch.arange(B, device=local.device)
        return F.cross_entropy(logits, labels)
