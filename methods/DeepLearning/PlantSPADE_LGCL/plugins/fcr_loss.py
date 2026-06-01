"""
Frequency Contrastive Regularization (FCR)
Paper: Efficient Frequency-Domain Image Deraining with Contrastive Regularization (ECCV 2024)
Reference: https://www.ecva.net/papers/eccv_2024/papers_ECCV/papers/05751.pdf

Adjusted for single-cell RNA-seq: applies FFT-based contrastive loss over cell embeddings.
The loss is logarithmically compressed to stay in a stable range (~0.1-10) compatible with BPR.
"""
import torch
import torch.nn as nn


class FCRLoss(nn.Module):
    """Frequency Contrastive Regularization for cell embeddings.

    Applies FFT to anchor/positive/negative embeddings, computes L1-based
    contrastive loss in the frequency domain. The ratio is logarithmically
    compressed to prevent extreme values from unit-norm embeddings.
    """

    def __init__(self, multi_n_num: int = 2, seq_len: int = 1, emb_dim: int = 32):
        super().__init__()
        self.multi_n_num = multi_n_num
        self.seq_len = seq_len
        self.emb_dim = emb_dim

    def forward(self, anchor: torch.Tensor, positive: torch.Tensor, negative: torch.Tensor) -> torch.Tensor:
        b = anchor.shape[0]
        if b <= self.multi_n_num:
            return torch.zeros((), dtype=anchor.dtype, device=anchor.device)

        # Reshape to (B, 1, emb_dim, seq_len) for 2D FFT
        a = anchor.view(b, 1, self.emb_dim, self.seq_len)
        p = positive.view(b, 1, self.emb_dim, self.seq_len)
        n = negative.view(b, 1, self.emb_dim, self.seq_len)

        # FFT magnitude
        a_fft = torch.fft.fft2(a).abs()
        p_fft = torch.fft.fft2(p).abs()
        n_fft = torch.fft.fft2(n).abs()

        # Compute contrastive ratio: d_ap / d_an
        # For normalized embeddings, raw values are ~5-10; log1p compresses to ~1.8-2.4
        ratio_sum = 0.0
        multi_n = self.multi_n_num
        for i in range(b):
            d_ap = torch.abs(a_fft[i] - p_fft[i]).mean()
            # Negative sampling: pick one or more negatives
            for j in range(i + 1, min(i + 1 + multi_n, b)):
                d_an = torch.abs(a_fft[i] - n_fft[j]).mean()
                # Apply log1p to compress the ratio to prevent extreme values
                ratio_sum += torch.log1p(d_ap) / (torch.log1p(d_an) + 1e-7)

        ratio_sum /= (multi_n * b)
        return ratio_sum
