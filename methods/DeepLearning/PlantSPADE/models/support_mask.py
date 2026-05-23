import torch
import torch.nn as nn


# ---------- 支持掩码预测网络 ----------
# 核心问题：scRNA-seq中零值可能是"基因未表达"（生物学真零），
# 也可能是"技术原因导致的dropout"（假零）。该网络预测每个(cell,gene)对的激活概率。


class SupportMaskNet(nn.Module):
    """预测每个(cell,gene)对的基因激活概率，以区分真零与dropout零.

    作用：使重建误差仅聚焦于真正可能出错的位置（高激活概率但实际为零的基因），
    避免把"基因确实未表达"的零值当作错误来惩罚。
    """

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
        # 输出：每个基因在该细胞中是否激活的概率（sigmoid激活）
        self.head = nn.Linear(in_dim, num_genes)

    def forward(self, x: torch.Tensor) -> dict:
        h = self.encoder(x)
        logits = self.head(h)
        return {
            "gene_activation_logits": logits,
            "gene_activation_prob": torch.sigmoid(logits),
            "encoder_output": h,
        }

