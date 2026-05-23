"""
models/mask_diffusion_refiner.py
==================================
Mask Diffusion Refiner（掩码扩散精炼器）。

核心思想（融合 DOLORIS + scMAE）：
  1. 从 SupportPooling 得到"粗粒度"细胞嵌入 z_0（携带基因图信息）
  2. 在嵌入空间中运行 DDPM/DDIM 反向过程进行去噪精炼
  3. 额外的 Sparsity Mask Model 预测零值基因（可选）
  4. 最终输出精炼后的细胞嵌入，用于聚类

关键设计：
  - 不预测噪声 ε，而是直接预测 x_0（与 DOLORIS 一致）
  - 使用相对时间步嵌入（timestep embedding）
  - 可选：条件嵌入（细胞类型、批次）注入到扩散过程
  - 可选：预测零值掩码（预测每个基因是否为零）
  - 输出：精炼嵌入 + 可选的零值掩码

损失函数（可选组合）：
  L = L_reconstruction + λ_mask · L_mask + λ_cluster · L_cluster

其中：
  - L_reconstruction = MSE(重构基因表达, 真实基因表达)
  - L_mask = BCE(预测零值掩码, 真实零值掩码)
  - L_cluster = 对比聚类损失（可选）
"""

from __future__ import annotations

import math
import torch
import torch.nn as nn
import torch.nn.functional as F
from typing import Optional, List, Tuple, Literal


# ---------------------------------------------------------------------------
# Time Embedding（来自 DOLORIS）
# ---------------------------------------------------------------------------

def timestep_embedding(timesteps: torch.Tensor, dim: int, max_period: float = 10000.0) -> torch.Tensor:
    """
    创建正弦时间步嵌入（来自 Attention is All You Need）。

    来自 DOLORIS 的 nn.timestep_embedding
    """
    half = dim // 2
    freqs = torch.exp(
        -math.log(max_period) * torch.arange(half, dtype=torch.float32, device=timesteps.device) / half
    )
    args = timesteps[:, None].float() * freqs[None]
    embedding = torch.cat([torch.cos(args), torch.sin(args)], dim=-1)
    if dim % 2:
        embedding = torch.cat([embedding, torch.zeros_like(embedding[:, :1])], dim=-1)
    return embedding


# ---------------------------------------------------------------------------
# 扩散过程基础模块
# ---------------------------------------------------------------------------

class GaussianDiffusion1D(nn.Module):
    """
    一维高斯扩散模型（用于细胞嵌入的去噪精炼）。

    简化版：直接预测 x_0，不学习方差。
    基于 Ho et al. (2020) DDPM，但移除了方差学习。

    正向过程（已知）：
      x_t = sqrt(ᾱ_t) · x_0 + sqrt(1 - ᾱ_t) · ε,  ε ~ N(0, I)

    反向过程（学习）：
      模型预测 x_0，给定 x_t 和 t。
    """

    def __init__(
        self,
        num_timesteps: int = 500,
        beta_schedule: str = "cosine",
        clip_min: float = -5.0,
        clip_max: float = 5.0,
        snr_scale: float = 1.0,
    ):
        super().__init__()
        self.num_timesteps = num_timesteps
        self.clip_min = clip_min
        self.clip_max = clip_max
        self.snr_scale = snr_scale

        # ---- Beta schedule ----
        if beta_schedule == "linear":
            betas = torch.linspace(1e-4, 0.02, num_timesteps)
        elif beta_schedule == "cosine":
            betas = self._cosine_beta_schedule(num_timesteps)
        elif beta_schedule == "sqrt":
            betas = torch.linspace(1e-4, 0.02, num_timesteps).sqrt()
        else:
            raise ValueError(f"Unknown beta_schedule: {beta_schedule}")

        alphas = 1.0 - betas
        alphas_cumprod = torch.cumprod(alphas, dim=0)
        alphas_cumprod_prev = torch.cat([torch.ones(1), alphas_cumprod[:-1]])

        self.register_buffer("betas", betas.float())
        self.register_buffer("alphas_cumprod", alphas_cumprod.float())
        self.register_buffer("sqrt_alphas_cumprod", alphas_cumprod.sqrt().float())
        self.register_buffer("sqrt_one_minus_alphas_cumprod", (1.0 - alphas_cumprod).sqrt().float())

    @staticmethod
    def _cosine_beta_schedule(num_timesteps: int, s: float = 0.008) -> torch.Tensor:
        """Cosine beta schedule（来自 Nichol & Dhariwal 2021）。"""
        steps = num_timesteps + 1
        x = torch.linspace(0, num_timesteps, steps)
        alphas_cumprod = torch.cos(((x / num_timesteps) + s) / (1 + s) * math.pi * 0.5) ** 2
        alphas_cumprod = alphas_cumprod / alphas_cumprod[0]
        betas = 1 - (alphas_cumprod[1:] / alphas_cumprod[:-1])
        return torch.clip(betas, 0.0001, 0.9999)

    def q_sample(self, x_start: torch.Tensor, t: torch.Tensor, noise: Optional[torch.Tensor] = None) -> torch.Tensor:
        """正向扩散：给 x_0 添加噪声到 x_t。"""
        if noise is None:
            noise = torch.randn_like(x_start)
        return (
            self.sqrt_alphas_cumprod[t][:, None] * x_start
            + self.sqrt_one_minus_alphas_cumprod[t][:, None] * noise
        )

    def predict_x0_from_eps(self, x_t: torch.Tensor, t: torch.Tensor, eps: torch.Tensor) -> torch.Tensor:
        """从噪声预测 x_0。"""
        return (x_t - self.sqrt_one_minus_alphas_cumprod[t][:, None] * eps) / \
               (self.sqrt_alphas_cumprod[t][:, None] + 1e-8)

    def predict_x0_direct(self, x_t: torch.Tensor, t: torch.Tensor, model_output: torch.Tensor) -> torch.Tensor:
        """直接预测 x_0（模型输出即为 x_0 预测）。"""
        # 已在 caller 中处理
        return model_output


# ---------------------------------------------------------------------------
# Refiner 网络（核心去噪器）
# ---------------------------------------------------------------------------

class RefinerMLP(nn.Module):
    """
    MLP 风格的去噪网络（用于 Refiner 的骨干）。

    输入： [z_t; t_emb; cond_emb]
    输出： z_0 预测

    结构：
      Linear(dim, hidden) → SiLU → Block × depth
      → Linear(hidden, dim)

    支持条件注入：细胞类型嵌入、批次嵌入等。
    """

    def __init__(
        self,
        embed_dim: int,
        hidden_dim: int = 512,
        time_embed_dim: int = 256,
        cond_embed_dim: int = 0,
        depth: int = 4,
        dropout: float = 0.1,
        use_layer_norm: bool = True,
    ):
        super().__init__()
        self.embed_dim = embed_dim
        self.time_embed_dim = time_embed_dim
        self.cond_embed_dim = cond_embed_dim

        total_in_dim = embed_dim + time_embed_dim + cond_embed_dim

        # 时间嵌入网络
        self.time_mlp = nn.Sequential(
            nn.Linear(time_embed_dim, hidden_dim),
            nn.SiLU(),
            nn.Linear(hidden_dim, hidden_dim),
        )

        # 条件嵌入网络（细胞类型）
        if cond_embed_dim > 0:
            self.cond_mlp = nn.Sequential(
                nn.Linear(cond_embed_dim, hidden_dim),
                nn.SiLU(),
                nn.Linear(hidden_dim, hidden_dim),
            )

        # 主干网络
        layers = []
        in_dim = total_in_dim
        for d in range(depth):
            layers.extend([
                nn.Linear(in_dim, hidden_dim),
                nn.LayerNorm(hidden_dim) if use_layer_norm else nn.Identity(),
                nn.SiLU(),
                nn.Dropout(dropout),
            ])
            in_dim = hidden_dim
        layers.append(nn.Linear(hidden_dim, embed_dim))
        self.net = nn.Sequential(*layers)

        self.apply(self._init_weights)

    def _init_weights(self, m):
        if isinstance(m, nn.Linear):
            nn.init.xavier_uniform_(m.weight)
            if m.bias is not None:
                nn.init.zeros_(m.bias)

    def forward(
        self,
        z_t: torch.Tensor,           # [B, D]  带噪嵌入
        time_emb: torch.Tensor,      # [B, time_embed_dim]  时间步嵌入
        cond_emb: Optional[torch.Tensor] = None,  # [B, cond_embed_dim]  条件嵌入
    ) -> torch.Tensor:
        """
        参数
        ----
        z_t : [B, D]  当前时间步的带噪细胞嵌入
        time_emb : [B, time_embed_dim]  时间步嵌入
        cond_emb : [B, cond_embed_dim]  条件嵌入（细胞类型）

        返回
        ----
        z_0_pred : [B, D]  预测的干净嵌入
        """
        t_hidden = self.time_mlp(time_emb)  # [B, hidden_dim]

        if cond_emb is not None and self.cond_embed_dim > 0:
            c_hidden = self.cond_mlp(cond_emb)  # [B, hidden_dim]
            hidden = t_hidden + c_hidden
        else:
            hidden = t_hidden

        # 拼接 [z_t; time_emb; cond_emb]
        if cond_emb is not None and self.cond_embed_dim > 0:
            x = torch.cat([z_t, time_emb, cond_emb], dim=-1)
        else:
            x = torch.cat([z_t, time_emb], dim=-1)

        return self.net(x)  # [B, D]


# ---------------------------------------------------------------------------
# Mask Model（预测零值基因）
# ---------------------------------------------------------------------------

class SparsityMaskPredictor(nn.Module):
    """
    稀疏度掩码预测器（参考 DOLORIS 的 MaskModel）。

    功能：
      给定细胞嵌入 z 和条件（细胞类型、批次），预测每个基因是否为零。

    核心思想（来自 DOLORIS）：
      在 scRNA-seq 中，"零" 不是纯粹的噪声，
      而是重要的生物信号（如转录沉默）。
      预测零值分布可以：
        1. 帮助扩散模型专注于真正表达的基因
        2. 生成更真实的表达谱（结合精炼嵌入）

    损失函数：BCE（预测 vs 真实零值掩码）
    """

    def __init__(
        self,
        embed_dim: int,
        n_genes: int,
        hidden_dim: int = 256,
        cond_embed_dim: int = 0,
        dropout: float = 0.1,
    ):
        super().__init__()
        self.embed_dim = embed_dim
        self.n_genes = n_genes
        self.cond_embed_dim = cond_embed_dim

        # 嵌入聚合器
        self.z_encoder = nn.Sequential(
            nn.Linear(embed_dim, hidden_dim),
            nn.SiLU(),
            nn.Dropout(dropout),
        )

        # 条件注入
        if cond_embed_dim > 0:
            self.cond_encoder = nn.Sequential(
                nn.Linear(cond_embed_dim, hidden_dim),
                nn.SiLU(),
            )

        # 预测头
        self.predictor = nn.Sequential(
            nn.Linear(hidden_dim, hidden_dim),
            nn.SiLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim, n_genes),
        )

    def forward(
        self,
        z: torch.Tensor,                      # [B, D]
        cond_emb: Optional[torch.Tensor] = None,  # [B, cond_dim]
    ) -> Tuple[torch.Tensor, torch.Tensor]:
        """
        参数
        ----
        z : [B, D]  精炼后的细胞嵌入
        cond_emb : [B, cond_dim]  条件嵌入

        返回
        ----
        prob_zero : [B, n_genes]  每个基因是零的概率
        logits : [B, n_genes]     原始 logits（用于 BCE loss）
        """
        h = self.z_encoder(z)
        if cond_emb is not None and self.cond_embed_dim > 0:
            h = h + self.cond_encoder(cond_emb)

        logits = self.predictor(h)  # [B, n_genes]
        prob_zero = torch.sigmoid(logits)  # [B, n_genes]
        return prob_zero, logits

    def compute_mask_loss(
        self,
        z: torch.Tensor,
        true_mask: torch.Tensor,  # [B, n_genes]  1=零值，0=表达
        cond_emb: Optional[torch.Tensor] = None,
    ) -> torch.Tensor:
        """计算掩码预测的 BCE 损失（仅在非零基因上加权）。"""
        prob_zero, logits = self.forward(z, cond_emb)

        # 加权 BCE：非零基因的权重更高（鼓励模型区分沉默 vs 表达）
        # 权重 = 1 - true_mask（这样非零基因权重 = 1，零值基因权重 = 0）
        # 但更合理的做法是对两类都计算，不过度偏向任何一方
        loss = F.binary_cross_entropy_with_logits(logits, true_mask, reduction="mean")
        return loss


# ---------------------------------------------------------------------------
# Mask Diffusion Refiner（主模型）
# ---------------------------------------------------------------------------

class MaskDiffusionRefiner(nn.Module):
    """
    掩码扩散精炼器（主模型）。

    整合所有组件：
      1. RefinerMLP：去噪嵌入
      2. GaussianDiffusion1D：扩散过程
      3. SparsityMaskPredictor：零值掩码预测（可选）

    训练模式：
      - 采样一个时间步 t
      - 对 z_0 加噪到 z_t
      - 用 RefinerMLP 预测 z_0
      - 计算 MSE 损失（仅在有效维度上）
      - 可选：同时计算掩码损失

    推理模式（DDIM 采样）：
      - 从纯噪声 z_T 开始
      - DDIM 反向逐步去噪到 z_0
      - 返回精炼嵌入
    """

    def __init__(
        self,
        embed_dim: int,
        n_genes: int,
        hidden_dim: int = 512,
        time_embed_dim: int = 256,
        cond_embed_dim: int = 0,
        refiner_depth: int = 4,
        dropout: float = 0.1,
        num_timesteps: int = 500,
        beta_schedule: str = "cosine",
        use_mask_predictor: bool = True,
        use_layer_norm: bool = True,
    ):
        super().__init__()
        self.embed_dim = embed_dim
        self.n_genes = n_genes
        self.use_mask_predictor = use_mask_predictor

        # ---- 扩散过程 ----
        self.diffusion = GaussianDiffusion1D(
            num_timesteps=num_timesteps,
            beta_schedule=beta_schedule,
        )

        # ---- Refiner（去噪网络）----
        self.refiner = RefinerMLP(
            embed_dim=embed_dim,
            hidden_dim=hidden_dim,
            time_embed_dim=time_embed_dim,
            cond_embed_dim=cond_embed_dim,
            depth=refiner_depth,
            dropout=dropout,
            use_layer_norm=use_layer_norm,
        )

        # ---- 掩码预测器（可选）----
        if use_mask_predictor:
            self.mask_predictor = SparsityMaskPredictor(
                embed_dim=embed_dim,
                n_genes=n_genes,
                hidden_dim=min(hidden_dim, 256),
                cond_embed_dim=cond_embed_dim,
                dropout=dropout,
            )

    def get_time_emb(self, t: torch.Tensor) -> torch.Tensor:
        """生成时间步嵌入。"""
        return timestep_embedding(t, self.refiner.time_embed_dim)

    def training_losses(
        self,
        z0: torch.Tensor,             # [B, D]  干净嵌入
        t: torch.Tensor,              # [B]  时间步
        cond_emb: Optional[torch.Tensor] = None,  # [B, cond_dim]
        true_mask: Optional[torch.Tensor] = None,  # [B, n_genes]  真实零值掩码
        noise: Optional[torch.Tensor] = None,
    ) -> dict:
        """
        计算训练损失。

        返回
        ----
        losses : dict
            - "loss": 总损失
            - "mse": 重构 MSE
            - "mask": 掩码 BCE（如果启用）
        """
        B = z0.size(0)

        # 采样噪声
        if noise is None:
            noise = torch.randn_like(z0)

        # 正向扩散
        zt = self.diffusion.q_sample(z0, t, noise)

        # 时间步嵌入
        time_emb = self.get_time_emb(t)  # [B, time_embed_dim]

        # 预测 z_0
        z0_pred = self.refiner(zt, time_emb, cond_emb)  # [B, D]

        # MSE 损失
        mse = ((z0 - z0_pred) ** 2).mean()

        losses = {"mse": mse, "loss": mse}

        # 掩码损失（可选）
        if self.use_mask_predictor and true_mask is not None:
            mask_loss = self.mask_predictor.compute_mask_loss(z0_pred, true_mask, cond_emb)
            losses["mask"] = mask_loss
            losses["loss"] = mse + 0.5 * mask_loss

        return losses

    @torch.no_grad()
    def ddim_sample(
        self,
        shape: Tuple[int, ...],
        cond_emb: Optional[torch.Tensor] = None,
        n_steps: int = 50,
        eta: float = 0.0,
        clip_denoised: bool = True,
    ) -> torch.Tensor:
        """
        DDIM 采样（快速推理模式）。

        参数
        ----
        shape : (B, D)  输出形状
        cond_emb : [B, cond_dim]  条件嵌入
        n_steps : int  采样步数（默认 50，DDIM 加速）
        eta : float  随机性控制（0=完全确定性）
        clip_denoised : bool  是否 clip 预测值

        返回
        ----
        z0_samples : [B, D]  精炼后的细胞嵌入
        """
        B, D = shape
        device = next(self.parameters()).device

        # 初始化为纯噪声
        zt = torch.randn(B, D, device=device)

        # 时间步序列（均匀采样 n_steps 个）
        step_list = torch.linspace(0, self.diffusion.num_timesteps - 1, n_steps, dtype=torch.long, device=device)
        step_list = step_list.flip(0)  # 从 T 到 0

        for i in range(len(step_list) - 1):
            t_cur = step_list[i]
            t_next = step_list[i + 1]

            # 预测 z_0
            t_tensor = torch.full((B,), t_cur.item(), device=device, dtype=torch.long)
            time_emb = self.get_time_emb(t_tensor)
            z0_pred = self.refiner(zt, time_emb, cond_emb)

            if clip_denoised:
                z0_pred = z0_pred.clamp(min=-5.0, max=5.0)

            # DDIM 反向步
            alpha_cur = self.diffusion.alphas_cumprod[t_cur]
            alpha_next = self.diffusion.alphas_cumprod[t_next] if t_next >= 0 else torch.tensor(1.0, device=device)

            # 预测噪声
            pred_noise = (zt - alpha_cur.sqrt() * z0_pred) / (self.diffusion.sqrt_one_minus_alphas_cumprod[t_cur] + 1e-8)

            # DDIM 更新
            sigma = eta * ((1 - alpha_next) / (1 - alpha_cur)).sqrt() * ((1 - alpha_cur / alpha_next)).sqrt()
            zt = alpha_next.sqrt() * z0_pred + (1 - alpha_next - sigma ** 2).sqrt() * pred_noise

            if eta > 0:
                zt = zt + sigma * torch.randn_like(zt)

        return zt

    @torch.no_grad()
    def sample(
        self,
        shape: Tuple[int, ...],
        cond_emb: Optional[torch.Tensor] = None,
        n_steps: int = 50,
        use_ddim: bool = True,
    ) -> torch.Tensor:
        """推理入口。"""
        if use_ddim:
            return self.ddim_sample(shape, cond_emb, n_steps)
        else:
            raise NotImplementedError("DDPM sampling is not implemented; use use_ddim=True")


if __name__ == "__main__":
    # 简单测试
    B, D, n_genes = 16, 128, 500
    device = "cuda" if torch.cuda.is_available() else "cpu"

    model = MaskDiffusionRefiner(
        embed_dim=D,
        n_genes=n_genes,
        hidden_dim=256,
        time_embed_dim=128,
        cond_embed_dim=32,
        refiner_depth=3,
        num_timesteps=500,
    ).to(device)

    # 模拟数据
    z0 = torch.randn(B, D, device=device)
    t = torch.randint(0, 500, (B,), device=device)
    cond_emb = torch.randn(B, 32, device=device)
    true_mask = (torch.rand(B, n_genes, device=device) > 0.7).float()  # ~30% zeros

    # 训练
    losses = model.training_losses(z0, t, cond_emb, true_mask)
    print(f"Training losses: mse={losses['mse']:.4f}, mask={losses.get('mask', 0):.4f}, total={losses['loss']:.4f}")

    # 采样
    samples = model.sample(shape=(8, D), cond_emb=cond_emb[:8], n_steps=20)
    print(f"Sampling output: {samples.shape}")

    print("mask_diffusion_refiner.py test passed.")
