import math

import torch
import torch.nn as nn
import torch.nn.functional as F


# ---------- 时间步编码（正弦位置编码）----------
# 将扩散时间步 t ∈ [0, T] 映射为周期特征向量，让去噪网络感知当前处于哪个噪声水平


def timestep_embedding(timesteps: torch.Tensor, dim: int) -> torch.Tensor:
    """正弦位置编码：将扩散时间步映射为周期性特征，用于去噪网络的时序条件."""
    half = dim // 2
    freqs = torch.exp(
        -math.log(10000) * torch.arange(half, device=timesteps.device, dtype=torch.float32) / max(half - 1, 1)
    )
    args = timesteps[:, None].float() * freqs[None, :]
    emb = torch.cat([torch.sin(args), torch.cos(args)], dim=-1)
    if dim % 2:
        emb = F.pad(emb, (0, 1))
    return emb


# ---------- 扩散残差块 ----------
# DDPM中每个去噪步骤由一个带时间条件的大门控残差块实现


class DiffusionBlock(nn.Module):
    """带时间条件的大门控残差块：时间信息通过可学习的缩放因子注入特征流."""
    def __init__(self, hidden_dim: int, time_dim: int, dropout: float):
        super().__init__()
        self.time_proj = nn.Sequential(
            # 正弦时间编码 → hidden_dim
            nn.Linear(time_dim, hidden_dim),
            nn.Mish(inplace=True),
            nn.Linear(hidden_dim, hidden_dim),
        )
        # 门控残差：shortcut由学习到的门控比例调制
        self.net = nn.Sequential(
            nn.Linear(hidden_dim, hidden_dim * 2),
            nn.LayerNorm(hidden_dim * 2),
            nn.Mish(inplace=True),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim * 2, hidden_dim),
        )
        self.time_scale = nn.Parameter(torch.tensor(0.1))

    def forward(self, x: torch.Tensor, time_emb: torch.Tensor) -> torch.Tensor:
        h = x + self.time_scale * self.time_proj(time_emb)
        return h + self.net(h)


# ---------- U-Net风格去噪器 ----------
# 核心思想：从含噪隐变量 z_t 预测干净隐变量 z_0


class LatentDenoiser(nn.Module):
    """U-Net风格去噪器：从带噪隐变量 z_t 预测干净隐变量 z_0."""
    def __init__(
        self,
        latent_dim: int,
        hidden_dim: int = 256,
        num_layers: int = 3,
        time_dim: int = 128,
        dropout: float = 0.1,
    ):
        super().__init__()
        self.time_dim = time_dim
        self.time_mlp = nn.Sequential(
            # 正弦时间编码 → 展平特征
            nn.Linear(time_dim, hidden_dim),
            nn.Mish(inplace=True),
            nn.Linear(hidden_dim, time_dim),
        )
        self.input_proj = nn.Linear(latent_dim, hidden_dim)
        # 多层扩散块逐步去噪，每层都注入当前时间步信息
        self.blocks = nn.ModuleList([DiffusionBlock(hidden_dim, time_dim, dropout) for _ in range(num_layers)])
        self.output_proj = nn.Linear(hidden_dim, latent_dim)

    def forward(self, z_t: torch.Tensor, timesteps: torch.Tensor) -> torch.Tensor:
        t_emb = timestep_embedding(timesteps, self.time_dim)
        t_emb = self.time_mlp(t_emb)
        h = self.input_proj(z_t)
        for block in self.blocks:
            h = block(h, t_emb)
        return self.output_proj(h)


# ---------- 带扩散的稀疏重建自编码器 ----------
# 核心：编码 → 扩散隐空间 → 去噪 → 解码重建


class LatentDiffusionAE(nn.Module):
    """稀疏重建自编码器 + DDPM式隐空间去噪扩散."""
    def __init__(
        self,
        num_genes: int,
        latent_dim: int = 32,
        hidden_dims=None,
        diffusion_steps: int = 100,
        dropout: float = 0.1,
    ):
        super().__init__()
        if hidden_dims is None:
            hidden_dims = [512, 256]

        self.num_genes = num_genes
        self.latent_dim = latent_dim
        self.diffusion_steps = diffusion_steps

        # 编码器：基因表达向量 → 低维隐变量 z
        enc_layers = []
        in_dim = num_genes
        for hidden_dim in hidden_dims:
            enc_layers.extend(
                [
                    nn.Dropout(dropout),
                    nn.Linear(in_dim, hidden_dim),
                    nn.LayerNorm(hidden_dim),
                    nn.Mish(inplace=True),
                ]
            )
            in_dim = hidden_dim
        enc_layers.append(nn.Linear(in_dim, latent_dim))
        self.encoder = nn.Sequential(*enc_layers)

        # 解码器：隐变量 z → 重建基因表达向量
        dec_layers = []
        in_dim = latent_dim
        for hidden_dim in reversed(hidden_dims):
            dec_layers.extend(
                [
                    nn.Linear(in_dim, hidden_dim),
                    nn.LayerNorm(hidden_dim),
                    nn.Mish(inplace=True),
                    nn.Dropout(dropout),
                ]
            )
            in_dim = hidden_dim
        dec_layers.append(nn.Linear(in_dim, num_genes))
        self.decoder = nn.Sequential(*dec_layers)

        # 隐空间去噪器：学习逆向扩散过程（DDPM）
        self.denoiser = LatentDenoiser(
            latent_dim=latent_dim,
            hidden_dim=max(128, latent_dim * 4),
            num_layers=3,
            time_dim=128,
            dropout=dropout,
        )

        # 余弦噪声调度（比线性调度更平缓，更适合隐空间扩散）
        betas = self._cosine_beta_schedule(diffusion_steps)
        alphas = 1.0 - betas
        alphas_cumprod = torch.cumprod(alphas, dim=0)
        self.register_buffer("betas", betas)
        self.register_buffer("alphas_cumprod", alphas_cumprod)
        self.register_buffer("sqrt_alphas_cumprod", torch.sqrt(alphas_cumprod))
        self.register_buffer("sqrt_one_minus_alphas_cumprod", torch.sqrt(1.0 - alphas_cumprod))

    def _cosine_beta_schedule(self, timesteps: int, s: float = 0.008) -> torch.Tensor:
        """余弦调度：末期噪声水平更低，更好地匹配隐变量分布."""
        steps = timesteps + 1
        x = torch.linspace(0, timesteps, steps)
        alphas_cumprod = torch.cos(((x / timesteps) + s) / (1 + s) * torch.pi * 0.5) ** 2
        alphas_cumprod = alphas_cumprod / alphas_cumprod[0]
        betas = 1.0 - (alphas_cumprod[1:] / alphas_cumprod[:-1])
        return torch.clamp(betas, 0.0001, 0.02)

    def _extract(self, values: torch.Tensor, timesteps: torch.Tensor, shape: torch.Size) -> torch.Tensor:
        """根据每个样本的时间步收集对应的扩散系数标量."""
        out = values.to(timesteps.device).gather(0, timesteps)
        return out.view(timesteps.shape[0], *((1,) * (len(shape) - 1)))

    def encode(self, x: torch.Tensor) -> torch.Tensor:
        """编码：基因表达 → 隐变量 z（截断防止极端激活）."""
        return torch.clamp(self.encoder(x), min=-10.0, max=10.0)

    def decode(self, z: torch.Tensor) -> torch.Tensor:
        """解码：隐变量 z → 重建基因表达."""
        return self.decoder(z)

    def q_sample(self, z0: torch.Tensor, timesteps: torch.Tensor, noise: torch.Tensor = None) -> torch.Tensor:
        """前向扩散：给定干净隐变量 z0，按调度添加噪声得到 z_t."""
        if noise is None:
            noise = torch.randn_like(z0)
        sqrt_alpha = self._extract(self.sqrt_alphas_cumprod, timesteps, z0.shape)
        sqrt_one_minus = self._extract(self.sqrt_one_minus_alphas_cumprod, timesteps, z0.shape)
        return sqrt_alpha * z0 + sqrt_one_minus * noise

    def masked_mse(self, pred: torch.Tensor, target: torch.Tensor, mask: torch.Tensor = None) -> torch.Tensor:
        """掩码MSE：仅在观测位置（有支持掩码的区域）计算重建误差."""
        loss = F.mse_loss(pred, target, reduction="none")
        if mask is None:
            return loss.mean()
        denom = mask.sum().clamp_min(1.0)
        return (loss * mask).sum() / denom

    def forward(
        self,
        x: torch.Tensor,
        mask: torch.Tensor = None,
        return_recon: bool = True,
        recon_from_denoised: bool = False,
    ) -> dict:
        """前向传播：编码 → 加噪 → 去噪预测 z_0 → 解码重建."""
        z = self.encode(x)
        batch_size = x.shape[0]
        # 随机采样时间步（用于训练去噪器）
        timesteps = torch.randint(0, self.diffusion_steps, (batch_size,), device=x.device)
        z_noisy = self.q_sample(z, timesteps)          # 前向加噪
        z_denoised = self.denoiser(z_noisy, timesteps)  # 预测 z_0
        # 去噪器训练目标：预测值尽量接近原始干净隐变量
        diffusion_loss = F.mse_loss(z_denoised, z)

        x_recon = None
        recon_loss = torch.tensor(0.0, device=x.device)
        if return_recon:
            # 可选：从去噪后的 z_0 而非原始 z 重建（更强的条件生成信号）
            z_for_recon = z_denoised if recon_from_denoised else z
            x_recon = self.decode(z_for_recon)
            # 掩码重建误差：只惩罚观测位置的重构错误
            recon_loss = self.masked_mse(x_recon, x, mask)

        return {
            "z": z,
            "z_denoised": z_denoised,
            "x_recon": x_recon,
            "losses": {
                "diffusion": diffusion_loss,
                "recon": recon_loss,
            },
        }

    @torch.no_grad()
    def denoise_embedding(self, z: torch.Tensor, start_frac: float = 0.35) -> torch.Tensor:
        """确定性隐变量精炼：用DDPM逆向过程从 z 出发做部分去噪，而非从纯噪声采样."""
        self.eval()
        # 从扩散过程的中后段开始（跳过早期高噪声步），减少计算量
        start_t = int(round((self.diffusion_steps - 1) * start_frac))
        start_t = max(0, min(start_t, self.diffusion_steps - 1))
        # 注入零噪声将 z "推送"到 start_t 的含噪分布
        t = torch.full((z.shape[0],), start_t, device=z.device, dtype=torch.long)
        z_t = self.q_sample(z, t, noise=torch.zeros_like(z))

        # DDPM逆向采样：从 start_t 逐步回退到 0
        for step in reversed(range(start_t + 1)):
            step_t = torch.full((z.shape[0],), step, device=z.device, dtype=torch.long)
            pred_x0 = self.denoiser(z_t, step_t)
            if step == 0:
                z_t = pred_x0
                break
            # 从预测的 x_0 反推噪声 ε
            prev_t = torch.full((z.shape[0],), step - 1, device=z.device, dtype=torch.long)
            alpha_t = self._extract(self.alphas_cumprod, step_t, z.shape)
            alpha_prev = self._extract(self.alphas_cumprod, prev_t, z.shape)
            eps = (z_t - torch.sqrt(alpha_t) * pred_x0) / torch.sqrt(1.0 - alpha_t + 1e-8)
            # DDPM逆向步公式：从 z_t 采样 z_{t-1}
            z_t = torch.sqrt(alpha_prev) * pred_x0 + torch.sqrt(1.0 - alpha_prev) * eps
            z_t = torch.clamp(z_t, -10.0, 10.0)
        return torch.clamp(z_t, -10.0, 10.0)

