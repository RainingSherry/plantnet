import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from sklearn.cluster import KMeans

from ..models import LatentDiffusionAE, SupportMaskNet


# ---------- 指标记录器 ----------
class AverageMeter:
    """滚动平均指标，用于跟踪训练中各损失项的均值."""
    def __init__(self):
        self.reset()

    def reset(self):
        self.sum = 0.0
        self.count = 0
        self.avg = 0.0

    def update(self, value: float, n: int):
        self.sum += value * n
        self.count += n
        self.avg = self.sum / max(self.count, 1)


# ---------- DEC风格聚类损失 ----------
# 通过可学习的聚类中心，将隐空间几何正则化为K个紧凑簇
class ClusterLoss(nn.Module):
    """DEC风格隐空间几何正则化：使细胞嵌入聚集在可学习的聚类中心周围.

    工作原理：用t-分布核计算每个细胞属于各聚类的软分配概率 Q，
    再构造目标分布 P（Q的平方归一化），通过 KL(P‖Q) 让 Q 逼近 P，
    从而使嵌入向聚类中心聚拢（类似 Deep Embedded Clustering）。
    """

    def __init__(self, n_clusters: int, latent_dim: int, alpha: float = 1.0):
        super().__init__()
        self.n_clusters = n_clusters
        self.alpha = alpha
        self.cluster_centers = nn.Parameter(torch.randn(n_clusters, latent_dim) * 0.05)
        self.initialized = False

    def soft_assign(self, z: torch.Tensor) -> torch.Tensor:
        """t-分布核：计算细胞嵌入与每个聚类中心的软分配概率 Q."""
        dist = torch.cdist(z, self.cluster_centers).pow(2)
        q = 1.0 / (1.0 + dist / self.alpha)
        q = q.pow((self.alpha + 1.0) / 2.0)
        return q / q.sum(dim=1, keepdim=True).clamp_min(1e-8)

    def forward(self, z: torch.Tensor) -> torch.Tensor:
        """KL(P‖Q)：Q为软分配，P为辅助目标分布，使隐空间向聚类中心收紧."""
        q = self.soft_assign(z)
        weight = q.pow(2) / q.sum(dim=0, keepdim=True).clamp_min(1e-8)
        p = weight / weight.sum(dim=1, keepdim=True).clamp_min(1e-8)
        return F.kl_div(torch.log(q.clamp_min(1e-8)), p.detach(), reduction="batchmean")


# ---------- 完整ScSpade模型 ----------
# 三大核心组件联合训练：支持掩码预测 + 扩散自编码器 + 可选聚类正则化
class ScSpade(nn.Module):
    """支持掩码驱动的隐扩散自编码器，用于单细胞表达数据.

    三大任务联合学习：
    1. Support Mask Net：预测哪些零值是"可能表达但被dropout"（而非生物学真零）
    2. Latent Diffusion AE：在隐空间做去噪扩散，增强嵌入的鲁棒性
    3. Cluster Loss（可选）：将隐空间几何正则化为K个紧凑簇
    """

    def __init__(
        self,
        num_genes: int,
        n_clusters: int,
        latent_dim: int = 32,
        mask_hidden_dims=None,
        diffusion_hidden_dims=None,
        diffusion_steps: int = 100,
        dropout: float = 0.1,
        mask_coupling: str = "weighted_observed",
    ):
        super().__init__()
        self.mask_coupling = mask_coupling
        self.mask_net = SupportMaskNet(num_genes=num_genes, hidden_dims=mask_hidden_dims, dropout=dropout)
        self.diffusion_ae = LatentDiffusionAE(
            num_genes=num_genes,
            latent_dim=latent_dim,
            hidden_dims=diffusion_hidden_dims,
            diffusion_steps=diffusion_steps,
            dropout=dropout,
        )
        self.cluster_loss = ClusterLoss(n_clusters=n_clusters, latent_dim=latent_dim)

    def _soft_reconstruction_mask(self, support: torch.Tensor, mask_prob: torch.Tensor) -> torch.Tensor:
        """将观测支持掩码与预测的基因激活概率融合，作为重建加权掩码.

        三种策略：
        - prob：观测支持 AND 预测激活概率（最严格）
        - observed：仅用观测支持（最宽松）
        - weighted_observed（默认）：两者加权平均（平衡策略）
        """
        if self.mask_coupling == "prob":
            return torch.clamp(support * mask_prob.detach(), 0.0, 1.0)
        if self.mask_coupling == "observed":
            return support
        return support * (0.5 + 0.5 * mask_prob.detach())

    def forward(
        self,
        x: torch.Tensor,
        support: torch.Tensor = None,
        return_recon: bool = True,
        recon_from_denoised: bool = False,
    ) -> dict:
        if support is None:
            support = (x > 0).float()

        # 步骤1：预测基因激活概率（区分真零与dropout零）
        mask_output = self.mask_net(x)
        mask_logits = mask_output["gene_activation_logits"]
        mask_prob = torch.nan_to_num(mask_output["gene_activation_prob"], nan=0.5, posinf=1.0, neginf=0.0)
        mask_prob = mask_prob.clamp(1e-6, 1.0 - 1e-6)

        # 步骤2：融合观测支持与预测支持，生成用于重建加权的软掩码
        soft_mask = self._soft_reconstruction_mask(support, mask_prob)

        # 步骤3：编码 → 加噪 → 去噪 → 解码（全量扩散AE前向）
        ae_result = self.diffusion_ae(
            x,
            mask=soft_mask,
            return_recon=return_recon,
            recon_from_denoised=recon_from_denoised,
        )

        z = ae_result["z"]
        losses = {
            # mask_loss：预测的支持概率 vs 实际观测支持（BCE损失）
            # diffusion_loss：去噪后隐变量 vs 原始干净隐变量（MSE）
            # recon_loss：掩码重建误差（仅在观测位置计算）
            # cluster_loss：DEC聚类正则化（初始化后才生效，默认关闭）
            "mask": F.binary_cross_entropy_with_logits(mask_logits, support, reduction="mean"),
            "diffusion": ae_result["losses"]["diffusion"],
            "recon": ae_result["losses"]["recon"],
            "cluster": torch.tensor(0.0, device=x.device),
        }
        if self.training and self.cluster_loss.initialized:
            losses["cluster"] = self.cluster_loss(z)

        return {
            "mask_prob": mask_prob,
            "soft_mask": soft_mask,
            "z": z,
            "z_denoised": ae_result["z_denoised"],
            "x_recon": ae_result["x_recon"],
            "losses": losses,
        }

    @torch.no_grad()
    def get_embedding(self, x: torch.Tensor, use_diffusion: bool = False, diffusion_start_frac: float = 0.35):
        """提取细胞嵌入：use_diffusion=True 时额外做DDPM逆向去噪以精炼嵌入."""
        self.eval()
        z = self.diffusion_ae.encode(x)
        if not use_diffusion:
            return z
        return self.diffusion_ae.denoise_embedding(z, start_frac=diffusion_start_frac)

    @torch.no_grad()
    def predict_mask(self, x: torch.Tensor) -> torch.Tensor:
        """推理时预测每个(cell,gene)的基因激活概率."""
        self.eval()
        return self.mask_net(x)["gene_activation_prob"].clamp(1e-6, 1.0 - 1e-6)


# ---------- 训练计划函数 ----------
# 三阶段课程学习：纯AE热启动 → 扩散损失渐进引入 → 全量联合训练（含可选聚类损失）


def _scheduled_weights(
    epoch: int,
    warmup_epochs: int,
    diffusion_ramp_epochs: int,
    cluster_warmup_epochs: int,
    mask_weight: float,
    diffusion_weight: float,
    recon_weight: float,
    cluster_weight: float,
    diffusion_warmup_weight: float,
) -> dict:
    """根据当前epoch返回各损失项的权重及训练阶段名称.

    阶段1（AE_WARMUP）：仅训练编解码器重建，去噪器未训练，不加扩散损失
    阶段2（DIFF_RAMP）：逐步引入扩散损失，去噪器学习逆向过程
    阶段3（DIFF_FULL / JOINT）：全量联合训练，聚类损失在cluster_warmup后激活
    """
    if epoch < warmup_epochs:
        return {
            "mask": mask_weight,
            "diffusion": 0.0,
            "recon": recon_weight,
            "cluster": 0.0,
            "phase": "AE_WARMUP",
            "recon_from_denoised": False,
        }
    if epoch < warmup_epochs + diffusion_ramp_epochs:
        return {
            "mask": mask_weight,
            "diffusion": min(diffusion_weight, diffusion_warmup_weight),
            "recon": recon_weight,
            "cluster": 0.0,
            "phase": "DIFF_RAMP",
            "recon_from_denoised": False,
        }
    cluster_start = warmup_epochs + diffusion_ramp_epochs + cluster_warmup_epochs
    return {
        "mask": mask_weight,
        "diffusion": diffusion_weight,
        "recon": recon_weight,
        "cluster": cluster_weight if epoch >= cluster_start else 0.0,
        "phase": "JOINT" if epoch >= cluster_start and cluster_weight > 0 else "DIFF_FULL",
        "recon_from_denoised": False,
    }


# ---------- 单轮训练 ----------
def train_epoch(
    model: nn.Module,
    dataloader,
    optimizer,
    device: torch.device,
    epoch: int,
    warmup_epochs: int = 30,
    diffusion_ramp_epochs: int = 50,
    cluster_warmup_epochs: int = 0,
    mask_weight: float = 0.2,
    diffusion_weight: float = 0.05,
    recon_weight: float = 0.8,
    cluster_weight: float = 0.0,
    diffusion_warmup_weight: float = 0.05,
    grad_clip: float = 1.0,
) -> dict:
    """执行一轮前向-反向-优化步骤，返回各损失项的滑动平均值."""
    model.train()
    # 查询当前epoch对应的训练阶段和损失权重
    weights = _scheduled_weights(
        epoch,
        warmup_epochs,
        diffusion_ramp_epochs,
        cluster_warmup_epochs,
        mask_weight,
        diffusion_weight,
        recon_weight,
        cluster_weight,
        diffusion_warmup_weight,
    )

    meters = {name: AverageMeter() for name in ["loss", "mask", "diffusion", "recon", "cluster"]}
    for x, support, _ in dataloader:
        x = x.to(device)
        support = support.to(device)

        optimizer.zero_grad(set_to_none=True)
        result = model(
            x,
            support=support,
            return_recon=True,
            recon_from_denoised=weights["recon_from_denoised"],
        )
        losses = result["losses"]
        # 加权组合四个损失项
        total_loss = (
            weights["mask"] * losses["mask"]
            + weights["diffusion"] * losses["diffusion"]
            + weights["recon"] * losses["recon"]
            + weights["cluster"] * losses["cluster"]
        )

        total_loss.backward()
        # 梯度裁剪防止梯度爆炸
        if grad_clip and grad_clip > 0:
            torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=grad_clip)
        optimizer.step()

        n = x.shape[0]
        meters["loss"].update(float(total_loss.detach().cpu()), n)
        for name in ["mask", "diffusion", "recon", "cluster"]:
            meters[name].update(float(losses[name].detach().cpu()), n)

    out = {name if name == "loss" else f"{name}_loss": meter.avg for name, meter in meters.items()}
    out.update({f"weight_{key}": value for key, value in weights.items() if key in {"mask", "diffusion", "recon", "cluster"}})
    out["phase"] = weights["phase"]
    return out


# ---------- 嵌入提取 ----------
# 支持直接提取（编码器输出）和扩散精炼提取（额外经过DDPM逆向过程）
@torch.no_grad()
def extract_embeddings(
    model: nn.Module,
    dataloader,
    device: torch.device,
    diffusion_start_frac: float = 0.35,
    return_masks: bool = False,
) -> dict:
    """提取全体细胞的直接嵌入和扩散精炼嵌入，用于下游聚类评估."""
    model.eval()
    direct, diffusion, labels, mask_probs = [], [], [], []
    for x, _, y in dataloader:
        x = x.to(device)
        # 直接嵌入：编码器输出 z（不经过去噪）
        z_direct = model.get_embedding(x, use_diffusion=False, diffusion_start_frac=diffusion_start_frac)
        # 扩散嵌入：额外经过DDPM逆向过程精炼后的 z
        z_diff = model.get_embedding(x, use_diffusion=True, diffusion_start_frac=diffusion_start_frac)
        direct.append(torch.nan_to_num(z_direct, nan=0.0, posinf=0.0, neginf=0.0).cpu())
        diffusion.append(torch.nan_to_num(z_diff, nan=0.0, posinf=0.0, neginf=0.0).cpu())
        labels.append(y.cpu())
        if return_masks:
            mask_probs.append(model.predict_mask(x).cpu())

    result = {
        "direct": torch.cat(direct, dim=0).numpy(),
        "diffusion": torch.cat(diffusion, dim=0).numpy(),
        "labels": torch.cat(labels, dim=0).numpy(),
    }
    if return_masks:
        result["mask_probs"] = torch.cat(mask_probs, dim=0).numpy()
    return result


# ---------- 聚类中心初始化 ----------
# 用KMeans在直接嵌入空间中初始化DEC的可学习聚类中心
@torch.no_grad()
def initialize_cluster_centers(
    model: ScSpade,
    dataloader,
    device: torch.device,
    n_clusters: int,
    random_state: int = 42,
) -> np.ndarray:
    """先用KMeans在直接嵌入上做聚类，用聚类中心初始化DEC的cluster_centers参数."""
    model.eval()
    embeddings = []
    for x, _, _ in dataloader:
        x = x.to(device)
        embeddings.append(model.get_embedding(x, use_diffusion=False).cpu())
    embeddings = torch.cat(embeddings, dim=0).numpy()
    kmeans = KMeans(n_clusters=n_clusters, random_state=random_state, n_init=20)
    pred = kmeans.fit_predict(embeddings)
    centers = torch.tensor(kmeans.cluster_centers_, dtype=torch.float32, device=device)
    model.cluster_loss.cluster_centers.data.copy_(centers)
    model.cluster_loss.initialized = True
    return pred

