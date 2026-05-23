import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from sklearn.cluster import KMeans

from ..models import LatentDiffusionAE, SupportMaskNet


class AverageMeter:
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


class ClusterLoss(nn.Module):
    """DEC-style latent geometry regularizer with externally initialized centers."""

    def __init__(self, n_clusters: int, latent_dim: int, alpha: float = 1.0):
        super().__init__()
        self.n_clusters = n_clusters
        self.alpha = alpha
        self.cluster_centers = nn.Parameter(torch.randn(n_clusters, latent_dim) * 0.05)
        self.initialized = False

    def soft_assign(self, z: torch.Tensor) -> torch.Tensor:
        dist = torch.cdist(z, self.cluster_centers).pow(2)
        q = 1.0 / (1.0 + dist / self.alpha)
        q = q.pow((self.alpha + 1.0) / 2.0)
        return q / q.sum(dim=1, keepdim=True).clamp_min(1e-8)

    def forward(self, z: torch.Tensor) -> torch.Tensor:
        q = self.soft_assign(z)
        weight = q.pow(2) / q.sum(dim=0, keepdim=True).clamp_min(1e-8)
        p = weight / weight.sum(dim=1, keepdim=True).clamp_min(1e-8)
        return F.kl_div(torch.log(q.clamp_min(1e-8)), p.detach(), reduction="batchmean")


class ScSpade(nn.Module):
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

        mask_output = self.mask_net(x)
        mask_logits = mask_output["gene_activation_logits"]
        mask_prob = torch.nan_to_num(mask_output["gene_activation_prob"], nan=0.5, posinf=1.0, neginf=0.0)
        mask_prob = mask_prob.clamp(1e-6, 1.0 - 1e-6)

        soft_mask = self._soft_reconstruction_mask(support, mask_prob)
        ae_result = self.diffusion_ae(
            x,
            mask=soft_mask,
            return_recon=return_recon,
            recon_from_denoised=recon_from_denoised,
        )

        z = ae_result["z"]
        losses = {
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
        self.eval()
        z = self.diffusion_ae.encode(x)
        if not use_diffusion:
            return z
        return self.diffusion_ae.denoise_embedding(z, start_frac=diffusion_start_frac)

    @torch.no_grad()
    def predict_mask(self, x: torch.Tensor) -> torch.Tensor:
        self.eval()
        return self.mask_net(x)["gene_activation_prob"].clamp(1e-6, 1.0 - 1e-6)


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
    model.train()
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
        total_loss = (
            weights["mask"] * losses["mask"]
            + weights["diffusion"] * losses["diffusion"]
            + weights["recon"] * losses["recon"]
            + weights["cluster"] * losses["cluster"]
        )

        total_loss.backward()
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


@torch.no_grad()
def extract_embeddings(
    model: nn.Module,
    dataloader,
    device: torch.device,
    diffusion_start_frac: float = 0.35,
    return_masks: bool = False,
) -> dict:
    model.eval()
    direct, diffusion, labels, mask_probs = [], [], [], []
    for x, _, y in dataloader:
        x = x.to(device)
        z_direct = model.get_embedding(x, use_diffusion=False, diffusion_start_frac=diffusion_start_frac)
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


@torch.no_grad()
def initialize_cluster_centers(
    model: ScSpade,
    dataloader,
    device: torch.device,
    n_clusters: int,
    random_state: int = 42,
) -> np.ndarray:
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

