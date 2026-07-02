from __future__ import annotations

import torch
from torch import nn
import torch.nn.functional as F


class DINOv2KoLeoLoss(nn.Module):
    """scMAE reconstruction plus DINOv2 class/patch distillation and KoLeo spreading."""

    def __init__(
        self,
        cls_out_dim: int = 128,
        patch_out_dim: int = 128,
        masked_data_weight: float = 0.75,
        mask_weight: float = 0.60,
        cls_distill_weight: float = 0.08,
        patch_distill_weight: float = 0.12,
        koleo_weight: float = 0.02,
        variance_weight: float = 0.01,
        student_temp: float = 0.10,
        teacher_temp: float = 0.07,
        teacher_patch_temp: float = 0.07,
        center_momentum: float = 0.90,
        patch_center_momentum: float = 0.90,
    ):
        super().__init__()
        self.masked_data_weight = float(masked_data_weight)
        self.mask_weight = float(mask_weight)
        self.cls_distill_weight = float(cls_distill_weight)
        self.patch_distill_weight = float(patch_distill_weight)
        self.koleo_weight = float(koleo_weight)
        self.variance_weight = float(variance_weight)
        self.student_temp = float(student_temp)
        self.teacher_temp = float(teacher_temp)
        self.teacher_patch_temp = float(teacher_patch_temp)
        self.center_momentum = float(center_momentum)
        self.patch_center_momentum = float(patch_center_momentum)
        self.register_buffer("center", torch.zeros(1, int(cls_out_dim)))
        self.register_buffer("patch_center", torch.zeros(1, 1, int(patch_out_dim)))

    @staticmethod
    def variance_loss(z: torch.Tensor) -> torch.Tensor:
        std = torch.sqrt(z.var(dim=0) + 1e-4)
        return F.relu(0.5 - std).mean()

    @staticmethod
    def koleo_loss(z: torch.Tensor) -> torch.Tensor:
        z = F.normalize(z, dim=-1)
        dots = z @ z.t()
        n = dots.shape[0]
        dots.view(-1)[:: n + 1].fill_(-1.0)
        nn_idx = torch.argmax(dots, dim=1)
        dist = torch.norm(z - z[nn_idx], dim=1)
        return -torch.log(dist + 1e-8).mean()

    @staticmethod
    def soft_ce(student_logits: torch.Tensor, teacher_prob: torch.Tensor, temp: float) -> torch.Tensor:
        logp = F.log_softmax(student_logits / float(temp), dim=-1)
        return torch.sum(-teacher_prob * logp, dim=-1)

    @torch.no_grad()
    def update_center(self, teacher_cls: torch.Tensor, teacher_patch: torch.Tensor) -> None:
        cls_center = teacher_cls.mean(dim=0, keepdim=True)
        patch_center = teacher_patch.mean(dim=(0, 1), keepdim=True)
        self.center.mul_(self.center_momentum).add_(cls_center, alpha=1.0 - self.center_momentum)
        self.patch_center.mul_(self.patch_center_momentum).add_(patch_center, alpha=1.0 - self.patch_center_momentum)

    def forward(
        self,
        student: dict[str, torch.Tensor],
        teacher: dict[str, torch.Tensor],
        target_expr: torch.Tensor,
        gene_mask: torch.Tensor,
        module_mask: torch.Tensor,
        distill_weight_scale: float = 1.0,
    ) -> tuple[torch.Tensor, dict[str, float]]:
        weights = gene_mask * self.masked_data_weight + (1.0 - gene_mask) * (1.0 - self.masked_data_weight)
        rec = (weights * F.smooth_l1_loss(student["reconstruction"], target_expr, reduction="none")).mean()
        mask_loss = F.binary_cross_entropy_with_logits(student["mask_logits"], gene_mask.float())
        scmae = (1.0 - self.mask_weight) * rec + self.mask_weight * mask_loss

        teacher_cls = teacher["class_logits"].detach()
        teacher_patch = teacher["patch_logits"].detach()
        cls_prob = F.softmax((teacher_cls - self.center) / self.teacher_temp, dim=-1).detach()
        patch_prob = F.softmax((teacher_patch - self.patch_center) / self.teacher_patch_temp, dim=-1).detach()
        cls_loss = self.soft_ce(student["class_logits"], cls_prob, self.student_temp).mean()
        patch_ce = self.soft_ce(student["patch_logits"], patch_prob, self.student_temp)
        patch_weight = module_mask.float()
        patch_loss = (patch_ce * patch_weight).sum() / patch_weight.sum().clamp_min(1.0)
        var_loss = self.variance_loss(student["latent"])
        koleo = self.koleo_loss(student["latent"])
        distill_scale = float(distill_weight_scale)
        total = (
            scmae
            + distill_scale * self.cls_distill_weight * cls_loss
            + distill_scale * self.patch_distill_weight * patch_loss
            + self.koleo_weight * koleo
            + self.variance_weight * var_loss
        )
        self.update_center(teacher_cls, teacher_patch)
        with torch.no_grad():
            cls_entropy = -(cls_prob * cls_prob.clamp_min(1e-8).log()).sum(dim=1).mean()
            patch_entropy = -(patch_prob * patch_prob.clamp_min(1e-8).log()).sum(dim=-1).mean()
            cls_conf = cls_prob.max(dim=1).values.mean()
            patch_conf = patch_prob.max(dim=-1).values.mean()
        return total, {
            "loss": float(total.detach().cpu()),
            "scmae_loss": float(scmae.detach().cpu()),
            "reconstruction_loss": float(rec.detach().cpu()),
            "mask_loss": float(mask_loss.detach().cpu()),
            "cls_distill_loss": float(cls_loss.detach().cpu()),
            "patch_distill_loss": float(patch_loss.detach().cpu()),
            "koleo_loss": float(koleo.detach().cpu()),
            "variance_loss": float(var_loss.detach().cpu()),
            "teacher_cls_entropy": float(cls_entropy.detach().cpu()),
            "teacher_patch_entropy": float(patch_entropy.detach().cpu()),
            "teacher_cls_confidence": float(cls_conf.detach().cpu()),
            "teacher_patch_confidence": float(patch_conf.detach().cpu()),
            "distill_weight_scale": distill_scale,
        }
