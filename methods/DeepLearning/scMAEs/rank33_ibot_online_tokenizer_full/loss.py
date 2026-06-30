from __future__ import annotations

from dataclasses import dataclass

import torch
import torch.nn.functional as F


@dataclass
class IBOTLossParts:
    total: torch.Tensor
    cls: torch.Tensor
    patch: torch.Tensor
    reconstruction: torch.Tensor
    mask: torch.Tensor
    mask_rate: torch.Tensor


def make_patch_mask(batch_size: int, num_patches: int, mask_ratio: float, device: torch.device) -> torch.Tensor:
    if batch_size <= 0 or num_patches <= 0:
        raise ValueError("batch_size and num_patches must be positive")
    if not 0.0 < float(mask_ratio) < 1.0:
        raise ValueError("mask_ratio must be in (0, 1)")
    mask = (torch.rand(batch_size, num_patches, device=device) < float(mask_ratio)).float()
    empty = mask.sum(dim=1) == 0
    if empty.any():
        chosen = torch.randint(0, num_patches, (int(empty.sum()),), device=device)
        mask[empty, chosen] = 1.0
    return mask


def patch_mask_to_gene_mask(patch_mask: torch.Tensor, patch_size: int, num_genes: int) -> torch.Tensor:
    if patch_mask.ndim != 2:
        raise ValueError(f"patch_mask must be [batch, patches], got {tuple(patch_mask.shape)}")
    return patch_mask.repeat_interleave(int(patch_size), dim=1)[:, : int(num_genes)]


def soft_cross_entropy(student_logits: torch.Tensor, teacher_probs: torch.Tensor) -> torch.Tensor:
    if student_logits.shape != teacher_probs.shape:
        raise ValueError("student_logits and teacher_probs must share shape")
    return -(teacher_probs.detach() * F.log_softmax(student_logits, dim=-1)).sum(dim=-1)


def masked_reconstruction_loss(reconstruction: torch.Tensor, clean: torch.Tensor, gene_mask: torch.Tensor) -> torch.Tensor:
    if reconstruction.shape != clean.shape or gene_mask.shape != clean.shape:
        raise ValueError("reconstruction, clean, and gene_mask must share [batch, genes] shape")
    denom = gene_mask.sum().clamp_min(1.0)
    return (F.smooth_l1_loss(reconstruction, clean, reduction="none") * gene_mask).sum() / denom


def ibot_distillation_loss(
    outputs: dict[str, dict[str, torch.Tensor]],
    clean: torch.Tensor,
    mask1: torch.Tensor,
    mask2: torch.Tensor,
    gene_mask1: torch.Tensor,
    gene_mask2: torch.Tensor,
    cls_center: torch.Tensor,
    patch_center: torch.Tensor,
    *,
    student_temp: float,
    teacher_temp: float,
    teacher_patch_temp: float,
    cls_weight: float,
    patch_weight: float,
    reconstruction_weight: float,
    mask_weight: float,
) -> IBOTLossParts:
    s1 = outputs["student1"]
    s2 = outputs["student2"]
    t1 = outputs["teacher1"]
    t2 = outputs["teacher2"]
    for name, student in (("student1", s1), ("student2", s2)):
        if student["cls_logits"].ndim != 2 or student["patch_logits"].ndim != 3:
            raise ValueError(f"{name} logits must be [batch, out_dim] and [batch, patches, out_dim]")
    if s1["patch_logits"].shape[:2] != mask1.shape or s2["patch_logits"].shape[:2] != mask2.shape:
        raise ValueError("patch masks must match student patch logits [batch, patches]")
    if t1["patch_logits"].shape != s1["patch_logits"].shape or t2["patch_logits"].shape != s2["patch_logits"].shape:
        raise ValueError("teacher and student patch logits must share [batch, patches, out_dim] shape")
    if cls_center.shape != (1, s1["cls_logits"].shape[-1]):
        raise ValueError("cls_center must be [1, out_dim]")
    if patch_center.shape != (1, 1, s1["patch_logits"].shape[-1]):
        raise ValueError("patch_center must be [1, 1, out_dim]")

    tp_cls1 = F.softmax((t1["cls_logits"] - cls_center) / float(teacher_temp), dim=-1).detach()
    tp_cls2 = F.softmax((t2["cls_logits"] - cls_center) / float(teacher_temp), dim=-1).detach()
    tp_patch1 = F.softmax((t1["patch_logits"] - patch_center) / float(teacher_patch_temp), dim=-1).detach()
    tp_patch2 = F.softmax((t2["patch_logits"] - patch_center) / float(teacher_patch_temp), dim=-1).detach()

    cls_loss = 0.5 * (
        soft_cross_entropy(s1["cls_logits"] / float(student_temp), tp_cls2).mean()
        + soft_cross_entropy(s2["cls_logits"] / float(student_temp), tp_cls1).mean()
    )

    patch_loss1 = soft_cross_entropy(s1["patch_logits"] / float(student_temp), tp_patch1) * mask1.to(dtype=s1["patch_logits"].dtype)
    patch_loss2 = soft_cross_entropy(s2["patch_logits"] / float(student_temp), tp_patch2) * mask2.to(dtype=s2["patch_logits"].dtype)
    patch_denom = (mask1.sum() + mask2.sum()).clamp_min(1.0)
    patch_loss = (patch_loss1.sum() + patch_loss2.sum()) / patch_denom

    reconstruction = 0.5 * (
        masked_reconstruction_loss(s1["reconstruction"], clean, gene_mask1)
        + masked_reconstruction_loss(s2["reconstruction"], clean, gene_mask2)
    )
    mask_loss = 0.5 * (
        F.binary_cross_entropy_with_logits(s1["mask_logits"], gene_mask1)
        + F.binary_cross_entropy_with_logits(s2["mask_logits"], gene_mask2)
    )
    total = (
        float(cls_weight) * cls_loss
        + float(patch_weight) * patch_loss
        + float(reconstruction_weight) * reconstruction
        + float(mask_weight) * mask_loss
    )
    return IBOTLossParts(
        total=total,
        cls=cls_loss,
        patch=patch_loss,
        reconstruction=reconstruction,
        mask=mask_loss,
        mask_rate=0.5 * (gene_mask1.mean() + gene_mask2.mean()),
    )
