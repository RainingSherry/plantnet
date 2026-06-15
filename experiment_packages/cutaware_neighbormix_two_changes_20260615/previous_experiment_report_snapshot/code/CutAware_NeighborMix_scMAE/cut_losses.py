from __future__ import annotations

import torch
import torch.nn.functional as F


def sinkhorn_balanced_assignment(
    logits: torch.Tensor,
    temperature: float = 0.2,
    iterations: int = 3,
    eps: float = 1e-6,
) -> torch.Tensor:
    """Return approximately balanced soft assignments with Sinkhorn-Knopp."""
    if logits.ndim != 2:
        raise ValueError(f"logits must be [batch, clusters], got {tuple(logits.shape)}")
    q = torch.exp(logits / max(float(temperature), eps)).t()
    q = q / q.sum().clamp_min(eps)
    k, b = q.shape
    for _ in range(max(1, int(iterations))):
        q = q / q.sum(dim=1, keepdim=True).clamp_min(eps)
        q = q / float(k)
        q = q / q.sum(dim=0, keepdim=True).clamp_min(eps)
        q = q / float(b)
    q = (q * float(b)).t()
    return q.clamp_min(eps)


def ot_self_training_loss(
    logits: torch.Tensor,
    temperature: float = 0.2,
    iterations: int = 3,
) -> tuple[torch.Tensor, torch.Tensor, dict]:
    """Balanced self-training target inspired by scCDCG/DCGC OT modules."""
    with torch.no_grad():
        target = sinkhorn_balanced_assignment(logits.detach(), temperature=temperature, iterations=iterations)
        target = target / target.sum(dim=1, keepdim=True).clamp_min(1e-8)
    log_prob = F.log_softmax(logits, dim=1)
    loss = -(target * log_prob).sum(dim=1).mean()
    mass = target.mean(dim=0)
    stats = {
        "ot_assignment_entropy": float((-(target * target.clamp_min(1e-8).log()).sum(dim=1).mean()).detach().cpu()),
        "ot_mass_min": float(mass.min().detach().cpu()),
        "ot_mass_max": float(mass.max().detach().cpu()),
    }
    return loss, target, stats


def graph_cut_loss(
    q_src: torch.Tensor,
    q_dst: torch.Tensor,
    src_rep: torch.Tensor,
    edge_weight: torch.Tensor,
    eps: float = 1e-8,
) -> tuple[torch.Tensor, dict]:
    """
    Mini-batch normalized-cut surrogate.

    This is not a message-passing layer: it never averages embeddings. It only
    asks whether the current soft cluster assignment cuts weak/ambiguous graph
    edges while keeping balanced non-collapsed assignments.
    """
    if q_src.ndim != 2 or q_dst.ndim != 2:
        raise ValueError("q_src and q_dst must be [n, clusters] tensors.")
    edge_weight = edge_weight.to(dtype=q_src.dtype, device=q_src.device).view(-1)
    q_anchor = q_src[src_rep]
    same_prob = (q_anchor * q_dst).sum(dim=1).clamp(0.0, 1.0)
    cut_per_edge = edge_weight * (1.0 - same_prob)
    cut_num = cut_per_edge.mean()

    degree = torch.zeros(q_src.shape[0], dtype=q_src.dtype, device=q_src.device)
    degree.scatter_add_(0, src_rep, edge_weight)
    assoc = (degree[:, None] * q_src).sum(dim=0).clamp_min(eps)
    cut_by_cluster = (edge_weight[:, None] * q_anchor * (1.0 - q_dst)).sum(dim=0)
    ncut = (cut_by_cluster / assoc).mean()

    mean_q = q_src.mean(dim=0)
    uniform = torch.full_like(mean_q, 1.0 / float(q_src.shape[1]))
    balance = F.mse_loss(mean_q, uniform)
    entropy = -(q_src * q_src.clamp_min(eps).log()).sum(dim=1).mean()
    stats = {
        "cut_loss": float(ncut.detach().cpu()),
        "cut_num": float(cut_num.detach().cpu()),
        "cut_balance_mse": float(balance.detach().cpu()),
        "assignment_entropy": float(entropy.detach().cpu()),
        "cluster_mass_min": float(mean_q.min().detach().cpu()),
        "cluster_mass_max": float(mean_q.max().detach().cpu()),
    }
    return ncut + balance, stats


def attention_fusion_probe(
    q_expr: torch.Tensor,
    q_graph: torch.Tensor,
    reliability: torch.Tensor,
) -> tuple[torch.Tensor, dict]:
    """
    AttentionAE-sc style diagnostic fusion without making attention the main mechanism.

    The graph assignment is only trusted in proportion to reliability. This is
    intentionally a probe because attention alone cannot cut bad edges.
    """
    reliability = reliability.to(dtype=q_expr.dtype, device=q_expr.device).view(-1, 1).clamp(0.0, 1.0)
    fused = (1.0 - reliability) * q_expr + reliability * q_graph
    fused = fused / fused.sum(dim=1, keepdim=True).clamp_min(1e-8)
    disagreement = torch.mean(torch.sum(torch.abs(q_expr - q_graph), dim=1))
    stats = {
        "attention_probe_reliability": float(reliability.mean().detach().cpu()),
        "attention_probe_disagreement": float(disagreement.detach().cpu()),
    }
    return fused, stats
