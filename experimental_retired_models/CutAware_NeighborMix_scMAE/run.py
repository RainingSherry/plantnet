#!/usr/bin/env python
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import torch
import torch.nn.functional as F
from torch.utils.data import DataLoader

CURRENT_DIR = Path(__file__).resolve().parent
ROOT = next(parent for parent in [CURRENT_DIR, *CURRENT_DIR.parents] if (parent / "methods" / "DeepLearning" / "scMAE_family.py").exists())
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from methods.DeepLearning import scMAE_family as family
from experimental_retired_models.CutAware_NeighborMix_scMAE.cut_losses import (
    attention_fusion_probe,
    graph_cut_loss,
    ot_self_training_loss,
)
from experimental_retired_models.CutAware_NeighborMix_scMAE.diagnostics import (
    embedding_geometry,
    mapped_predictions,
    pairwise_similarity_digest,
    per_cell_type_metrics,
)
from experimental_retired_models.CutAware_NeighborMix_scMAE.mixing import make_gated_neighbor_mixed_batch, make_neighbor_mixed_batch
from experimental_retired_models.CutAware_NeighborMix_scMAE.model import CutAwareAutoEncoder
from experimental_retired_models.CutAware_NeighborMix_scMAE.neighbor_graph import (
    apply_cluster_cut_reweight,
    build_embedding_knn_graph,
    build_pca_knn_graph,
    compute_edge_weights,
    graph_cut_diagnostics,
    neighbor_tensors_for_batch,
)
from methods.shared_utils import ensure_dir, save_json, sanitize_anndata_for_write


def str2bool(value):
    return family.str2bool(value)


def parse_args():
    parser = argparse.ArgumentParser(description="Cut-aware NeighborMix-scMAE experimental runner")
    parser.add_argument("--data_path", required=True)
    parser.add_argument("--save_dir", required=True)
    parser.add_argument("--dataset_name", default=None)
    parser.add_argument("--method_name", default="CutAware_NeighborMix_scMAE")
    parser.add_argument(
        "--variant_name",
        default="canm_cut_ot",
        choices=[
            "canm_diagnostic_only",
            "canm_cut_ot",
            "canm_mix_plus_cut",
            "canm_cut_ot_warm",
            "canm_mix_plus_cut_warm",
            "canm_cut_reweighted_mix",
            "canm_cut_reweighted_mix_contrast",
            "canm_gated_cut_mix",
            "canm_gated_cut_warm",
            "canm_attention_fusion_probe",
        ],
    )
    parser.add_argument("--label_key", default="auto")
    parser.add_argument("--input_mode", default="auto", choices=["auto", "raw", "log1p"])
    parser.add_argument("--n_top_genes", type=int, default=1000)
    parser.add_argument("--target_sum", type=float, default=10000.0)
    parser.add_argument("--scale_input", type=str2bool, default=True)
    parser.add_argument("--n_clusters", type=int, default=0)
    parser.add_argument("--hidden_size", type=int, default=128)
    parser.add_argument("--hidden_dim", type=int, default=None)
    parser.add_argument("--latent_dim", type=int, default=None)
    parser.add_argument("--dropout", type=float, default=0.0)
    parser.add_argument("--masked_data_weight", type=float, default=0.75)
    parser.add_argument("--mask_loss_weight", type=float, default=0.7)
    parser.add_argument("--epochs", type=int, default=80)
    parser.add_argument("--batch_size", type=int, default=256)
    parser.add_argument("--lr", type=float, default=1e-3)
    parser.add_argument("--mask_ratio", type=float, default=0.4)
    parser.add_argument("--neighbor_k", type=int, default=10)
    parser.add_argument("--cut_neighbors", type=int, default=10)
    parser.add_argument("--mix_neighbors", type=int, default=4)
    parser.add_argument("--mix_alpha", type=float, default=0.8)
    parser.add_argument("--tau", type=float, default=0.2)
    parser.add_argument("--knn_pca_dim", type=int, default=50)
    parser.add_argument("--edge_reliability_mode", default="sim_mutual_snn_distance")
    parser.add_argument("--gamma_sim", type=float, default=1.0)
    parser.add_argument("--gamma_mutual", type=float, default=1.0)
    parser.add_argument("--gamma_snn", type=float, default=1.0)
    parser.add_argument("--gamma_distance", type=float, default=1.0)
    parser.add_argument("--edge_prune_quantile", type=float, default=0.0)
    parser.add_argument("--cut_cross_weight", type=float, default=0.05)
    parser.add_argument("--pseudo_weight", type=float, default=0.2)
    parser.add_argument("--cut_weight", type=float, default=0.2)
    parser.add_argument("--ot_weight", type=float, default=0.1)
    parser.add_argument("--gate_prior_weight", type=float, default=0.02)
    parser.add_argument("--gate_entropy_weight", type=float, default=0.001)
    parser.add_argument("--gate_cluster_weight", type=float, default=0.02)
    parser.add_argument("--gate_start_epoch", type=int, default=1)
    parser.add_argument("--gate_cluster_start_epoch", type=int, default=10)
    parser.add_argument("--gate_temperature", type=float, default=1.0)
    parser.add_argument("--gate_min", type=float, default=0.05)
    parser.add_argument("--attention_probe_weight", type=float, default=0.0)
    parser.add_argument("--contrast_weight", type=float, default=0.05)
    parser.add_argument("--contrast_temperature", type=float, default=0.2)
    parser.add_argument("--contrast_start_epoch", type=int, default=5)
    parser.add_argument("--contrast_neighbor_positive_weight", type=float, default=0.5)
    parser.add_argument("--contrast_hard_negative_weight", type=float, default=1.0)
    parser.add_argument("--contrast_projection_dim", type=int, default=128)
    parser.add_argument("--ot_temperature", type=float, default=0.2)
    parser.add_argument("--ot_iterations", type=int, default=3)
    parser.add_argument("--cluster_temperature", type=float, default=1.0)
    parser.add_argument("--cut_start_epoch", type=int, default=1)
    parser.add_argument("--ot_start_epoch", type=int, default=1)
    parser.add_argument("--graph_refresh_interval", type=int, default=0)
    parser.add_argument("--graph_refresh_start_epoch", type=int, default=20)
    parser.add_argument("--grad_clip", type=float, default=5.0)
    parser.add_argument("--skip_eval", type=str2bool, default=False)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--gpu", type=int, default=1)
    parser.add_argument("--no_cuda", action="store_true")
    parser.add_argument("--no_save_h5ad", action="store_true")
    return parser.parse_args()


def make_loaders(data_np: np.ndarray, labels: np.ndarray, batch_size: int, seed: int):
    dataset = family.IndexedExpressionDataset(data_np, labels)
    generator = torch.Generator()
    generator.manual_seed(seed)
    train_loader = DataLoader(dataset, batch_size=batch_size, shuffle=True, drop_last=False, generator=generator)
    eval_loader = DataLoader(dataset, batch_size=max(batch_size * 4, 512), shuffle=False, drop_last=False)
    return train_loader, eval_loader


@torch.no_grad()
def extract_cluster_probs(model: CutAwareAutoEncoder, loader: DataLoader, device: torch.device, temperature: float) -> tuple[np.ndarray, np.ndarray]:
    model.eval()
    probs = []
    labels = []
    for _, x, y in loader:
        latent = model.encoder(x.to(device))
        probs.append(model.cluster_probs(latent, temperature=temperature).detach().cpu().numpy())
        labels.append(y.numpy())
    return np.concatenate(probs, axis=0).astype(np.float32), np.concatenate(labels, axis=0).astype(np.int64)


def save_graph_outputs(save_dir: Path, graph, edge_reliability, edge_weights, edge_summary):
    np.save(save_dir / "neighbor_indices.npy", graph.indices)
    np.save(save_dir / "neighbor_base_probs.npy", graph.probs)
    np.save(save_dir / "neighbor_similarity.npy", graph.similarity)
    np.save(save_dir / "neighbor_distance.npy", graph.distance)
    np.save(save_dir / "edge_reliability.npy", edge_reliability)
    save_json(graph.profile, str(save_dir / "neighbor_graph_profile.json"))
    save_json(edge_summary, str(save_dir / "edge_weight_summary.json"))
    save_json({"edge_reliability_mode": edge_summary.get("edge_reliability_mode", ""), **edge_summary}, str(save_dir / "neighbor_reliability_summary.json"))


def maybe_apply_cut_reweight(args, save_dir: Path, graph, edge_weights, edge_summary, n_clusters: int, suffix: str = ""):
    cut_reweight_summary = {}
    cut_labels = None
    if args.variant_name in {"canm_cut_reweighted_mix", "canm_cut_reweighted_mix_contrast", "canm_gated_cut_mix", "canm_gated_cut_warm"}:
        cut_labels, edge_weights, cut_reweight_summary = apply_cluster_cut_reweight(
            graph=graph,
            edge_weights=edge_weights,
            n_clusters=n_clusters,
            cross_weight=args.cut_cross_weight,
            seed=args.seed,
        )
        np.save(save_dir / f"graph_cut_seed_clusters{suffix}.npy", cut_labels)
        save_json(cut_reweight_summary, str(save_dir / f"cut_reweight_summary{suffix}.json"))
        edge_summary = {**edge_summary, **cut_reweight_summary}
    return edge_weights, edge_summary, cut_reweight_summary, cut_labels


def cut_reweighted_contrastive_loss(
    model: CutAwareAutoEncoder,
    anchor_latent: torch.Tensor,
    mixed_latent: torch.Tensor,
    batch_indices: np.ndarray,
    graph,
    cut_labels: np.ndarray | None,
    temperature: float,
    neighbor_positive_weight: float,
    hard_negative_weight: float,
) -> tuple[torch.Tensor, dict]:
    """Masked weighted InfoNCE for the cut-reweighted NeighborMix view.

    Positives are the same cell across original/mixed views plus high-confidence
    same-partition KNN cells. Different-partition non-neighbors are negatives,
    and cross-partition KNN cells receive the hard-negative weight.
    """
    device = anchor_latent.device
    dtype = anchor_latent.dtype
    batch_size = int(anchor_latent.shape[0])
    zero = torch.zeros((), dtype=dtype, device=device)
    empty_stats = {
        "contrast_loss": 0.0,
        "mean_positives": 0.0,
        "mean_negatives": 0.0,
        "mean_hard_negatives": 0.0,
        "valid_batch_fraction": 0.0,
    }
    if batch_size <= 1:
        return zero, empty_stats

    idx = np.asarray(batch_indices, dtype=np.int64)
    neighbor_mask_np = np.zeros((batch_size, batch_size), dtype=bool)
    if graph.indices.shape[1] > 0:
        batch_pos = {int(cell): pos for pos, cell in enumerate(idx.tolist())}
        for row_pos, cell in enumerate(idx.tolist()):
            for nb in graph.indices[int(cell)].tolist():
                col_pos = batch_pos.get(int(nb))
                if col_pos is not None:
                    neighbor_mask_np[row_pos, col_pos] = True

    max_idx = int(idx.max()) if idx.size else -1
    if cut_labels is not None and len(cut_labels) > max_idx:
        part = np.asarray(cut_labels, dtype=np.int64)[idx]
        same_partition_np = part[:, None] == part[None, :]
    else:
        same_partition_np = np.eye(batch_size, dtype=bool)
    eye_np = np.eye(batch_size, dtype=bool)

    pos_mask_np = eye_np | (neighbor_mask_np & same_partition_np)
    hard_neg_np = neighbor_mask_np & (~same_partition_np)
    neg_mask_np = (~same_partition_np) & (~eye_np)

    pos_weight_np = np.where(eye_np, 1.0, np.where(pos_mask_np, float(neighbor_positive_weight), 0.0)).astype(np.float32)
    neg_weight_np = np.where(hard_neg_np, float(hard_negative_weight), np.where(neg_mask_np, 1.0, 0.0)).astype(np.float32)
    candidate_weight_np = pos_weight_np + neg_weight_np
    valid_np = (pos_weight_np.sum(axis=1) > 0.0) & (neg_weight_np.sum(axis=1) > 0.0)
    if not np.any(valid_np):
        return zero, empty_stats

    z_anchor = F.normalize(model.contrast_projection(anchor_latent), dim=1, p=2)
    z_mixed = F.normalize(model.contrast_projection(mixed_latent), dim=1, p=2)
    logits = torch.matmul(z_anchor, z_mixed.t()) / max(float(temperature), 1e-6)
    logits = logits - logits.max(dim=1, keepdim=True).values.detach()
    exp_logits = torch.exp(logits)

    pos_weight = torch.as_tensor(pos_weight_np, dtype=dtype, device=device)
    candidate_weight = torch.as_tensor(candidate_weight_np, dtype=dtype, device=device)
    valid = torch.as_tensor(valid_np, dtype=torch.bool, device=device)
    numerator = (exp_logits * pos_weight).sum(dim=1).clamp_min(1e-12)
    denominator = (exp_logits * candidate_weight).sum(dim=1).clamp_min(1e-12)
    row_loss = -torch.log(numerator / denominator)
    loss = row_loss[valid].mean()

    stats = {
        "contrast_loss": float(loss.detach().cpu()),
        "mean_positives": float(pos_weight_np[valid_np].sum(axis=1).mean()),
        "mean_negatives": float(neg_weight_np[valid_np].sum(axis=1).mean()),
        "mean_hard_negatives": float(hard_neg_np[valid_np].sum(axis=1).mean()),
        "valid_batch_fraction": float(np.mean(valid_np)),
    }
    return loss, stats


def main():
    args = parse_args()
    if args.hidden_dim and not args.hidden_size:
        args.hidden_size = args.hidden_dim
    if args.latent_dim:
        args.hidden_size = args.latent_dim

    if args.variant_name == "canm_diagnostic_only":
        args.cut_weight = 0.0
        args.ot_weight = 0.0
        args.pseudo_weight = 0.0
        args.attention_probe_weight = 0.0
    elif args.variant_name == "canm_cut_ot":
        args.pseudo_weight = 0.0
        args.attention_probe_weight = 0.0
    elif args.variant_name == "canm_mix_plus_cut":
        args.attention_probe_weight = 0.0
    elif args.variant_name == "canm_cut_ot_warm":
        args.pseudo_weight = 0.0
        args.cut_weight = min(float(args.cut_weight), 0.05)
        args.ot_weight = min(float(args.ot_weight), 0.02)
        args.cut_start_epoch = max(int(args.cut_start_epoch), 20)
        args.ot_start_epoch = max(int(args.ot_start_epoch), 20)
        args.attention_probe_weight = 0.0
    elif args.variant_name == "canm_mix_plus_cut_warm":
        args.pseudo_weight = min(float(args.pseudo_weight), 0.1)
        args.cut_weight = min(float(args.cut_weight), 0.05)
        args.ot_weight = min(float(args.ot_weight), 0.02)
        args.cut_start_epoch = max(int(args.cut_start_epoch), 20)
        args.ot_start_epoch = max(int(args.ot_start_epoch), 20)
        args.attention_probe_weight = 0.0
    elif args.variant_name == "canm_cut_reweighted_mix":
        args.pseudo_weight = max(float(args.pseudo_weight), 0.2)
        args.cut_weight = 0.0
        args.ot_weight = 0.0
        args.attention_probe_weight = 0.0
        args.contrast_weight = 0.0
    elif args.variant_name == "canm_cut_reweighted_mix_contrast":
        args.pseudo_weight = max(float(args.pseudo_weight), 0.2)
        args.cut_weight = 0.0
        args.ot_weight = 0.0
        args.attention_probe_weight = 0.0
        args.contrast_weight = max(float(args.contrast_weight), 0.0)
    elif args.variant_name == "canm_gated_cut_mix":
        args.pseudo_weight = max(float(args.pseudo_weight), 0.2)
        args.cut_weight = 0.0
        args.ot_weight = 0.0
        args.gate_prior_weight = min(max(float(args.gate_prior_weight), 0.0), 0.05)
        args.gate_entropy_weight = min(max(float(args.gate_entropy_weight), 0.0), 0.005)
        args.gate_cluster_weight = min(max(float(args.gate_cluster_weight), 0.0), 0.05)
        args.attention_probe_weight = 0.0
    elif args.variant_name == "canm_gated_cut_warm":
        args.pseudo_weight = max(float(args.pseudo_weight), 0.2)
        args.cut_weight = 0.0
        args.ot_weight = 0.0
        args.gate_prior_weight = min(max(float(args.gate_prior_weight), 0.0), 0.05)
        args.gate_entropy_weight = 0.0
        args.gate_cluster_weight = 0.0
        args.gate_start_epoch = max(int(args.gate_start_epoch), 20)
        args.gate_cluster_start_epoch = max(int(args.gate_cluster_start_epoch), 999999)
        args.gate_min = max(float(args.gate_min), 0.2)
        args.attention_probe_weight = 0.0
    elif args.variant_name == "canm_attention_fusion_probe":
        args.pseudo_weight = 0.0
        args.attention_probe_weight = max(float(args.attention_probe_weight), 0.05)

    family.set_seed(args.seed)
    save_dir = Path(ensure_dir(args.save_dir))
    save_json(vars(args), str(save_dir / "args.json"))
    device = family.get_device(args.gpu, args.no_cuda)
    print(f"Using device: {device}", flush=True)

    bundle = family.load_scmae_dataset(
        file_path=args.data_path,
        input_mode=args.input_mode,
        n_top_genes=args.n_top_genes,
        target_sum=args.target_sum,
        scale_input=args.scale_input,
        label_key=args.label_key,
        seed=args.seed,
    )
    data_np = bundle.data
    labels = bundle.labels
    dataset_name = args.dataset_name or Path(args.data_path).stem
    n_clusters = int(args.n_clusters) if args.n_clusters and args.n_clusters > 0 else int(len(np.unique(labels)))
    save_json(bundle.profile, str(save_dir / "dataset_profile.json"))
    save_json(bundle.preprocess_config, str(save_dir / "preprocess_config.json"))
    np.save(save_dir / "gene_names.npy", bundle.gene_names.astype(str))

    graph = build_pca_knn_graph(data_np, k=args.neighbor_k, pca_dim=args.knn_pca_dim, tau=args.tau, seed=args.seed)
    edge_reliability, edge_weights, edge_summary = compute_edge_weights(
        graph,
        mode=args.edge_reliability_mode,
        gamma_sim=args.gamma_sim,
        gamma_mutual=args.gamma_mutual,
        gamma_snn=args.gamma_snn,
        gamma_distance=args.gamma_distance,
        prune_quantile=args.edge_prune_quantile,
    )
    edge_summary = {"edge_reliability_mode": args.edge_reliability_mode, **edge_summary}
    edge_weights, edge_summary, cut_reweight_summary, cut_labels = maybe_apply_cut_reweight(
        args=args,
        save_dir=save_dir,
        graph=graph,
        edge_weights=edge_weights,
        edge_summary=edge_summary,
        n_clusters=n_clusters,
    )
    save_graph_outputs(save_dir, graph, edge_reliability, edge_weights, edge_summary)

    train_loader, eval_loader = make_loaders(data_np, labels, args.batch_size, args.seed)
    model = CutAwareAutoEncoder(
        num_genes=data_np.shape[1],
        n_clusters=n_clusters,
        hidden_size=args.hidden_size,
        dropout=args.dropout,
        masked_data_weight=args.masked_data_weight,
        mask_loss_weight=args.mask_loss_weight,
        contrast_projection_dim=args.contrast_projection_dim if args.variant_name == "canm_cut_reweighted_mix_contrast" else 0,
    ).to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=args.lr)
    gated_mix_enabled = args.variant_name in {"canm_gated_cut_mix", "canm_gated_cut_warm"}
    pseudo_enabled = args.variant_name in {"canm_mix_plus_cut", "canm_mix_plus_cut_warm", "canm_cut_reweighted_mix", "canm_cut_reweighted_mix_contrast", "canm_gated_cut_mix", "canm_gated_cut_warm"} and float(args.pseudo_weight) > 0.0
    cut_enabled = float(args.cut_weight) > 0.0 and int(args.cut_neighbors) > 0
    ot_enabled = float(args.ot_weight) > 0.0
    attention_probe_enabled = float(args.attention_probe_weight) > 0.0
    contrast_enabled = args.variant_name == "canm_cut_reweighted_mix_contrast" and float(args.contrast_weight) > 0.0

    history = {
        "loss": [],
        "scmae_loss": [],
        "pseudo_loss": [],
        "contrast_loss": [],
        "cut_loss": [],
        "ot_loss": [],
        "attention_probe_loss": [],
        "gate_prior_loss": [],
        "gate_entropy_loss": [],
        "gate_cluster_loss": [],
        "real_reconstruction_loss": [],
        "real_mask_loss": [],
        "real_mask_rate": [],
        "mean_mix_delta": [],
        "mean_gate": [],
        "gate_effective_neighbors": [],
        "contrast_mean_positives": [],
        "contrast_mean_negatives": [],
        "contrast_hard_negatives": [],
        "contrast_valid_batch_fraction": [],
        "cluster_mass_min": [],
        "cluster_mass_max": [],
        "assignment_entropy": [],
    }

    for epoch in range(1, max(1, int(args.epochs)) + 1):
        model.train()
        totals = {key: 0.0 for key in history}
        n_batches = 0
        for idx_t, x_cpu, _ in train_loader:
            idx_np = idx_t.numpy().astype(np.int64, copy=False)
            x = x_cpu.to(device)
            x_corrupt, real_mask = family.apply_scmae_noise(x, args.mask_ratio)
            latent, scmae_loss, real_parts = model.loss_mask_weighted(x_corrupt, x, real_mask)
            logits = model.cluster_logits(latent)
            q = torch.softmax(logits / max(float(args.cluster_temperature), 1e-6), dim=1)
            loss = scmae_loss

            pseudo_loss = torch.zeros((), dtype=scmae_loss.dtype, device=device)
            contrast_loss_value = torch.zeros((), dtype=scmae_loss.dtype, device=device)
            contrast_stats = {
                "mean_positives": 0.0,
                "mean_negatives": 0.0,
                "mean_hard_negatives": 0.0,
                "valid_batch_fraction": 0.0,
            }
            mix_info = {"mean_mix_delta": 0.0}
            if pseudo_enabled:
                gate_losses = {
                    "gate_prior_loss": torch.zeros((), dtype=scmae_loss.dtype, device=device),
                    "gate_entropy_loss": torch.zeros((), dtype=scmae_loss.dtype, device=device),
                    "gate_cluster_loss": torch.zeros((), dtype=scmae_loss.dtype, device=device),
                }
                if gated_mix_enabled and epoch >= int(args.gate_start_epoch):
                    x_prime, sample_weight, mix_info, gate_losses = make_gated_neighbor_mixed_batch(
                        data_np=data_np,
                        batch_indices=idx_np,
                        batch_x=x,
                        graph=graph,
                        edge_weights=edge_weights,
                        model=model,
                        anchor_latent=latent,
                        anchor_probs=q,
                        mix_neighbors=args.mix_neighbors,
                        alpha=args.mix_alpha,
                        gate_temperature=args.gate_temperature,
                        gate_min=args.gate_min,
                        cluster_temperature=args.cluster_temperature,
                    )
                else:
                    x_prime, sample_weight, mix_info = make_neighbor_mixed_batch(
                        data_np=data_np,
                        batch_indices=idx_np,
                        batch_x=x,
                        graph=graph,
                        edge_weights=edge_weights,
                        mix_neighbors=args.mix_neighbors,
                        alpha=args.mix_alpha,
                    )
                xp_corrupt, pseudo_mask = family.apply_scmae_noise(x_prime, args.mask_ratio)
                _, pseudo_loss, _ = model.loss_mask_weighted(
                    xp_corrupt,
                    x,
                    pseudo_mask,
                    sample_weight=sample_weight,
                )
                loss = loss + float(args.pseudo_weight) * pseudo_loss
                if contrast_enabled and epoch >= int(args.contrast_start_epoch):
                    mixed_latent = model.encoder(x_prime)
                    contrast_loss_value, contrast_stats = cut_reweighted_contrastive_loss(
                        model=model,
                        anchor_latent=latent,
                        mixed_latent=mixed_latent,
                        batch_indices=idx_np,
                        graph=graph,
                        cut_labels=cut_labels,
                        temperature=args.contrast_temperature,
                        neighbor_positive_weight=args.contrast_neighbor_positive_weight,
                        hard_negative_weight=args.contrast_hard_negative_weight,
                    )
                    loss = loss + float(args.contrast_weight) * contrast_loss_value
                if gated_mix_enabled:
                    loss = loss + float(args.gate_prior_weight) * gate_losses["gate_prior_loss"]
                    loss = loss + float(args.gate_entropy_weight) * gate_losses["gate_entropy_loss"]
                    if epoch >= int(args.gate_cluster_start_epoch):
                        loss = loss + float(args.gate_cluster_weight) * gate_losses["gate_cluster_loss"]
            else:
                gate_losses = {
                    "gate_prior_loss": torch.zeros((), dtype=scmae_loss.dtype, device=device),
                    "gate_entropy_loss": torch.zeros((), dtype=scmae_loss.dtype, device=device),
                    "gate_cluster_loss": torch.zeros((), dtype=scmae_loss.dtype, device=device),
                }

            cut_loss_value = torch.zeros((), dtype=scmae_loss.dtype, device=device)
            cut_stats = {"cut_loss": 0.0, "cluster_mass_min": 0.0, "cluster_mass_max": 0.0, "assignment_entropy": 0.0}
            nb_x, src_rep, edge_w = neighbor_tensors_for_batch(
                data_np=data_np,
                batch_indices=idx_np,
                graph=graph,
                edge_weights=edge_weights,
                max_neighbors=args.cut_neighbors,
                device=device,
            )
            if cut_enabled and epoch >= int(args.cut_start_epoch) and nb_x is not None:
                nb_latent = model.encoder(nb_x)
                q_dst = torch.softmax(model.cluster_logits(nb_latent) / max(float(args.cluster_temperature), 1e-6), dim=1)
                cut_loss_value, cut_stats = graph_cut_loss(q, q_dst, src_rep, edge_w)
                loss = loss + float(args.cut_weight) * cut_loss_value

            ot_loss_value = torch.zeros((), dtype=scmae_loss.dtype, device=device)
            ot_stats = {}
            if ot_enabled and epoch >= int(args.ot_start_epoch):
                ot_loss_value, ot_target, ot_stats = ot_self_training_loss(
                    logits,
                    temperature=args.ot_temperature,
                    iterations=args.ot_iterations,
                )
                loss = loss + float(args.ot_weight) * ot_loss_value

            attention_loss_value = torch.zeros((), dtype=scmae_loss.dtype, device=device)
            if attention_probe_enabled and nb_x is not None:
                with torch.no_grad():
                    nb_latent = model.encoder(nb_x)
                    q_dst = torch.softmax(model.cluster_logits(nb_latent), dim=1)
                    q_graph = torch.zeros_like(q)
                    q_graph.index_add_(0, src_rep, q_dst * edge_w[:, None])
                    q_graph = q_graph / q_graph.sum(dim=1, keepdim=True).clamp_min(1e-8)
                    rel = q_graph.max(dim=1).values
                fused, _ = attention_fusion_probe(q, q_graph, rel)
                attention_loss_value = -(fused.detach() * torch.log(q.clamp_min(1e-8))).sum(dim=1).mean()
                loss = loss + float(args.attention_probe_weight) * attention_loss_value

            if not torch.isfinite(loss):
                print(f"[warn] non-finite loss at epoch={epoch}; skipping batch", flush=True)
                continue

            optimizer.zero_grad(set_to_none=True)
            loss.backward()
            if float(args.grad_clip) > 0.0:
                torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=float(args.grad_clip))
            optimizer.step()

            mean_q = q.detach().mean(dim=0)
            totals["loss"] += float(loss.detach().cpu())
            totals["scmae_loss"] += float(scmae_loss.detach().cpu())
            totals["pseudo_loss"] += float(pseudo_loss.detach().cpu())
            totals["contrast_loss"] += float(contrast_loss_value.detach().cpu())
            totals["cut_loss"] += float(cut_loss_value.detach().cpu())
            totals["ot_loss"] += float(ot_loss_value.detach().cpu())
            totals["attention_probe_loss"] += float(attention_loss_value.detach().cpu())
            totals["gate_prior_loss"] += float(gate_losses["gate_prior_loss"].detach().cpu())
            totals["gate_entropy_loss"] += float(gate_losses["gate_entropy_loss"].detach().cpu())
            totals["gate_cluster_loss"] += float(gate_losses["gate_cluster_loss"].detach().cpu())
            totals["real_reconstruction_loss"] += float(real_parts["reconstruction_loss"].cpu())
            totals["real_mask_loss"] += float(real_parts["mask_loss"].cpu())
            totals["real_mask_rate"] += float(real_mask.mean().detach().cpu())
            totals["mean_mix_delta"] += float(mix_info["mean_mix_delta"])
            totals["mean_gate"] += float(mix_info.get("mean_gate", 0.0))
            totals["gate_effective_neighbors"] += float(mix_info.get("gate_effective_neighbors", 0.0))
            totals["contrast_mean_positives"] += float(contrast_stats.get("mean_positives", 0.0))
            totals["contrast_mean_negatives"] += float(contrast_stats.get("mean_negatives", 0.0))
            totals["contrast_hard_negatives"] += float(contrast_stats.get("mean_hard_negatives", 0.0))
            totals["contrast_valid_batch_fraction"] += float(contrast_stats.get("valid_batch_fraction", 0.0))
            totals["cluster_mass_min"] += float(cut_stats.get("cluster_mass_min", mean_q.min().detach().cpu()))
            totals["cluster_mass_max"] += float(cut_stats.get("cluster_mass_max", mean_q.max().detach().cpu()))
            totals["assignment_entropy"] += float(cut_stats.get("assignment_entropy", (-(q * q.clamp_min(1e-8).log()).sum(dim=1).mean()).detach().cpu()))
            n_batches += 1

        for key, value in totals.items():
            history[key].append(value / max(1, n_batches))
        if epoch == 1 or epoch == args.epochs or epoch % 10 == 0:
            print(
                f"Epoch {epoch:03d}/{args.epochs} loss={history['loss'][-1]:.4f} "
                f"scmae={history['scmae_loss'][-1]:.4f} pseudo={history['pseudo_loss'][-1]:.4f} "
                f"contrast={history['contrast_loss'][-1]:.4f} cut={history['cut_loss'][-1]:.4f} "
                f"ot={history['ot_loss'][-1]:.4f} gate={history['mean_gate'][-1]:.3f} "
                f"mass=[{history['cluster_mass_min'][-1]:.3f},{history['cluster_mass_max'][-1]:.3f}]",
                flush=True,
            )

        if (
            args.graph_refresh_interval > 0
            and epoch >= int(args.graph_refresh_start_epoch)
            and epoch % int(args.graph_refresh_interval) == 0
            and epoch != args.epochs
        ):
            embedding, _ = family.extract_embedding(model, eval_loader, device)
            graph = build_embedding_knn_graph(embedding, k=args.neighbor_k, tau=args.tau, source=f"embedding_epoch_{epoch}")
            edge_reliability, edge_weights, edge_summary = compute_edge_weights(
                graph,
                mode=args.edge_reliability_mode,
                gamma_sim=args.gamma_sim,
                gamma_mutual=args.gamma_mutual,
                gamma_snn=args.gamma_snn,
                gamma_distance=args.gamma_distance,
                prune_quantile=args.edge_prune_quantile,
            )
            edge_summary = {"edge_reliability_mode": args.edge_reliability_mode, **edge_summary}
            edge_weights, edge_summary, cut_reweight_summary, cut_labels = maybe_apply_cut_reweight(
                args=args,
                save_dir=save_dir,
                graph=graph,
                edge_weights=edge_weights,
                edge_summary=edge_summary,
                n_clusters=n_clusters,
                suffix=f"_epoch_{epoch}",
            )

    embedding, labels_out = family.extract_embedding(model, eval_loader, device)
    q_all, _ = extract_cluster_probs(model, eval_loader, device, temperature=args.cluster_temperature)
    cut_diag = graph_cut_diagnostics(q_all, graph, edge_weights)
    sim_diag = pairwise_similarity_digest(embedding, seed=args.seed)
    np.save(save_dir / "embedding_final.npy", embedding.astype(np.float32))
    np.save(save_dir / "embeddings_base.npy", embedding.astype(np.float32))
    np.save(save_dir / "labels.npy", labels_out.astype(np.int64))
    np.save(save_dir / "cluster_probs.npy", q_all.astype(np.float32))
    family.save_embedding_h5(save_dir / "embedding.h5", embedding, labels_out)
    save_json(history, str(save_dir / "training_history.json"))
    save_json(cut_diag, str(save_dir / "cut_diagnostics.json"))
    save_json(sim_diag, str(save_dir / "embedding_similarity_diagnostics.json"))
    contrast_diag = {
        "contrast_enabled": bool(contrast_enabled),
        "contrast_weight": float(args.contrast_weight),
        "contrast_temperature": float(args.contrast_temperature),
        "contrast_start_epoch": int(args.contrast_start_epoch),
        "contrast_neighbor_positive_weight": float(args.contrast_neighbor_positive_weight),
        "contrast_hard_negative_weight": float(args.contrast_hard_negative_weight),
        "contrast_projection_dim": int(args.contrast_projection_dim) if contrast_enabled else 0,
        "mean_contrast_loss": float(np.mean(history["contrast_loss"])) if history["contrast_loss"] else 0.0,
        "final_contrast_loss": float(history["contrast_loss"][-1]) if history["contrast_loss"] else 0.0,
        "mean_positives": float(np.mean(history["contrast_mean_positives"])) if history["contrast_mean_positives"] else 0.0,
        "mean_negatives": float(np.mean(history["contrast_mean_negatives"])) if history["contrast_mean_negatives"] else 0.0,
        "mean_hard_negatives": float(np.mean(history["contrast_hard_negatives"])) if history["contrast_hard_negatives"] else 0.0,
        "mean_valid_batch_fraction": float(np.mean(history["contrast_valid_batch_fraction"])) if history["contrast_valid_batch_fraction"] else 0.0,
    }
    save_json(contrast_diag, str(save_dir / "contrast_diagnostics.json"))
    save_graph_outputs(save_dir, graph, edge_reliability, edge_weights, edge_summary)
    torch.save(
        {
            "model_state": model.state_dict(),
            "args": vars(args),
            "gene_names": bundle.gene_names.astype(str),
            "neighbor_graph_profile": graph.profile,
            "edge_weight_summary": edge_summary,
            "cut_diagnostics": cut_diag,
            "embedding_similarity_diagnostics": sim_diag,
            "contrast_diagnostics": contrast_diag,
        },
        save_dir / "model.pt",
    )

    result = None
    pred = None
    if not args.skip_eval:
        result = family.write_kmeans_known_k_outputs(
            output_dir=save_dir,
            dataset=dataset_name,
            method=args.method_name,
            seed=args.seed,
            embedding=embedding,
            labels=labels_out,
            n_clusters=n_clusters,
            extra={
                "variant": args.variant_name,
                "cut_weight": float(args.cut_weight),
                "ot_weight": float(args.ot_weight),
                "pseudo_weight": float(args.pseudo_weight),
                "contrast_enabled": bool(contrast_enabled),
                "contrast_weight": float(args.contrast_weight),
                "contrast_temperature": float(args.contrast_temperature),
                "cut_cross_weight": float(args.cut_cross_weight),
                "gate_prior_weight": float(args.gate_prior_weight),
                "gate_entropy_weight": float(args.gate_entropy_weight),
                "gate_cluster_weight": float(args.gate_cluster_weight),
                "gate_start_epoch": int(args.gate_start_epoch),
                "neighbor_k": int(args.neighbor_k),
                "edge_prune_quantile": float(args.edge_prune_quantile),
                "graph_refresh_interval": int(args.graph_refresh_interval),
            },
        )
        pred = result["preds"]["kmeans_known_k"]
        save_json(result["fixed"], str(save_dir / "metrics.json"))
    else:
        save_json({}, str(save_dir / "metrics.json"))

    pd.DataFrame([embedding_geometry(embedding, labels_out)]).to_csv(save_dir / "embedding_geometry_summary.csv", index=False)
    if pred is not None:
        pct = per_cell_type_metrics(labels_out, pred)
        pct.to_csv(save_dir / "per_cell_type_metrics.csv", index=False)
        mapped = mapped_predictions(labels_out, pred)
        pd.crosstab(pd.Series(labels_out, name="true"), pd.Series(pred, name="pred")).to_csv(save_dir / "confusion_matrix_raw.csv")
        pd.crosstab(pd.Series(labels_out, name="true"), pd.Series(mapped, name="mapped")).to_csv(save_dir / "confusion_matrix_mapped.csv")
        rare = pct.groupby("is_rare_lt_50")[["precision", "recall", "f1"]].mean().reset_index()
        rare.to_csv(save_dir / "rare_cell_effect_summary.csv", index=False)

    summary = {
        "dataset": dataset_name,
        "method": args.method_name,
        "variant": args.variant_name,
        "seed": int(args.seed),
        "n_cells": int(data_np.shape[0]),
        "n_genes": int(data_np.shape[1]),
        "n_clusters": int(n_clusters),
        "pseudo_enabled": bool(pseudo_enabled),
        "gated_mix_enabled": bool(gated_mix_enabled),
        "cut_enabled": bool(cut_enabled),
        "ot_enabled": bool(ot_enabled),
        "attention_probe_enabled": bool(attention_probe_enabled),
        "contrast_enabled": bool(contrast_enabled),
        "contrast_weight": float(args.contrast_weight),
        "contrast_temperature": float(args.contrast_temperature),
        "contrast_start_epoch": int(args.contrast_start_epoch),
        "contrast_neighbor_positive_weight": float(args.contrast_neighbor_positive_weight),
        "contrast_hard_negative_weight": float(args.contrast_hard_negative_weight),
        "contrast_projection_dim": int(args.contrast_projection_dim) if contrast_enabled else 0,
        "edge_weight_summary": edge_summary,
        "cut_reweight_summary": cut_reweight_summary,
        "cut_diagnostics": cut_diag,
        "embedding_similarity_diagnostics": sim_diag,
        "contrast_diagnostics": contrast_diag,
        "fixed_metrics": result["fixed"] if result is not None else {},
        "label_leakage": False,
    }
    save_json(summary, str(save_dir / "summary.json"))

    if not args.no_save_h5ad:
        bundle.adata.obsm["X_cutaware_neighbormix_scmae"] = embedding
        bundle.adata.uns["cutaware_neighbormix_scmae"] = summary
        sanitize_anndata_for_write(bundle.adata)
        bundle.adata.write_h5ad(save_dir / "adata_cutaware_neighbormix_scmae.h5ad", compression="gzip")

    print(f"Results saved to: {save_dir}", flush=True)


if __name__ == "__main__":
    main()
