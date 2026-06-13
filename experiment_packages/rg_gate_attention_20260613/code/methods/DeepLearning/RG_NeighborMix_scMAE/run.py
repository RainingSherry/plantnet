#!/usr/bin/env python
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from torch.utils.data import DataLoader

CURRENT_DIR = Path(__file__).resolve().parent
ROOT = CURRENT_DIR.parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from methods.DeepLearning import scMAE_family as family
from methods.DeepLearning.RG_NeighborMix_scMAE.diagnostics import (
    embedding_geometry,
    mapped_predictions,
    per_cell_type_metrics,
)
from methods.DeepLearning.RG_NeighborMix_scMAE.mixing import compute_node_gate, make_pseudo_batch
from methods.DeepLearning.RG_NeighborMix_scMAE.model import AutoEncoder
from methods.DeepLearning.RG_NeighborMix_scMAE.neighbor_graph import (
    build_far_neighbors,
    build_pca_knn_graph,
    build_random_neighbors,
    compute_edge_reliability,
)
from methods.shared_utils import ensure_dir, save_json, sanitize_anndata_for_write


def str2bool(value):
    return family.str2bool(value)


def parse_args():
    parser = argparse.ArgumentParser(description="Reliability-Gated NeighborMix-scMAE runner")
    parser.add_argument("--data_path", required=True)
    parser.add_argument("--save_dir", required=True)
    parser.add_argument("--dataset_name", default=None)
    parser.add_argument("--method_name", default="RG_NeighborMix_scMAE")
    parser.add_argument("--variant_name", default="rg_nm_v1_reliability")
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
    parser.add_argument("--mix_mode", default="reliability", choices=["none", "fixed", "mutual", "reliability", "random", "far", "attention"])
    parser.add_argument("--allow_attention_phase2", type=str2bool, default=False)
    parser.add_argument("--gate_mode", default="topology", choices=["none", "constant", "topology", "learned"])
    parser.add_argument("--edge_reliability_mode", default="sim_mutual_snn_distance")
    parser.add_argument("--neighbor_k", type=int, default=10)
    parser.add_argument("--mix_neighbors", type=int, default=4)
    parser.add_argument("--tau", type=float, default=0.2)
    parser.add_argument("--knn_pca_dim", type=int, default=50)
    parser.add_argument("--pseudo_weight", type=float, default=0.3)
    parser.add_argument("--gate_max", type=float, default=0.15)
    parser.add_argument("--gate_min", type=float, default=0.0)
    parser.add_argument("--gamma_sim", type=float, default=1.0)
    parser.add_argument("--gamma_mutual", type=float, default=1.0)
    parser.add_argument("--gamma_snn", type=float, default=1.0)
    parser.add_argument("--gamma_distance", type=float, default=1.0)
    parser.add_argument("--beta_mutual", type=float, default=1.0)
    parser.add_argument("--beta_snn", type=float, default=1.0)
    parser.add_argument("--beta_perturb", type=float, default=2.0)
    parser.add_argument("--beta_uncertainty", type=float, default=1.0)
    parser.add_argument("--save_neighbor_diagnostics", type=str2bool, default=True)
    parser.add_argument("--save_gate_diagnostics", type=str2bool, default=True)
    parser.add_argument("--save_embedding_geometry", type=str2bool, default=True)
    parser.add_argument("--save_per_cell_type_metrics", type=str2bool, default=True)
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


def main():
    args = parse_args()
    if args.mix_mode == "attention" and not args.allow_attention_phase2:
        raise ValueError("--mix_mode attention is reserved for phase 2. Pass --allow_attention_phase2 true only after phase-1 criteria pass.")
    if args.hidden_dim and not args.hidden_size:
        args.hidden_size = args.hidden_dim
    if args.latent_dim:
        args.hidden_size = args.latent_dim
    if args.mix_mode == "none":
        args.gate_mode = "none"
        args.pseudo_weight = 0.0
        args.neighbor_k = 0
    if args.mix_mode == "fixed":
        args.gate_mode = "constant"
        args.edge_reliability_mode = "none"

    family.set_seed(args.seed)
    rng = np.random.default_rng(args.seed + 3089)
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
    edge_reliability, edge_weights, edge_summary = compute_edge_reliability(
        graph,
        mode=args.edge_reliability_mode,
        gamma_sim=args.gamma_sim,
        gamma_mutual=args.gamma_mutual,
        gamma_snn=args.gamma_snn,
        gamma_distance=args.gamma_distance,
    )
    node_gate, normalized_gate, gate_summary = compute_node_gate(
        graph,
        edge_weights=edge_weights,
        gate_mode=args.gate_mode,
        gate_min=args.gate_min,
        gate_max=args.gate_max,
        beta_mutual=args.beta_mutual,
        beta_snn=args.beta_snn,
        beta_perturb=args.beta_perturb,
        beta_uncertainty=args.beta_uncertainty,
        uncertainty=None,
    )
    random_neighbors = build_random_neighbors(data_np.shape[0], max(1, min(args.mix_neighbors, max(1, data_np.shape[0] - 1))), rng, graph.indices)
    far_neighbors = build_far_neighbors(graph.embedding, max(1, min(args.mix_neighbors, max(1, data_np.shape[0] - 1))), rng)
    perturb_proxy = (1.0 - np.sum(graph.probs * graph.similarity, axis=1)).astype(np.float32) if graph.probs.size else np.zeros(data_np.shape[0], dtype=np.float32)

    if args.save_neighbor_diagnostics:
        np.save(save_dir / "neighbor_indices.npy", graph.indices)
        np.save(save_dir / "neighbor_base_probs.npy", graph.probs)
        np.save(save_dir / "neighbor_similarity.npy", graph.similarity)
        np.save(save_dir / "neighbor_distance.npy", graph.distance)
        save_json(graph.profile, str(save_dir / "neighbor_graph_profile.json"))
        save_json(graph.profile, str(save_dir / "train_knn_diagnostics.json"))
    np.save(save_dir / "edge_reliability.npy", edge_reliability)
    np.save(save_dir / "node_gate.npy", node_gate)
    np.save(save_dir / "pseudo_perturbation.npy", perturb_proxy)
    save_json({"edge_reliability_mode": args.edge_reliability_mode, **edge_summary}, str(save_dir / "edge_weight_summary.json"))
    save_json({"edge_reliability_mode": args.edge_reliability_mode, **edge_summary}, str(save_dir / "neighbor_reliability_summary.json"))
    save_json(gate_summary, str(save_dir / "gate_summary.json"))
    save_json(
        {
            "mean_pseudo_perturbation": float(np.mean(perturb_proxy)) if perturb_proxy.size else 0.0,
            "p95_pseudo_perturbation": float(np.percentile(perturb_proxy, 95)) if perturb_proxy.size else 0.0,
        },
        str(save_dir / "pseudo_perturbation_summary.json"),
    )

    train_loader, eval_loader = make_loaders(data_np, labels, args.batch_size, args.seed)
    model = AutoEncoder(
        num_genes=data_np.shape[1],
        hidden_size=args.hidden_size,
        dropout=args.dropout,
        masked_data_weight=args.masked_data_weight,
        mask_loss_weight=args.mask_loss_weight,
    ).to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=args.lr)
    pseudo_enabled = args.mix_mode != "none" and float(args.pseudo_weight) > 0 and int(args.mix_neighbors) > 0

    history = {
        "loss": [],
        "real_loss": [],
        "real_reconstruction_loss": [],
        "real_mask_loss": [],
        "pseudo_loss": [],
        "pseudo_reconstruction_loss": [],
        "pseudo_mask_loss": [],
        "mean_node_gate": [],
        "mean_pseudo_perturbation": [],
        "real_mask_rate": [],
        "pseudo_mask_rate": [],
        "mix_mode": args.mix_mode,
        "pseudo_enabled": bool(pseudo_enabled),
    }

    for epoch in range(1, max(1, int(args.epochs)) + 1):
        model.train()
        totals = {key: 0.0 for key in history if isinstance(history[key], list)}
        n_batches = 0
        for idx_t, x_cpu, _ in train_loader:
            idx_np = idx_t.numpy().astype(np.int64, copy=False)
            x = x_cpu.to(device)
            x_corrupt, real_mask = family.apply_scmae_noise(x, args.mask_ratio)
            _, real_loss, real_parts = model.loss_mask_weighted(x_corrupt, x, real_mask, mask_loss_scale=1.0)
            loss = real_loss
            pseudo_loss = torch.zeros((), dtype=real_loss.dtype, device=device)
            pseudo_parts = {"reconstruction_loss": pseudo_loss, "mask_loss": pseudo_loss, "mask_positive_rate": pseudo_loss}
            mix_info = {"mean_node_gate": 0.0, "mean_perturb_norm": 0.0}
            if pseudo_enabled:
                x_prime, sample_weight, mix_info = make_pseudo_batch(
                    data_np=data_np,
                    batch_indices=idx_np,
                    batch_x=x,
                    mix_mode=args.mix_mode,
                    graph=graph,
                    edge_weights=edge_weights,
                    node_gate=node_gate,
                    mix_neighbors=args.mix_neighbors,
                    rng=rng,
                    random_neighbors=random_neighbors,
                    far_neighbors=far_neighbors,
                )
                xp_corrupt, pseudo_mask = family.apply_scmae_noise(x_prime, args.mask_ratio)
                _, pseudo_loss, pseudo_parts = model.loss_mask_weighted(
                    xp_corrupt,
                    x,
                    pseudo_mask,
                    sample_weight=sample_weight,
                    mask_loss_scale=1.0,
                )
                loss = loss + float(args.pseudo_weight) * pseudo_loss

            optimizer.zero_grad(set_to_none=True)
            loss.backward()
            optimizer.step()

            totals["loss"] += float(loss.detach().cpu())
            totals["real_loss"] += float(real_loss.detach().cpu())
            totals["real_reconstruction_loss"] += float(real_parts["reconstruction_loss"].cpu())
            totals["real_mask_loss"] += float(real_parts["mask_loss"].cpu())
            totals["pseudo_loss"] += float(pseudo_loss.detach().cpu())
            totals["pseudo_reconstruction_loss"] += float(pseudo_parts["reconstruction_loss"].cpu())
            totals["pseudo_mask_loss"] += float(pseudo_parts["mask_loss"].cpu())
            totals["mean_node_gate"] += float(mix_info["mean_node_gate"])
            totals["mean_pseudo_perturbation"] += float(mix_info["mean_perturb_norm"])
            totals["real_mask_rate"] += float(real_mask.mean().detach().cpu())
            totals["pseudo_mask_rate"] += float(pseudo_parts["mask_positive_rate"].cpu()) if pseudo_enabled else 0.0
            n_batches += 1

        for key, value in totals.items():
            history[key].append(value / max(1, n_batches))
        if epoch == 1 or epoch == args.epochs or epoch % 10 == 0:
            print(
                f"Epoch {epoch:03d}/{args.epochs} loss={history['loss'][-1]:.4f} "
                f"real={history['real_loss'][-1]:.4f} pseudo={history['pseudo_loss'][-1]:.4f} "
                f"gate={history['mean_node_gate'][-1]:.4f}",
                flush=True,
            )

    embedding, labels_out = family.extract_embedding(model, eval_loader, device)
    np.save(save_dir / "embedding_final.npy", embedding.astype(np.float32))
    np.save(save_dir / "embeddings_base.npy", embedding.astype(np.float32))
    np.save(save_dir / "labels.npy", labels_out.astype(np.int64))
    family.save_embedding_h5(save_dir / "embedding.h5", embedding, labels_out)
    save_json(history, str(save_dir / "training_history.json"))
    torch.save(
        {
            "model_state": model.state_dict(),
            "args": vars(args),
            "gene_names": bundle.gene_names.astype(str),
            "neighbor_graph_profile": graph.profile,
            "edge_weight_summary": edge_summary,
            "gate_summary": gate_summary,
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
                "mix_mode": args.mix_mode,
                "gate_mode": args.gate_mode,
                "edge_reliability_mode": args.edge_reliability_mode,
                "gate_max": float(args.gate_max),
                "pseudo_weight": float(args.pseudo_weight),
            },
        )
        pred = result["preds"]["kmeans_known_k"]
        save_json(result["fixed"], str(save_dir / "metrics.json"))
    else:
        save_json({}, str(save_dir / "metrics.json"))

    if args.save_embedding_geometry:
        pd.DataFrame([embedding_geometry(embedding, labels_out)]).to_csv(save_dir / "embedding_geometry_summary.csv", index=False)
    if pred is not None and args.save_per_cell_type_metrics:
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
        "mix_mode": args.mix_mode,
        "pseudo_enabled": bool(pseudo_enabled),
        "edge_weight_summary": edge_summary,
        "gate_summary": gate_summary,
        "fixed_metrics": result["fixed"] if result is not None else {},
        "label_leakage": False,
    }
    save_json(summary, str(save_dir / "summary.json"))

    if not args.no_save_h5ad:
        bundle.adata.obsm["X_rg_neighbormix_scmae"] = embedding
        bundle.adata.uns["rg_neighbormix_scmae"] = summary
        sanitize_anndata_for_write(bundle.adata)
        bundle.adata.write_h5ad(save_dir / "adata_rg_neighbormix_scmae.h5ad", compression="gzip")

    print(f"Results saved to: {save_dir}", flush=True)


if __name__ == "__main__":
    main()
