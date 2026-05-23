import argparse
import os
import sys

import numpy as np
import torch

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
PARENT = os.path.dirname(CURRENT_DIR)
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(PARENT)))
METHODS_DIR = os.path.join(ROOT, "methods")
for path in [PARENT, ROOT, METHODS_DIR, CURRENT_DIR]:
    if path not in sys.path:
        sys.path.insert(0, path)

from methods.utils import save as save_benchmark
from maskdiffusion.data import load_sc_dataset, make_dataloader
from maskdiffusion.eval import cluster_and_evaluate
from maskdiffusion.train import ScSpade, extract_embeddings, initialize_cluster_centers, train_epoch
from maskdiffusion.utils import ensure_dir, get_device, save_json, set_seed


def parse_args():
    parser = argparse.ArgumentParser(description="Codex ScSpade: support-masked latent diffusion AE")
    parser.add_argument("--data_path", type=str, required=True)
    parser.add_argument("--save_dir", type=str, required=True)
    parser.add_argument("--input_mode", type=str, default="auto", choices=["auto", "raw", "log1p"])
    parser.add_argument("--n_top_genes", type=int, default=2000)
    parser.add_argument("--n_clusters", type=int, default=0)
    parser.add_argument("--latent_dim", type=int, default=32)
    parser.add_argument("--mask_hidden_dims", type=str, default="512,256,128")
    parser.add_argument("--diffusion_hidden_dims", type=str, default="1024,512")
    parser.add_argument("--diffusion_steps", type=int, default=100)
    parser.add_argument("--dropout", type=float, default=0.15)
    parser.add_argument("--epochs", type=int, default=150)
    parser.add_argument("--warmup_epochs", type=int, default=30)
    parser.add_argument("--diffusion_ramp_epochs", type=int, default=50)
    parser.add_argument("--cluster_warmup_epochs", type=int, default=0)
    parser.add_argument("--batch_size", type=int, default=256)
    parser.add_argument("--lr", type=float, default=1e-3)
    parser.add_argument("--weight_decay", type=float, default=1e-4)
    parser.add_argument("--mask_loss_weight", type=float, default=0.2)
    parser.add_argument("--recon_loss_weight", type=float, default=0.8)
    parser.add_argument("--diffusion_loss_weight", type=float, default=0.05)
    parser.add_argument("--diffusion_warmup_weight", type=float, default=0.05)
    parser.add_argument("--cluster_loss_weight", type=float, default=0.0)
    parser.add_argument("--mask_coupling", type=str, default="weighted_observed", choices=["weighted_observed", "prob", "observed"])
    parser.add_argument("--diffusion_start_frac", type=float, default=0.35)
    parser.add_argument("--eval_interval", type=int, default=10)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--gpu", type=int, default=0)
    parser.add_argument("--no_cuda", action="store_true")
    return parser.parse_args()


def _parse_dims(value: str):
    return [int(item) for item in value.split(",") if item.strip()]


def _save_embedding_outputs(save_dir: str, name: str, embedding: np.ndarray, labels: np.ndarray, result: dict) -> None:
    np.save(os.path.join(save_dir, f"embeddings_{name}.npy"), embedding.astype(np.float32))
    np.save(os.path.join(save_dir, f"pred_labels_{name}.npy"), result["pred_labels"])
    np.save(os.path.join(save_dir, f"pred_labels_{name}_mapped.npy"), result["pred_labels_mapped"])
    save_json(result["metrics"], os.path.join(save_dir, f"summary_{name}.json"))

    bench_dir = ensure_dir(os.path.join(save_dir, f"benchmark_{name}"))
    save_benchmark(bench_dir, labels, result["pred_labels"], epoch="final", embedding=embedding)


def main():
    args = parse_args()
    set_seed(args.seed)
    save_dir = ensure_dir(args.save_dir)
    save_json(vars(args), os.path.join(save_dir, "args.json"))

    device = get_device(gpu=args.gpu, no_cuda=args.no_cuda)
    print(f"Using device: {device}")
    if device.type == "cuda" and args.gpu == 0:
        print("Warning: --gpu 0 was requested. Use --gpu 2/3/4/5 for the requested non-zero GPU run.")

    print("=" * 70)
    print("Step 1: loading data and preserving support/value separation")
    print("=" * 70)
    bundle = load_sc_dataset(
        args.data_path,
        input_mode=args.input_mode,
        n_top_genes=args.n_top_genes,
    )
    if bundle.labels is None:
        raise ValueError("No labels found in h5ad obs; benchmark evaluation requires cell-type labels.")
    n_clusters = args.n_clusters if args.n_clusters > 0 else len(np.unique(bundle.labels))
    print(f"Input mode: {bundle.input_mode}")
    print(f"Cells: {bundle.values.shape[0]} | genes: {bundle.values.shape[1]} | clusters: {n_clusters}")
    print(f"Observed support density: {float(bundle.support.mean()):.4f}")

    train_loader = make_dataloader(
        bundle.values,
        bundle.support,
        bundle.labels,
        batch_size=args.batch_size,
        shuffle=True,
        drop_last=False,
    )
    eval_loader = make_dataloader(
        bundle.values,
        bundle.support,
        bundle.labels,
        batch_size=args.batch_size,
        shuffle=False,
        drop_last=False,
    )

    print("=" * 70)
    print("Step 2: building model")
    print("=" * 70)
    model = ScSpade(
        num_genes=bundle.values.shape[1],
        n_clusters=n_clusters,
        latent_dim=args.latent_dim,
        mask_hidden_dims=_parse_dims(args.mask_hidden_dims),
        diffusion_hidden_dims=_parse_dims(args.diffusion_hidden_dims),
        diffusion_steps=args.diffusion_steps,
        dropout=args.dropout,
        mask_coupling=args.mask_coupling,
    ).to(device)
    optimizer = torch.optim.AdamW(model.parameters(), lr=args.lr, weight_decay=args.weight_decay)
    scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=max(args.epochs, 1))
    print(f"Parameters: {sum(p.numel() for p in model.parameters()):,}")

    print("=" * 70)
    print("Step 3: training with strong AE warmup, delayed diffusion, optional late cluster loss")
    print("=" * 70)
    best_state = None
    best_direct_nmi = -1.0
    cluster_initialized = False
    cluster_start_epoch = args.warmup_epochs + args.diffusion_ramp_epochs + args.cluster_warmup_epochs

    for epoch in range(args.epochs):
        if (
            args.cluster_loss_weight > 0
            and not cluster_initialized
            and epoch >= cluster_start_epoch
        ):
            initialize_cluster_centers(model, eval_loader, device, n_clusters=n_clusters, random_state=args.seed)
            cluster_initialized = True
            print(f"Initialized cluster centers with KMeans at epoch {epoch + 1}")

        train_metrics = train_epoch(
            model,
            train_loader,
            optimizer,
            device,
            epoch,
            warmup_epochs=args.warmup_epochs,
            diffusion_ramp_epochs=args.diffusion_ramp_epochs,
            cluster_warmup_epochs=args.cluster_warmup_epochs,
            mask_weight=args.mask_loss_weight,
            diffusion_weight=args.diffusion_loss_weight,
            recon_weight=args.recon_loss_weight,
            cluster_weight=args.cluster_loss_weight,
            diffusion_warmup_weight=args.diffusion_warmup_weight,
        )
        scheduler.step()

        should_eval = (epoch + 1) % args.eval_interval == 0 or epoch == args.epochs - 1
        if should_eval:
            extracted = extract_embeddings(
                model,
                eval_loader,
                device,
                diffusion_start_frac=args.diffusion_start_frac,
                return_masks=False,
            )
            direct_result = cluster_and_evaluate(extracted["direct"], extracted["labels"], n_clusters, seed=args.seed)
            direct_nmi = direct_result["metrics"]["nmi"]
            if direct_nmi > best_direct_nmi:
                best_direct_nmi = direct_nmi
                best_state = {key: value.detach().cpu().clone() for key, value in model.state_dict().items()}
                torch.save(best_state, os.path.join(save_dir, "best_model_direct.pt"))
            print(
                f"Epoch {epoch + 1:03d}/{args.epochs} [{train_metrics['phase']}] "
                f"loss={train_metrics['loss']:.4f} recon={train_metrics['recon_loss']:.4f} "
                f"mask={train_metrics['mask_loss']:.4f} diff={train_metrics['diffusion_loss']:.4f} "
                f"cluster={train_metrics['cluster_loss']:.4f} direct_nmi={direct_nmi:.4f} "
                f"direct_ari={direct_result['metrics']['ari']:.4f}"
            )

    if best_state is not None:
        model.load_state_dict(best_state)

    print("=" * 70)
    print("Step 4: final direct and diffusion embedding evaluation")
    print("=" * 70)
    extracted = extract_embeddings(
        model,
        eval_loader,
        device,
        diffusion_start_frac=args.diffusion_start_frac,
        return_masks=True,
    )
    labels = extracted["labels"]
    direct_embedding = extracted["direct"]
    diffusion_embedding = extracted["diffusion"]
    mask_probs = extracted["mask_probs"]

    direct_result = cluster_and_evaluate(direct_embedding, labels, n_clusters, seed=args.seed)
    diffusion_result = cluster_and_evaluate(diffusion_embedding, labels, n_clusters, seed=args.seed)
    results = {
        "direct": direct_result["metrics"],
        "diffusion": diffusion_result["metrics"],
        "best_direct_nmi": best_direct_nmi,
        "n_cells": int(bundle.values.shape[0]),
        "n_genes": int(bundle.values.shape[1]),
        "n_clusters": int(n_clusters),
        "input_mode": bundle.input_mode,
        "support_density": float(bundle.support.mean()),
    }
    save_json(results, os.path.join(save_dir, "summary.json"))

    _save_embedding_outputs(save_dir, "direct", direct_embedding, labels, direct_result)
    _save_embedding_outputs(save_dir, "diffusion", diffusion_embedding, labels, diffusion_result)
    np.save(os.path.join(save_dir, "labels.npy"), labels)
    np.save(os.path.join(save_dir, "mask_probs.npy"), mask_probs.astype(np.float32))

    primary_name = "diffusion" if diffusion_result["metrics"]["nmi"] >= direct_result["metrics"]["nmi"] else "direct"
    primary_embedding = diffusion_embedding if primary_name == "diffusion" else direct_embedding
    primary_pred = diffusion_result["pred_labels"] if primary_name == "diffusion" else direct_result["pred_labels"]
    np.save(os.path.join(save_dir, "embedding_final.npy"), primary_embedding.astype(np.float32))
    np.save(os.path.join(save_dir, "pred_labels.npy"), primary_pred.astype(np.int64))

    print("Final direct metrics:", direct_result["metrics"])
    print("Final diffusion metrics:", diffusion_result["metrics"])
    print(f"Primary embedding selected by NMI: {primary_name}")
    print(f"Results saved to: {save_dir}")


if __name__ == "__main__":
    main()
