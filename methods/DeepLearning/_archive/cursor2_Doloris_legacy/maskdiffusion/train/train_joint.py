"""ScSpade joint training with staged loss schedule.

Training strategy follows the design doc: stable-first, not competition-first.

PHASE 1 — Warmup (epochs 0 to warmup_epochs):
    Only mask loss + reconstruction loss active.
    The model learns basic expression structure without interference from
    diffusion noise or random cluster centers.

PHASE 2 — Diffusion warmup (epochs warmup_epochs to warmup_epochs + 50):
    Mask loss + reconstruction loss + very weak diffusion loss.
    The denoiser learns latent geometry gently. Cluster loss stays OFF.

PHASE 3 — Full training (after warmup_epochs + 50):
    All losses active with configured weights.
    Cluster loss uses KMeans-initialized centers.

The staged schedule prevents the three failure modes identified in the analysis:
    1. Reconstruction signal destroyed by equal-weight multi-objective competition
    2. Latent space corrupted by random cluster centers
    3. Denoiser destabilized by noisy encoder output
"""

import os
import time
from typing import Optional, Tuple, Callable

import numpy as np
import torch
import torch.nn as nn
from tqdm import tqdm

from data import build_dataloader, SCDatasetBundle
from models import ScSpade


# ── Training epoch ─────────────────────────────────────────────────────────────


def train_scspade_epoch(
    model: nn.Module,
    dataloader: torch.utils.data.DataLoader,
    optimizer: torch.optim.Optimizer,
    epoch: int,
    *,
    warmup_epochs: int = 30,
    mask_loss_weight: float = 0.2,
    recon_loss_weight: float = 0.8,
    diffusion_loss_weight: float = 0.1,
    cluster_loss_weight: float = 0.0,
    cluster_warmup_epochs: int = 0,
    device: torch.device = None,
    progress_bar: bool = True,
) -> dict:
    """Run one training epoch with staged loss weights.

    The staged schedule:
        Epoch < warmup_epochs:
            effective_cluster_weight = 0.0
            effective_diffusion_weight = 0.0
            Only mask + reconstruction active.

        Epoch < warmup_epochs + 50:
            effective_diffusion_weight = min(diffusion_weight, 0.05)
            effective_cluster_weight = 0.0

        Epoch >= warmup_epochs + 50:
            effective_diffusion_weight = diffusion_weight
            effective_cluster_weight = cluster_weight (if centers initialized)
    """
    model.train()
    device = device or torch.device("cuda" if torch.cuda.is_available() else "cpu")

    epoch_losses = {"mask": 0.0, "recon": 0.0, "diffusion": 0.0, "cluster": 0.0, "total": 0.0}
    n_batches = 0

    iterator = tqdm(dataloader, desc=f"Epoch {epoch}") if progress_bar else dataloader

    for batch in iterator:
        x = batch["values"].to(device)
        support = batch.get("support", (x > 0).float()).to(device)

        # ── Forward pass ─────────────────────────────────────────────────────
        result = model(x, support=support, return_recon=True, sample_diffusion=False)
        losses = result["losses"]

        # ── Staged loss weighting ────────────────────────────────────────────
        if epoch < warmup_epochs:
            # Phase 1: Only mask + reconstruction
            effective_cluster_weight = 0.0
            effective_diffusion_weight = 0.0
        elif epoch < warmup_epochs + 50:
            # Phase 2: Gentle diffusion introduction
            effective_cluster_weight = 0.0
            effective_diffusion_weight = min(diffusion_loss_weight, 0.05)
        else:
            # Phase 3: Full training
            effective_cluster_weight = cluster_loss_weight
            effective_diffusion_weight = diffusion_loss_weight

        total_loss = (
            mask_loss_weight * losses["mask"]
            + recon_loss_weight * losses["recon"]
            + effective_diffusion_weight * losses["diffusion"]
            + effective_cluster_weight * losses["cluster"]
        )

        # ── Backward ────────────────────────────────────────────────────────
        optimizer.zero_grad()
        total_loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)
        optimizer.step()

        # Accumulate
        epoch_losses["mask"] += losses["mask"].item()
        epoch_losses["recon"] += losses["recon"].item()
        epoch_losses["diffusion"] += losses["diffusion"].item()
        epoch_losses["cluster"] += losses["cluster"].item()
        epoch_losses["total"] += total_loss.item()
        n_batches += 1

        if progress_bar:
            iterator.set_postfix(
                mask=f"{losses['mask'].item():.4f}",
                recon=f"{losses['recon'].item():.4f}",
                diff=f"{losses['diffusion'].item():.4f}",
                total=f"{total_loss.item():.4f}",
            )

    # Average over batches
    for k in epoch_losses:
        epoch_losses[k] /= max(n_batches, 1)

    return epoch_losses


# ── Embedding extraction ────────────────────────────────────────────────────────


def extract_embeddings(
    model: nn.Module,
    dataloader: torch.utils.data.DataLoader,
    device: torch.device,
    has_labels: bool = True,
) -> Tuple[np.ndarray, Optional[np.ndarray]]:
    """Extract both direct and diffusion embeddings from the model.

    Args:
        model: Trained ScSpade model.
        dataloader: DataLoader providing batches.
        device: Compute device.
        has_labels: Whether labels are available in the dataloader.

    Returns:
        (embeddings_direct, embeddings_diffusion, labels)
        embeddings_direct: (n_cells, latent_dim), raw encoder output
        embeddings_diffusion: (n_cells, latent_dim), full diffusion-sampled
        labels: (n_cells,) or None
    """
    model.eval()
    embeddings_direct_list = []
    embeddings_diffusion_list = []
    labels_list = []

    with torch.no_grad():
        for batch in dataloader:
            x = batch["values"].to(device)

            # Direct embedding: raw encoder output (default for clustering)
            z_direct = model.get_embedding(x, use_diffusion=False)
            embeddings_direct_list.append(z_direct.cpu().numpy())

            # Diffusion embedding: full reverse diffusion
            z_diff = model.get_embedding(x, use_diffusion=True)
            embeddings_diffusion_list.append(z_diff.cpu().numpy())

            if has_labels and "labels" in batch:
                labels_list.append(batch["labels"].cpu().numpy())

    embeddings_direct = np.concatenate(embeddings_direct_list, axis=0)
    embeddings_diffusion = np.concatenate(embeddings_diffusion_list, axis=0)

    if labels_list:
        labels = np.concatenate(labels_list, axis=0)
    else:
        labels = None

    return embeddings_direct, embeddings_diffusion, labels


def extract_pretrained_embeddings(
    model: nn.Module,
    dataloader: torch.utils.data.DataLoader,
    device: torch.device,
) -> np.ndarray:
    """Extract raw encoder embeddings (no diffusion) for cluster center init.

    Uses the cleanest possible embedding: direct encoder output.
    """
    model.eval()
    emb_list = []
    with torch.no_grad():
        for batch in dataloader:
            x = batch["values"].to(device)
            z = model.get_direct_embedding(x)
            emb_list.append(z.cpu().numpy())
    return np.concatenate(emb_list, axis=0)


# ── Main training loop ─────────────────────────────────────────────────────────


def run_scspade_training(
    bundle: SCDatasetBundle,
    n_clusters: int,
    *,
    latent_dim: int = 32,
    hidden_dim: int = 256,
    diffusion_hidden_dim: int = 256,
    diffusion_steps: int = 100,
    mask_hidden_dims: list = None,
    dropout: float = 0.1,
    mask_dropout: float = 0.1,
    epochs: int = 150,
    batch_size: int = 256,
    lr: float = 1e-3,
    weight_decay: float = 1e-4,
    warmup_epochs: int = 30,
    mask_loss_weight: float = 0.2,
    recon_loss_weight: float = 0.8,
    diffusion_loss_weight: float = 0.1,
    cluster_loss_weight: float = 0.0,
    eval_interval: int = 10,
    save_dir: str = None,
    device: torch.device = None,
    seed: int = 42,
    progress_bar: bool = True,
    eval_fn: Optional[Callable] = None,
    eval_fn_args: dict = None,
    label_encoder_fn: Optional[Callable] = None,
) -> dict:
    """Run the full ScSpade training pipeline.

    Args:
        bundle: Preprocessed dataset bundle.
        n_clusters: Number of cell-type clusters.
        latent_dim: Latent embedding dimension.
        hidden_dim: Encoder/decoder hidden dimension.
        diffusion_hidden_dim: Denoiser hidden dimension.
        diffusion_steps: Number of diffusion steps.
        mask_hidden_dims: Hidden dims for mask network.
        dropout: Dropout rate for encoder/decoder.
        mask_dropout: Dropout rate for mask network.
        epochs: Total training epochs.
        batch_size: Batch size.
        lr: Learning rate.
        weight_decay: L2 regularization.
        warmup_epochs: Phase 1 warmup (mask + recon only).
        mask_loss_weight: Weight for mask BCE loss.
        recon_loss_weight: Weight for reconstruction loss.
        diffusion_loss_weight: Weight for diffusion loss (gradually introduced).
        cluster_loss_weight: Weight for DEC clustering loss (activated late).
        eval_interval: Evaluate every N epochs.
        save_dir: Directory to save checkpoints and results.
        device: Compute device (auto-detected if None).
        seed: Random seed.
        progress_bar: Show tqdm progress bar.
        eval_fn: Optional evaluation function called every eval_interval epochs.
                 Signature: fn(embeddings_direct, embeddings_diff, labels) -> dict.
        eval_fn_args: Extra args passed to eval_fn.
        label_encoder_fn: Function to map predicted labels (for cluster_loss).
                          If None, uses KMeans directly on embeddings.

    Returns:
        dict with keys: model, embeddings_direct, embeddings_diffusion,
                        train_losses, eval_history
    """
    # ── Setup ─────────────────────────────────────────────────────────────────
    if device is None:
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    if save_dir:
        os.makedirs(save_dir, exist_ok=True)

    torch.manual_seed(seed)
    np.random.seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)

    n_genes = bundle.values.shape[1]
    print(f"\n[ScSpade Training] {bundle.values.shape[0]} cells × {n_genes} genes, "
          f"{n_clusters} clusters, latent={latent_dim}, device={device}")
    print(f"[ScSpade] Loss weights: mask={mask_loss_weight}, recon={recon_loss_weight}, "
          f"diff={diffusion_loss_weight}, cluster={cluster_loss_weight}")
    print(f"[ScSpade] Schedule: warmup={warmup_epochs}, "
          f"diffusion_start={warmup_epochs}, "
          f"cluster_start={warmup_epochs + 50}")

    # ── DataLoaders ──────────────────────────────────────────────────────────
    train_loader = build_dataloader(
        values=bundle.values,
        support=bundle.support,
        labels=bundle.labels,
        batch_size=batch_size,
        shuffle=True,
    )

    full_loader = build_dataloader(
        values=bundle.values,
        support=bundle.support,
        labels=bundle.labels,
        batch_size=batch_size * 2,
        shuffle=False,
    )

    # ── Model ────────────────────────────────────────────────────────────────
    model = ScSpade(
        n_genes=n_genes,
        n_clusters=n_clusters,
        latent_dim=latent_dim,
        hidden_dim=hidden_dim,
        diffusion_hidden_dim=diffusion_hidden_dim,
        diffusion_steps=diffusion_steps,
        dropout=dropout,
        mask_dropout=mask_dropout,
        mask_hidden_dims=mask_hidden_dims,
    ).to(device)

    optimizer = torch.optim.AdamW(
        model.parameters(), lr=lr, weight_decay=weight_decay
    )

    # Cosine annealing with warmup restarts
    scheduler = torch.optim.lr_scheduler.CosineAnnealingWarmRestarts(
        optimizer, T_0=warmup_epochs + 20, T_mult=2, eta_min=1e-5
    )

    # ── Training ─────────────────────────────────────────────────────────────
    train_losses = []
    eval_history = []
    best_eval_score = -1.0

    t_start = time.time()

    for epoch in range(epochs):
        epoch_losses = train_scspade_epoch(
            model=model,
            dataloader=train_loader,
            optimizer=optimizer,
            epoch=epoch,
            warmup_epochs=warmup_epochs,
            mask_loss_weight=mask_loss_weight,
            recon_loss_weight=recon_loss_weight,
            diffusion_loss_weight=diffusion_loss_weight,
            cluster_loss_weight=cluster_loss_weight,
            device=device,
            progress_bar=progress_bar,
        )
        train_losses.append(epoch_losses)
        scheduler.step()

        # ── Cluster center initialization (after warmup) ────────────────────
        if epoch == warmup_epochs - 1:
            print(f"\n[ScSpade] Epoch {epoch}: Initializing cluster centers with KMeans...")
            pretrained_emb = extract_pretrained_embeddings(model, full_loader, device)
            model.init_cluster_centers(
                torch.from_numpy(pretrained_emb).float().to(device)
            )
            print(f"[ScSpade] Cluster centers initialized. "
                  f"DEC clustering loss now ACTIVE at epoch {epoch + 1}.")

        # ── Evaluation ─────────────────────────────────────────────────────
        if (epoch + 1) % eval_interval == 0 or epoch == epochs - 1:
            print(f"\n[ScSpade] Epoch {epoch}: Extracting embeddings...")
            emb_direct, emb_diff, labels = extract_embeddings(
                model, full_loader, device, has_labels=True
            )

            print(f"[ScSpade] Direct embedding: mean={emb_direct.mean():.4f}, "
                  f"std={emb_direct.std():.4f}")
            print(f"[ScSpade] Diffusion embedding: mean={emb_diff.mean():.4f}, "
                  f"std={emb_diff.std():.4f}")

            eval_result = {
                "epoch": epoch,
                "emb_direct_stats": {
                    "mean": float(emb_direct.mean()),
                    "std": float(emb_direct.std()),
                    "min": float(emb_direct.min()),
                    "max": float(emb_direct.max()),
                },
                "emb_diff_stats": {
                    "mean": float(emb_diff.mean()),
                    "std": float(emb_diff.std()),
                    "min": float(emb_diff.min()),
                    "max": float(emb_diff.max()),
                },
            }

            if eval_fn is not None:
                fn_args = eval_fn_args or {}
                metrics_direct = eval_fn(emb_direct, bundle.labels, **fn_args)
                metrics_diff = eval_fn(emb_diff, bundle.labels, **fn_args)
                eval_result["metrics_direct"] = metrics_direct
                eval_result["metrics_diff"] = metrics_diff

                # Track best
                score = metrics_direct.get("nmi", metrics_direct.get("ari", -1))
                if score > best_eval_score:
                    best_eval_score = score
                    if save_dir:
                        torch.save(
                            model.state_dict(),
                            os.path.join(save_dir, "best_model.pt"),
                        )
                        np.save(os.path.join(save_dir, "best_embeddings_direct.npy"), emb_direct)
                        np.save(os.path.join(save_dir, "best_embeddings_diffusion.npy"), emb_diff)
                    print(f"[ScSpade] New best eval NMI={score:.4f} at epoch {epoch}")

            eval_history.append(eval_result)

        # ── Checkpoint ─────────────────────────────────────────────────────
        checkpoint_every = max(1, epochs // 3)
        if save_dir and (epoch + 1) % checkpoint_every == 0:
            torch.save(
                model.state_dict(),
                os.path.join(save_dir, f"checkpoint_epoch_{epoch}.pt"),
            )

    t_end = time.time()
    print(f"\n[ScSpade] Training complete in {t_end - t_start:.1f}s")

    # ── Final extraction ────────────────────────────────────────────────────
    print("[ScSpade] Final embedding extraction...")
    emb_direct, emb_diff, labels = extract_embeddings(
        model, full_loader, device, has_labels=True
    )

    # Final eval
    final_result = {
        "epoch": epochs - 1,
        "train_losses": train_losses,
        "eval_history": eval_history,
        "best_eval_nmi": best_eval_score,
    }

    if eval_fn is not None:
        fn_args = eval_fn_args or {}
        final_result["final_metrics_direct"] = eval_fn(emb_direct, bundle.labels, **fn_args)
        final_result["final_metrics_diff"] = eval_fn(emb_diff, bundle.labels, **fn_args)

    return {
        "model": model,
        "embeddings_direct": emb_direct,
        "embeddings_diffusion": emb_diff,
        "labels": labels,
        "train_losses": train_losses,
        "eval_history": eval_history,
        "final_result": final_result,
        "device": device,
    }
