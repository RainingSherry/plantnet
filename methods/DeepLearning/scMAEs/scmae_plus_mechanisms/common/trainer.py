from __future__ import annotations

import argparse
import os
import subprocess
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable

import numpy as np
import torch
from torch.utils.data import DataLoader

ROOT = Path(__file__).resolve().parents[5]
MECH_DIR = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")
os.environ.setdefault("NUMBA_CACHE_DIR", "/tmp/numba-cache")
os.makedirs(os.environ["MPLCONFIGDIR"], exist_ok=True)
os.makedirs(os.environ["NUMBA_CACHE_DIR"], exist_ok=True)

from methods.DeepLearning import scMAE_family as family
from methods.shared_utils import ensure_dir, sanitize_anndata_for_write, save_json

from . import diagnostics, masks, neighbormix, prototype


@dataclass
class VariantConfig:
    name: str
    method_name: str
    description: str
    defaults: dict = field(default_factory=dict)


def linear_warmup(epoch: int, start_epoch: int, max_weight: float, warmup_epochs: int) -> float:
    if epoch < int(start_epoch):
        return 0.0
    progress = (epoch - int(start_epoch) + 1) / max(1, int(warmup_epochs))
    return float(max_weight) * min(1.0, max(0.0, progress))


def resolve_input_path(data_path: str, dataset_name: str) -> str:
    path = Path(data_path).resolve()
    if path.suffix.lower() != ".h5":
        return str(path)
    out_dir = MECH_DIR.parent / "benchmark_data"
    out_dir.mkdir(parents=True, exist_ok=True)
    converted = out_dir / f"{dataset_name}.h5ad"
    if converted.exists():
        return str(converted)
    cmd = [
        sys.executable,
        str(ROOT / "scripts" / "prepare_dataset.py"),
        "--input_path",
        str(path),
        "--dataset_name",
        dataset_name,
        "--output_dir",
        str(out_dir),
        "--force",
    ]
    result = subprocess.run(cmd, cwd=str(ROOT), text=True, capture_output=True, timeout=3600)
    if result.returncode != 0:
        raise RuntimeError(
            "Failed to convert input .h5 to .h5ad\n"
            f"Command: {' '.join(cmd)}\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"
        )
    return str(converted)


def build_parser(config: VariantConfig) -> argparse.ArgumentParser:
    d = config.defaults
    parser = argparse.ArgumentParser(description=config.description, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--data_path", type=str, required=True)
    parser.add_argument("--save_dir", type=str, default="./results")
    parser.add_argument("--dataset_name", type=str, default=None)
    parser.add_argument("--method_name", type=str, default=config.method_name)
    parser.add_argument("--variant_name", type=str, default=config.name)
    parser.add_argument("--label_key", type=str, default="auto")
    parser.add_argument("--input_mode", type=str, default="auto", choices=["auto", "raw", "log1p"])
    parser.add_argument("--n_top_genes", type=int, default=1000)
    parser.add_argument("--target_sum", type=float, default=10000.0)
    parser.add_argument("--scale_input", type=family.str2bool, default=True)
    parser.add_argument("--n_clusters", type=int, default=0)
    parser.add_argument("--hidden_size", type=int, default=128)
    parser.add_argument("--dropout", type=float, default=0.0)
    parser.add_argument("--mask_prob", type=float, default=0.4)
    parser.add_argument("--masked_data_weight", type=float, default=0.75)
    parser.add_argument("--mask_loss_weight", type=float, default=0.7)
    parser.add_argument("--epochs", type=int, default=80)
    parser.add_argument("--batch_size", type=int, default=256)
    parser.add_argument("--lr", type=float, default=1e-3)
    parser.add_argument("--weight_decay", type=float, default=0.0)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--gpu", type=int, default=1)
    parser.add_argument("--no_cuda", action="store_true")
    parser.add_argument("--skip_eval", type=family.str2bool, default=False)
    parser.add_argument("--no_save_h5ad", action="store_true")

    parser.add_argument("--use_adaptive_mask", type=family.str2bool, default=bool(d.get("use_adaptive_mask", False)))
    parser.add_argument("--use_prototype", type=family.str2bool, default=bool(d.get("use_prototype", False)))
    parser.add_argument("--use_swav", type=family.str2bool, default=bool(d.get("use_swav", False)))
    parser.add_argument("--use_neighbormix", type=family.str2bool, default=bool(d.get("use_neighbormix", False)))
    parser.add_argument("--use_count_aux", type=family.str2bool, default=False)
    parser.add_argument("--use_fuzzy_boundary", type=family.str2bool, default=False)

    parser.add_argument("--mask_mode", type=str, default=d.get("mask_mode", "variance_adaptive"))
    parser.add_argument("--mask_p_min", type=float, default=0.05)
    parser.add_argument("--mask_p_max", type=float, default=0.75)
    parser.add_argument("--mask_strength", type=float, default=0.75)
    parser.add_argument("--module_size", type=int, default=32)

    parser.add_argument("--prototype_start_epoch", type=int, default=int(d.get("prototype_start_epoch", 20)))
    parser.add_argument("--prototype_warmup_epochs", type=int, default=10)
    parser.add_argument("--prototype_weight", type=float, default=float(d.get("prototype_weight", 0.1)))
    parser.add_argument("--prototype_confidence_threshold", type=float, default=0.6)
    parser.add_argument("--prototype_alpha", type=float, default=1.0)

    parser.add_argument("--swav_start_epoch", type=int, default=int(d.get("swav_start_epoch", 10)))
    parser.add_argument("--swav_warmup_epochs", type=int, default=10)
    parser.add_argument("--swav_weight", type=float, default=float(d.get("swav_weight", 0.1)))
    parser.add_argument("--swav_temperature", type=float, default=0.1)
    parser.add_argument("--swav_pred_temperature", type=float, default=0.1)
    parser.add_argument("--sinkhorn_iterations", type=int, default=3)

    parser.add_argument("--neighbor_start_epoch", type=int, default=int(d.get("neighbor_start_epoch", 20)))
    parser.add_argument("--neighbor_update_interval", type=int, default=10)
    parser.add_argument("--neighbor_k", type=int, default=15)
    parser.add_argument("--neighbor_metric", type=str, default="cosine")
    parser.add_argument(
        "--neighbor_mix_mode",
        type=str,
        default=d.get("neighbor_mix_mode", "first"),
        choices=["first", "soft_first", "mean", "weighted_mean"],
    )
    parser.add_argument("--neighbor_min_similarity", type=float, default=float(d.get("neighbor_min_similarity", -1.0)))
    parser.add_argument("--neighbor_min_shared_score", type=float, default=float(d.get("neighbor_min_shared_score", 0.0)))
    parser.add_argument("--neighbor_score_threshold", type=float, default=float(d.get("neighbor_score_threshold", 0.0)))
    parser.add_argument(
        "--neighbor_score_similarity_weight",
        type=float,
        default=float(d.get("neighbor_score_similarity_weight", 0.7)),
    )
    parser.add_argument(
        "--neighbor_pseudo_filter",
        type=str,
        default=d.get("neighbor_pseudo_filter", "none"),
        choices=["none", "same_cluster", "same_confident_cluster"],
    )
    parser.add_argument(
        "--neighbor_pseudo_confidence_quantile",
        type=float,
        default=float(d.get("neighbor_pseudo_confidence_quantile", 0.0)),
    )
    parser.add_argument(
        "--neighbor_graph_embedding",
        type=str,
        default=d.get("neighbor_graph_embedding", "current"),
        choices=["current", "ema"],
    )
    parser.add_argument(
        "--neighbor_embedding_ema_decay",
        type=float,
        default=float(d.get("neighbor_embedding_ema_decay", 0.8)),
    )
    parser.add_argument("--neighbor_consensus_window", type=int, default=int(d.get("neighbor_consensus_window", 1)))
    parser.add_argument("--neighbor_consensus_min_hits", type=int, default=int(d.get("neighbor_consensus_min_hits", 1)))
    parser.add_argument(
        "--neighbor_adaptive_consensus",
        type=family.str2bool,
        default=bool(d.get("neighbor_adaptive_consensus", False)),
    )
    parser.add_argument("--neighbor_adaptive_loose_hits", type=int, default=int(d.get("neighbor_adaptive_loose_hits", 2)))
    parser.add_argument("--neighbor_adaptive_strict_hits", type=int, default=int(d.get("neighbor_adaptive_strict_hits", 3)))
    parser.add_argument(
        "--neighbor_adaptive_score_threshold",
        type=float,
        default=float(d.get("neighbor_adaptive_score_threshold", 0.84)),
    )
    parser.add_argument(
        "--neighbor_boundary_protect",
        type=family.str2bool,
        default=bool(d.get("neighbor_boundary_protect", False)),
    )
    parser.add_argument(
        "--neighbor_boundary_confidence_quantile",
        type=float,
        default=float(d.get("neighbor_boundary_confidence_quantile", 0.20)),
    )
    parser.add_argument(
        "--neighbor_boundary_rare_quantile",
        type=float,
        default=float(d.get("neighbor_boundary_rare_quantile", 0.25)),
    )
    parser.add_argument(
        "--neighbor_boundary_score_threshold",
        type=float,
        default=float(d.get("neighbor_boundary_score_threshold", 0.84)),
    )
    parser.add_argument("--neighbor_soft_power", type=float, default=float(d.get("neighbor_soft_power", 1.0)))
    parser.add_argument("--mix_alpha", type=float, default=float(d.get("mix_alpha", 0.9)))
    parser.add_argument("--mix_weight", type=float, default=float(d.get("mix_weight", 0.3)))
    parser.add_argument("--mix_warmup_epochs", type=int, default=10)
    return parser


def _full_embedding(model, loader: DataLoader, device: torch.device) -> tuple[np.ndarray, np.ndarray]:
    return family.extract_embedding(model, loader, device)


def _maybe_initialize_prototypes(
    model,
    args,
    epoch: int,
    initialized: bool,
    loader: DataLoader,
    device: torch.device,
    n_clusters: int,
) -> bool:
    if initialized or not (args.use_prototype or args.use_swav) or not model.has_prototypes:
        return initialized
    start_epoch = min(
        args.prototype_start_epoch if args.use_prototype else 10**9,
        args.swav_start_epoch if args.use_swav else 10**9,
    )
    if epoch < start_epoch:
        return initialized
    embedding, _ = _full_embedding(model, loader, device)
    centers = prototype.kmeans_centers(embedding, n_clusters=n_clusters, seed=args.seed)
    model.set_prototypes(centers)
    model.normalize_prototypes()
    return True


def _maybe_build_neighbors(
    model,
    args,
    epoch: int,
    loader: DataLoader,
    device: torch.device,
    current_state: neighbormix.NeighborState | None,
    n_clusters: int,
) -> tuple[neighbormix.NeighborState | None, np.ndarray | None]:
    if not args.use_neighbormix or epoch < args.neighbor_start_epoch:
        return current_state, None
    if current_state is not None and (epoch - args.neighbor_start_epoch) % max(1, args.neighbor_update_interval) != 0:
        return current_state, None
    embedding, _ = _full_embedding(model, loader, device)
    return _build_neighbor_state_from_embedding(args, embedding, n_clusters), embedding


def _build_neighbor_state_from_embedding(
    args,
    graph_embedding: np.ndarray,
    n_clusters: int,
) -> neighbormix.NeighborState:
    state = neighbormix.build_neighbor_state(
        graph_embedding,
        k=args.neighbor_k,
        metric=args.neighbor_metric,
        min_similarity=args.neighbor_min_similarity,
        min_shared_score=args.neighbor_min_shared_score,
        score_threshold=args.neighbor_score_threshold,
        similarity_weight=args.neighbor_score_similarity_weight,
        pseudo_filter=args.neighbor_pseudo_filter,
        pseudo_n_clusters=n_clusters,
        pseudo_seed=args.seed,
        pseudo_confidence_quantile=args.neighbor_pseudo_confidence_quantile,
    )
    state.stats["neighbor_graph_embedding"] = str(args.neighbor_graph_embedding)
    state.stats["neighbor_embedding_ema_decay"] = float(args.neighbor_embedding_ema_decay)
    return state


def _loss_components_template() -> dict[str, float]:
    return {
        "loss": 0.0,
        "scmae_loss": 0.0,
        "prototype_loss": 0.0,
        "swav_loss": 0.0,
        "mix_loss": 0.0,
        "prototype_weight": 0.0,
        "swav_weight": 0.0,
        "mix_weight": 0.0,
        "effective_mask_rate": 0.0,
        "mix_used_fraction": 0.0,
        "prototype_used_fraction": 0.0,
    }


def run_training(config: VariantConfig, build_model: Callable) -> None:
    args = build_parser(config).parse_args()
    if args.use_count_aux:
        raise SystemExit("count-aware auxiliary is planned but not enabled in first-stage variants.")
    if args.use_fuzzy_boundary:
        raise SystemExit("fuzzy boundary is planned but not enabled in first-stage variants.")

    family.set_seed(args.seed)
    save_dir = Path(ensure_dir(args.save_dir))
    save_json(vars(args), str(save_dir / "args.json"))
    device = family.get_device(args.gpu, args.no_cuda)
    dataset_name = args.dataset_name or Path(args.data_path).stem

    data_path = resolve_input_path(args.data_path, dataset_name)
    bundle = family.load_scmae_dataset(
        file_path=data_path,
        input_mode=args.input_mode,
        n_top_genes=args.n_top_genes,
        target_sum=args.target_sum,
        scale_input=args.scale_input,
        label_key=args.label_key,
        seed=args.seed,
    )
    save_json(bundle.profile, str(save_dir / "dataset_profile.json"))
    save_json(bundle.preprocess_config, str(save_dir / "preprocess_config.json"))
    with open(save_dir / "selected_genes.txt", "w", encoding="utf-8") as handle:
        for gene in bundle.gene_names:
            handle.write(f"{gene}\n")

    data_np = bundle.data
    labels = bundle.labels
    n_clusters = args.n_clusters if args.n_clusters > 0 else int(len(np.unique(labels)))

    mask_state = None
    if args.use_adaptive_mask:
        mask_state = masks.build_mask_state(
            data=data_np,
            mode=args.mask_mode,
            base_rate=args.mask_prob,
            p_min=args.mask_p_min,
            p_max=args.mask_p_max,
            strength=args.mask_strength,
            module_size=args.module_size,
        )

    variant_payload = {
        "variant": config.name,
        "method_name": args.method_name,
        "description": config.description,
        "defaults": config.defaults,
        "active_mechanisms": {
            "adaptive_mask": bool(args.use_adaptive_mask),
            "prototype": bool(args.use_prototype),
            "swav": bool(args.use_swav),
            "neighbormix": bool(args.use_neighbormix),
            "count_aux": bool(args.use_count_aux),
            "fuzzy_boundary": bool(args.use_fuzzy_boundary),
        },
        "mask_state": None if mask_state is None else mask_state.stats,
    }
    save_json(variant_payload, str(save_dir / "variant_config.json"))

    dataset = family.IndexedExpressionDataset(data_np, labels)
    generator = torch.Generator()
    generator.manual_seed(args.seed)
    train_loader = DataLoader(
        dataset,
        batch_size=args.batch_size,
        shuffle=True,
        drop_last=True,
        generator=generator,
    )
    eval_loader = DataLoader(
        dataset,
        batch_size=max(args.batch_size * 4, 512),
        shuffle=False,
        drop_last=False,
    )
    full_data_cpu = torch.as_tensor(data_np, dtype=torch.float32)

    model = build_model(num_genes=data_np.shape[1], args=args, n_clusters=n_clusters).to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=args.lr, weight_decay=args.weight_decay)

    history = {"epoch": []}
    for key in _loss_components_template():
        history[key] = []
    history["configured_mask_rate"] = float(args.mask_prob)
    history["mask_state"] = None if mask_state is None else mask_state.stats

    proto_initialized = False
    neighbor_state = None
    neighbor_embedding_ema = None
    neighbor_state_history = []
    last_components = _loss_components_template()
    print(
        f"{config.name}: cells={data_np.shape[0]} genes={data_np.shape[1]} "
        f"clusters={n_clusters} device={device}",
        flush=True,
    )

    for epoch in range(1, max(1, args.epochs) + 1):
        proto_initialized = _maybe_initialize_prototypes(
            model, args, epoch, proto_initialized, eval_loader, device, n_clusters
        )
        maybe_state, current_neighbor_embedding = _maybe_build_neighbors(
            model, args, epoch, eval_loader, device, neighbor_state, n_clusters
        )
        if current_neighbor_embedding is not None and args.neighbor_graph_embedding == "ema":
            decay = float(np.clip(args.neighbor_embedding_ema_decay, 0.0, 0.999))
            if neighbor_embedding_ema is None:
                neighbor_embedding_ema = current_neighbor_embedding.astype(np.float32, copy=True)
            else:
                neighbor_embedding_ema = (
                    decay * neighbor_embedding_ema + (1.0 - decay) * current_neighbor_embedding.astype(np.float32)
                )
            neighbor_state = _build_neighbor_state_from_embedding(args, neighbor_embedding_ema, n_clusters)
        else:
            neighbor_state = maybe_state
        if current_neighbor_embedding is not None and neighbor_state is not None:
            window = max(1, int(args.neighbor_consensus_window))
            neighbor_state_history.append(neighbor_state)
            neighbor_state_history = neighbor_state_history[-window:]
            min_hits = max(1, int(args.neighbor_consensus_min_hits))
            if args.neighbor_adaptive_consensus:
                neighbor_state = neighbormix.adaptive_consensus_neighbor_state(
                    neighbor_state_history,
                    loose_hits=args.neighbor_adaptive_loose_hits,
                    strict_hits=args.neighbor_adaptive_strict_hits,
                    score_threshold=args.neighbor_adaptive_score_threshold,
                )
            elif window > 1 or min_hits > 1:
                neighbor_state = neighbormix.consensus_neighbor_state(neighbor_state_history, min_hits=min_hits)
            if args.neighbor_boundary_protect:
                graph_embedding_for_boundary = (
                    neighbor_embedding_ema
                    if args.neighbor_graph_embedding == "ema" and neighbor_embedding_ema is not None
                    else current_neighbor_embedding
                )
                neighbor_state = neighbormix.boundary_protected_neighbor_state(
                    neighbor_state,
                    graph_embedding_for_boundary,
                    n_clusters=n_clusters,
                    seed=args.seed,
                    confidence_quantile=args.neighbor_boundary_confidence_quantile,
                    rare_quantile=args.neighbor_boundary_rare_quantile,
                    score_threshold=args.neighbor_boundary_score_threshold,
                )
        model.train()
        totals = _loss_components_template()
        n_batches = 0

        for idx_cpu, x_cpu, _ in train_loader:
            x = x_cpu.to(device)
            optimizer.zero_grad()
            proto_weight = linear_warmup(
                epoch, args.prototype_start_epoch, args.prototype_weight, args.prototype_warmup_epochs
            )
            swav_weight = linear_warmup(epoch, args.swav_start_epoch, args.swav_weight, args.swav_warmup_epochs)
            mix_weight = linear_warmup(epoch, args.neighbor_start_epoch, args.mix_weight, args.mix_warmup_epochs)

            if args.use_swav:
                x_a, mask_a = masks.apply_replacement_noise(x, mask_state, args.mask_prob)
                x_b, mask_b = masks.apply_replacement_noise(x, mask_state, args.mask_prob)
                z_a, loss_a = model.scmae_loss(x_a, x, mask_a)
                z_b, loss_b = model.scmae_loss(x_b, x, mask_b)
                scmae_loss = 0.5 * (loss_a + loss_b)
                effective_mask = 0.5 * (mask_a.mean() + mask_b.mean())
                logits_a = model.prototype_logits(z_a)
                logits_b = model.prototype_logits(z_b)
                assign_a = prototype.sinkhorn_assignments(logits_a.detach(), args.swav_temperature, args.sinkhorn_iterations)
                assign_b = prototype.sinkhorn_assignments(logits_b.detach(), args.swav_temperature, args.sinkhorn_iterations)
                swav_loss = prototype.swapped_assignment_loss(
                    logits_a, logits_b, assign_a, assign_b, args.swav_pred_temperature
                )
                proto_loss = torch.zeros((), device=device)
                total_loss = scmae_loss + swav_weight * swav_loss
                proto_used_fraction = 0.0
            else:
                x_corrupted, mask = masks.apply_replacement_noise(x, mask_state, args.mask_prob)
                z, scmae_loss = model.scmae_loss(x_corrupted, x, mask)
                effective_mask = mask.mean()
                swav_loss = torch.zeros((), device=device)
                proto_loss = torch.zeros((), device=device)
                proto_used_fraction = 0.0
                if args.use_prototype and proto_initialized and proto_weight > 0.0:
                    proto_loss, proto_stats = prototype.dec_kl_loss(
                        z,
                        model.prototype_centers(),
                        alpha=args.prototype_alpha,
                        confidence_threshold=args.prototype_confidence_threshold,
                    )
                    proto_used_fraction = proto_stats["prototype_used_fraction"]
                total_loss = scmae_loss + proto_weight * proto_loss

            mix_loss = torch.zeros((), device=device)
            mix_used_fraction = 0.0
            if args.use_neighbormix and neighbor_state is not None and mix_weight > 0.0:
                x_mix, mix_used_fraction = neighbormix.mix_batch(
                    idx_cpu,
                    x,
                    full_data_cpu,
                    neighbor_state,
                    alpha=args.mix_alpha,
                    mode=args.neighbor_mix_mode,
                    soft_power=args.neighbor_soft_power,
                )
                x_mix_corrupted, mix_mask = masks.apply_replacement_noise(x_mix, mask_state, args.mask_prob)
                _, mix_loss = model.scmae_loss(x_mix_corrupted, x, mix_mask)
                total_loss = total_loss + mix_weight * mix_loss

            total_loss.backward()
            optimizer.step()
            if args.use_swav and model.has_prototypes:
                model.normalize_prototypes()

            totals["loss"] += float(total_loss.detach().cpu())
            totals["scmae_loss"] += float(scmae_loss.detach().cpu())
            totals["prototype_loss"] += float(proto_loss.detach().cpu())
            totals["swav_loss"] += float(swav_loss.detach().cpu())
            totals["mix_loss"] += float(mix_loss.detach().cpu())
            totals["prototype_weight"] += float(proto_weight if args.use_prototype else 0.0)
            totals["swav_weight"] += float(swav_weight if args.use_swav else 0.0)
            totals["mix_weight"] += float(mix_weight if args.use_neighbormix else 0.0)
            totals["effective_mask_rate"] += float(effective_mask.detach().cpu())
            totals["mix_used_fraction"] += float(mix_used_fraction)
            totals["prototype_used_fraction"] += float(proto_used_fraction)
            n_batches += 1

        history["epoch"].append(int(epoch))
        last_components = {}
        for key, value in totals.items():
            avg = value / max(1, n_batches)
            history[key].append(avg)
            last_components[key] = avg
        if epoch == 1 or epoch == args.epochs or epoch % 10 == 0:
            print(
                f"{config.name} epoch={epoch:03d}/{args.epochs} "
                f"loss={last_components['loss']:.4f} scmae={last_components['scmae_loss']:.4f} "
                f"proto={last_components['prototype_loss']:.4f} swav={last_components['swav_loss']:.4f} "
                f"mix={last_components['mix_loss']:.4f} mask={last_components['effective_mask_rate']:.4f}",
                flush=True,
            )

    embedding, labels_out = _full_embedding(model, eval_loader, device)
    np.save(save_dir / "embedding_final.npy", embedding.astype(np.float32))
    np.save(save_dir / "labels.npy", labels_out.astype(np.int64))
    np.save(save_dir / "gene_names.npy", bundle.gene_names.astype(str))
    family.save_embedding_h5(save_dir / "embedding.h5", embedding, labels_out)
    save_json(history, str(save_dir / "training_history.json"))
    torch.save(
        {
            "model": model.state_dict(),
            "optimizer": optimizer.state_dict(),
            "args": vars(args),
            "variant_config": variant_payload,
            "gene_names": bundle.gene_names.astype(str),
        },
        save_dir / "model_checkpoint.pth",
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
                "backbone": "original_scmae",
                "use_adaptive_mask": bool(args.use_adaptive_mask),
                "use_prototype": bool(args.use_prototype),
                "use_swav": bool(args.use_swav),
                "use_neighbormix": bool(args.use_neighbormix),
                "mask_mode": args.mask_mode if args.use_adaptive_mask else "random",
            },
        )
        pred = result["preds"]["kmeans_known_k"]
        save_json(result["fixed"], str(save_dir / "metrics.json"))

    prototype_conf = None
    if model.has_prototypes:
        prototype_conf = prototype.confidence_from_embedding(embedding, model.prototype_centers(), device)
    diag = diagnostics.build_diagnostics(
        embedding=embedding,
        pred_labels=pred,
        configured_mask_rate=args.mask_prob,
        effective_mask_rate=last_components.get("effective_mask_rate", 0.0),
        prototype_confidence_mean=prototype_conf,
        neighbor_stats=None if neighbor_state is None else neighbor_state.stats,
        loss_components=last_components,
    )
    save_json(diag, str(save_dir / "diagnostics.json"))

    if not args.no_save_h5ad:
        bundle.adata.obsm["X_scmae_plus"] = embedding
        bundle.adata.uns["scmae_plus"] = {
            "method": args.method_name,
            "variant": args.variant_name,
            "backbone": "original_scmae",
        }
        sanitize_anndata_for_write(bundle.adata)
        bundle.adata.write_h5ad(save_dir / "adata_scmae_plus.h5ad", compression="gzip")

    summary = {
        "dataset": dataset_name,
        "method": args.method_name,
        "variant": args.variant_name,
        "seed": int(args.seed),
        "n_cells": int(data_np.shape[0]),
        "n_genes": int(data_np.shape[1]),
        "n_clusters": int(n_clusters),
        "fixed_metrics": result["fixed"] if result is not None else {},
        "diagnostics": diag,
        "note": "Mechanism-search variant that preserves original scMAE backbone and shared KMeans-known-k protocol.",
    }
    save_json(summary, str(save_dir / "summary.json"))
    print(f"Training completed. Results saved to: {save_dir}", flush=True)


def main(config: VariantConfig, build_model: Callable) -> None:
    run_training(config, build_model)
