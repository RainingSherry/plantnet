from __future__ import annotations

import argparse
import copy
import json
import os
import subprocess
import sys
import time
from pathlib import Path
from typing import Dict, Optional, Tuple

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")
os.environ.setdefault("NUMBA_CACHE_DIR", "/tmp/numba-cache")
os.makedirs(os.environ["MPLCONFIGDIR"], exist_ok=True)
os.makedirs(os.environ["NUMBA_CACHE_DIR"], exist_ok=True)

import numpy as np
import pandas as pd
import torch
import torch.nn.functional as F
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import normalize
from torch.utils.data import DataLoader

CURRENT_DIR = Path(__file__).resolve().parent
SCMAES_DIR = CURRENT_DIR.parent
ROOT = SCMAES_DIR.parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from methods.DeepLearning import scMAE_family as family
from methods.shared_utils import ensure_dir, save_json

from .model import FlexibleScMAE
from .variant_defs import VARIANTS, VariantConfig, get_variant


def parse_args(variant_key: str) -> argparse.Namespace:
    cfg = get_variant(variant_key)
    parser = argparse.ArgumentParser(description=f"scMAEs variant runner: {variant_key}")
    parser.add_argument("--data_path", required=True)
    parser.add_argument("--save_dir", required=True)
    parser.add_argument("--dataset_name", default=None)
    parser.add_argument("--method_name", default=cfg.method_name)
    parser.add_argument("--variant_name", default=variant_key)
    parser.add_argument("--label_key", default="auto")
    parser.add_argument("--input_mode", default="auto", choices=["auto", "raw", "log1p"])
    parser.add_argument("--n_top_genes", type=int, default=1000)
    parser.add_argument("--target_sum", type=float, default=10000.0)
    parser.add_argument("--scale_input", type=family.str2bool, default=True)
    parser.add_argument("--n_clusters", type=int, default=0)
    parser.add_argument("--hidden_size", type=int, default=128)
    parser.add_argument("--dropout", type=float, default=None)
    parser.add_argument("--epochs", type=int, default=80)
    parser.add_argument("--batch_size", type=int, default=256)
    parser.add_argument("--lr", type=float, default=1e-3)
    parser.add_argument("--weight_decay", type=float, default=1e-5)
    parser.add_argument("--mask_ratio", type=float, default=None)
    parser.add_argument("--masked_data_weight", type=float, default=None)
    parser.add_argument("--mask_loss_weight", type=float, default=None)
    parser.add_argument("--neighbor_k", type=int, default=None)
    parser.add_argument("--knn_pca_dim", type=int, default=50)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--gpu", type=int, default=1)
    parser.add_argument("--no_cuda", action="store_true")
    parser.add_argument("--skip_eval", type=family.str2bool, default=False)
    parser.add_argument("--no_save_h5ad", action="store_true")
    parser.add_argument("--smoke", action="store_true")
    return parser.parse_args()


def build_profile(data_np: np.ndarray) -> Dict[str, np.ndarray]:
    mean = data_np.mean(axis=0).astype(np.float32)
    var = data_np.var(axis=0).astype(np.float32)
    zero_rate = (np.abs(data_np) < 1e-8).mean(axis=0).astype(np.float32)
    abs_mean = np.abs(data_np).mean(axis=0).astype(np.float32)
    dispersion = var / (np.abs(mean) + 1e-4)
    return {
        "mean": mean,
        "var": var,
        "zero_rate": zero_rate,
        "abs_mean": abs_mean,
        "dispersion": dispersion.astype(np.float32),
    }


def _scale_to_prob(values: np.ndarray, base: float) -> torch.Tensor:
    arr = np.asarray(values, dtype=np.float32)
    lo, hi = np.percentile(arr, [5, 95])
    scaled = (arr - lo) / max(float(hi - lo), 1e-6)
    scaled = np.clip(scaled, 0.0, 1.0)
    prob = base * (0.5 + scaled)
    return torch.as_tensor(np.clip(prob, 0.05, 0.85), dtype=torch.float32)


def mask_probability(config: VariantConfig, profile: Dict[str, np.ndarray], base: float, epoch: int, epochs: int) -> torch.Tensor:
    strategy = config.mask_strategy
    if strategy == "variance_adaptive":
        return _scale_to_prob(profile["var"], base)
    if strategy == "dropout_adaptive":
        score = 0.6 * profile["zero_rate"] + 0.4 * profile["dispersion"]
        return _scale_to_prob(score, base)
    if strategy == "marker_safe":
        score = profile["var"] - 0.35 * profile["abs_mean"]
        return _scale_to_prob(score, base)
    if strategy == "module_block":
        return _scale_to_prob(profile["dispersion"], base)
    if strategy == "high_mask_curriculum":
        frac = min(1.0, max(0.0, epoch / max(1, epochs)))
        return torch.full_like(torch.as_tensor(profile["var"], dtype=torch.float32), base + 0.2 * (1.0 - frac))
    if strategy in {"joao", "dropedge"}:
        blend = 0.5 * profile["var"] + 0.5 * profile["zero_rate"]
        return _scale_to_prob(blend, base)
    return torch.full_like(torch.as_tensor(profile["var"], dtype=torch.float32), base)


def apply_variant_noise(
    x: torch.Tensor,
    prob: torch.Tensor,
    strategy: str,
) -> Tuple[torch.Tensor, torch.Tensor]:
    prob = prob.to(device=x.device, dtype=x.dtype).view(1, -1).expand_as(x)
    should_swap = torch.bernoulli(prob).bool()
    if strategy == "module_block" and x.shape[1] >= 8:
        width = max(4, x.shape[1] // 40)
        starts = torch.randint(0, max(1, x.shape[1] - width), (x.shape[0],), device=x.device)
        block = torch.zeros_like(should_swap)
        for row, start in enumerate(starts.tolist()):
            block[row, start:start + width] = True
        should_swap = should_swap | block
    if x.shape[0] <= 1:
        replacement = x
    else:
        replacement = x[torch.randperm(x.shape[0], device=x.device)]
    corrupted = torch.where(should_swap, replacement, x)
    mask = (corrupted != x).to(dtype=x.dtype)
    return corrupted, mask


def build_token_bins(data_np: np.ndarray, n_bins: int) -> np.ndarray:
    qs = np.linspace(0, 100, n_bins + 1)[1:-1]
    return np.percentile(data_np.reshape(-1), qs).astype(np.float32)


def tokenize_tensor(x: torch.Tensor, bins: torch.Tensor) -> torch.Tensor:
    return torch.bucketize(x, bins.to(device=x.device, dtype=x.dtype)).long()


def target_gene_features(x: torch.Tensor) -> torch.Tensor:
    return torch.stack(
        [
            x.mean(dim=1),
            x.std(dim=1, unbiased=False),
            x.abs().mean(dim=1),
            (x.abs() < 1e-8).float().mean(dim=1),
        ],
        dim=1,
    )


def weighted_reconstruction_loss(reconstruction: torch.Tensor, target: torch.Tensor, mask: torch.Tensor, config: VariantConfig) -> torch.Tensor:
    weights = mask * config.masked_data_weight + (1.0 - mask) * (1.0 - config.masked_data_weight)
    if config.reconstruction == "huber":
        raw = F.smooth_l1_loss(reconstruction, target, reduction="none", beta=0.7)
    else:
        raw = F.mse_loss(reconstruction, target, reduction="none")
    denominator = weights.sum().clamp_min(1e-8)
    return (weights * raw).sum() / denominator


def prototype_loss(model: FlexibleScMAE, latent: torch.Tensor) -> torch.Tensor:
    q = model.soft_assignments(latent)
    with torch.no_grad():
        p = q.pow(2) / q.sum(dim=0, keepdim=True).clamp_min(1e-8)
        p = p / p.sum(dim=1, keepdim=True).clamp_min(1e-8)
    return F.kl_div((q + 1e-8).log(), p, reduction="batchmean")


def fuzzy_loss(model: FlexibleScMAE, latent: torch.Tensor) -> torch.Tensor:
    q = model.soft_assignments(latent)
    entropy = -(q * (q + 1e-8).log()).sum(dim=1).mean()
    balance = ((q.mean(dim=0) - 1.0 / q.shape[1]) ** 2).mean()
    return 0.2 * entropy + balance


def barlow_loss(z1: torch.Tensor, z2: torch.Tensor) -> torch.Tensor:
    z1 = (z1 - z1.mean(dim=0)) / z1.std(dim=0, unbiased=False).clamp_min(1e-4)
    z2 = (z2 - z2.mean(dim=0)) / z2.std(dim=0, unbiased=False).clamp_min(1e-4)
    c = torch.mm(z1.T, z2) / max(1, z1.shape[0])
    on_diag = torch.diagonal(c).add(-1).pow(2).sum()
    off_diag = (c - torch.diag(torch.diagonal(c))).pow(2).sum()
    return on_diag + 0.005 * off_diag


def build_neighbors(data_np: np.ndarray, k: int, pca_dim: int, seed: int) -> Optional[np.ndarray]:
    if k <= 0 or data_np.shape[0] <= 1:
        return None
    n_neighbors = min(k + 1, data_np.shape[0])
    dim = min(max(2, pca_dim), data_np.shape[1], data_np.shape[0] - 1)
    if dim < min(data_np.shape):
        emb = PCA(n_components=dim, random_state=seed).fit_transform(data_np)
    else:
        emb = data_np
    emb = normalize(emb, axis=1)
    nn = NearestNeighbors(n_neighbors=n_neighbors, metric="cosine")
    nn.fit(emb)
    indices = nn.kneighbors(emb, return_distance=False)
    return indices[:, 1:].astype(np.int64)


def neighbor_context(data_np: np.ndarray, neighbors: Optional[np.ndarray], batch_indices: np.ndarray, device: torch.device) -> Optional[torch.Tensor]:
    if neighbors is None or neighbors.shape[1] == 0:
        return None
    picked = neighbors[batch_indices]
    ctx = data_np[picked].mean(axis=1).astype(np.float32)
    return torch.as_tensor(ctx, dtype=torch.float32, device=device)


def update_ema(teacher: FlexibleScMAE, student: FlexibleScMAE, decay: float = 0.99) -> None:
    with torch.no_grad():
        for t_param, s_param in zip(teacher.parameters(), student.parameters()):
            t_param.data.mul_(decay).add_(s_param.data, alpha=1.0 - decay)


def loss_step(
    model: FlexibleScMAE,
    teacher: Optional[FlexibleScMAE],
    x: torch.Tensor,
    mask_prob: torch.Tensor,
    config: VariantConfig,
    token_bins: torch.Tensor,
    context: Optional[torch.Tensor],
    epoch: int,
    epochs: int,
) -> Tuple[torch.Tensor, Dict[str, torch.Tensor]]:
    x_corrupted, mask = apply_variant_noise(x, mask_prob, config.mask_strategy)
    latent, mask_logits, reconstruction, token_logits, gene_features = model.forward_mask(x_corrupted, context)
    rec = weighted_reconstruction_loss(reconstruction, x, mask, config)
    mask_loss = F.binary_cross_entropy_with_logits(mask_logits, mask)
    total = (1.0 - config.mask_loss_weight) * rec + config.mask_loss_weight * mask_loss
    parts: Dict[str, torch.Tensor] = {
        "rec": rec.detach(),
        "mask_loss": mask_loss.detach(),
        "mask_rate": mask.mean().detach(),
    }

    if config.token_weight > 0.0:
        tokens = tokenize_tensor(x, token_bins)
        flat_weight = (mask * config.masked_data_weight + (1.0 - mask) * (1.0 - config.masked_data_weight)).reshape(-1)
        ce = F.cross_entropy(token_logits.reshape(-1, config.token_bins), tokens.reshape(-1), reduction="none")
        tok_loss = (ce * flat_weight).sum() / flat_weight.sum().clamp_min(1e-8)
        total = total + config.token_weight * tok_loss
        parts["token"] = tok_loss.detach()

    if config.gene_feature_weight > 0.0:
        gf_loss = F.smooth_l1_loss(gene_features, target_gene_features(x))
        total = total + config.gene_feature_weight * gf_loss
        parts["gene_feature"] = gf_loss.detach()

    if config.consistency_weight > 0.0 or config.barlow_weight > 0.0:
        x_corrupted2, _ = apply_variant_noise(x, mask_prob, config.mask_strategy)
        latent2 = model.forward_mask(x_corrupted2, context)[0]
        if config.consistency_weight > 0.0:
            cons = F.mse_loss(F.normalize(latent, dim=1), F.normalize(latent2, dim=1))
            total = total + config.consistency_weight * cons
            parts["consistency"] = cons.detach()
        if config.barlow_weight > 0.0:
            bl = barlow_loss(latent, latent2)
            total = total + config.barlow_weight * bl
            parts["barlow"] = bl.detach()

    if teacher is not None and config.teacher_weight > 0.0:
        with torch.no_grad():
            teacher_latent = teacher.forward_mask(x, context)[0]
        distill = F.mse_loss(F.normalize(latent, dim=1), F.normalize(teacher_latent, dim=1))
        total = total + config.teacher_weight * distill
        parts["teacher"] = distill.detach()

    if context is not None and (config.graph_weight > 0.0 or config.retrieval_weight > 0.0):
        ctx_latent = model.forward_mask(context, None)[0].detach()
        graph = 1.0 - F.cosine_similarity(latent, ctx_latent, dim=1).mean()
        weight = config.graph_weight + config.retrieval_weight
        total = total + weight * graph
        parts["graph"] = graph.detach()

    if config.prototype_weight > 0.0:
        proto = prototype_loss(model, latent)
        total = total + config.prototype_weight * proto
        parts["prototype"] = proto.detach()

    if config.fuzzy_weight > 0.0:
        fz = fuzzy_loss(model, latent)
        total = total + config.fuzzy_weight * fz
        parts["fuzzy"] = fz.detach()

    if config.gate_l1_weight > 0.0:
        regs = model.regularization_terms()
        if "gene_gate_l1" in regs:
            total = total + config.gate_l1_weight * regs["gene_gate_l1"]
            parts["gene_gate_l1"] = regs["gene_gate_l1"].detach()

    if not torch.isfinite(total):
        raise FloatingPointError(f"Non-finite loss in {config.key}")
    parts["total"] = total.detach()
    return total, parts


def run_smoke(config: VariantConfig, device: torch.device) -> None:
    model = FlexibleScMAE(
        num_genes=64,
        hidden_size=32,
        dropout=config.dropout,
        encoder_kind=config.encoder_kind,
        token_bins=config.token_bins,
        n_prototypes=6,
        use_gene_gate=config.gate_l1_weight > 0.0,
    ).to(device)
    teacher = copy.deepcopy(model).eval() if config.teacher_weight > 0.0 else None
    x = torch.randn(16, 64, device=device)
    profile = build_profile(x.detach().cpu().numpy())
    prob = mask_probability(config, profile, config.mask_ratio, 1, 2)
    bins = torch.as_tensor(build_token_bins(x.detach().cpu().numpy(), config.token_bins), device=device)
    context = torch.randn_like(x) if (config.graph_weight + config.retrieval_weight) > 0.0 else None
    loss, _ = loss_step(model, teacher, x, prob, config, bins, context, 1, 2)
    loss.backward()
    grads = [p.grad for p in model.parameters() if p.requires_grad and p.grad is not None]
    if not grads:
        raise RuntimeError(f"No gradients produced for {config.key}")
    emb = model.feature(x)
    if emb.shape != (16, 32):
        raise RuntimeError(f"Unexpected embedding shape for {config.key}: {tuple(emb.shape)}")


def main(variant_key: str) -> None:
    args = parse_args(variant_key)
    base_cfg = get_variant(variant_key)
    config = base_cfg.with_cli_overrides(
        dropout=args.dropout,
        mask_ratio=args.mask_ratio,
        masked_data_weight=args.masked_data_weight,
        mask_loss_weight=args.mask_loss_weight,
        neighbor_k=args.neighbor_k,
    )
    family.set_seed(args.seed)
    save_dir = Path(ensure_dir(args.save_dir))
    save_json(vars(args), str(save_dir / "args.json"))
    save_json(config.__dict__, str(save_dir / "variant_config.json"))

    device = family.get_device(args.gpu, args.no_cuda)
    if args.smoke:
        run_smoke(config, device)
        save_json({"status": "smoke_pass", "variant": config.key}, str(save_dir / "smoke.json"))
        print(f"smoke_pass {config.key}")
        return

    start_time = time.time()
    dataset_name = args.dataset_name or Path(args.data_path).stem
    resolved_data_path = resolve_input_path(args.data_path, dataset_name)
    bundle = family.load_scmae_dataset(
        file_path=resolved_data_path,
        input_mode=args.input_mode,
        n_top_genes=args.n_top_genes,
        target_sum=args.target_sum,
        scale_input=args.scale_input,
        label_key=args.label_key,
        seed=args.seed,
    )
    save_json(bundle.profile, str(save_dir / "dataset_profile.json"))
    save_json(bundle.preprocess_config, str(save_dir / "preprocess_config.json"))
    data_np = bundle.data.astype(np.float32, copy=False)
    labels = bundle.labels
    n_clusters = int(args.n_clusters) if int(args.n_clusters) > 0 else int(len(np.unique(labels)))
    with open(save_dir / "selected_genes.txt", "w", encoding="utf-8") as handle:
        for gene in bundle.gene_names:
            handle.write(f"{gene}\n")

    profile = build_profile(data_np)
    token_bins_np = build_token_bins(data_np, config.token_bins)
    token_bins = torch.as_tensor(token_bins_np, dtype=torch.float32, device=device)
    use_neighbors = (config.graph_weight + config.retrieval_weight) > 0.0
    neighbors = build_neighbors(data_np, config.neighbor_k if use_neighbors else 0, args.knn_pca_dim, args.seed)

    dataset = family.IndexedExpressionDataset(data_np, labels)
    generator = torch.Generator()
    generator.manual_seed(args.seed)
    train_loader = DataLoader(dataset, batch_size=args.batch_size, shuffle=True, drop_last=True, generator=generator)
    eval_loader = DataLoader(dataset, batch_size=max(args.batch_size * 4, 512), shuffle=False, drop_last=False)

    model = FlexibleScMAE(
        num_genes=data_np.shape[1],
        hidden_size=args.hidden_size,
        dropout=config.dropout,
        encoder_kind=config.encoder_kind,
        token_bins=config.token_bins,
        n_prototypes=max(2, n_clusters),
        use_gene_gate=config.gate_l1_weight > 0.0,
    ).to(device)
    teacher = copy.deepcopy(model).eval() if config.teacher_weight > 0.0 else None
    optimizer = torch.optim.AdamW(model.parameters(), lr=args.lr, weight_decay=args.weight_decay)

    history = {
        "loss": [],
        "effective_mask_rate": [],
        "variant": config.key,
        "mask_semantics": "1 means actually replaced/corrupted and numerically changed",
    }

    for epoch in range(1, max(1, args.epochs) + 1):
        model.train()
        sums: Dict[str, float] = {}
        n_batches = 0
        prob = mask_probability(config, profile, config.mask_ratio, epoch, args.epochs)
        for batch_indices, x_cpu, _ in train_loader:
            x = x_cpu.to(device)
            ctx = neighbor_context(data_np, neighbors, batch_indices.numpy(), device) if use_neighbors else None
            optimizer.zero_grad(set_to_none=True)
            loss, parts = loss_step(model, teacher, x, prob, config, token_bins, ctx, epoch, args.epochs)
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 5.0)
            optimizer.step()
            if teacher is not None:
                update_ema(teacher, model)
            for key, value in parts.items():
                sums[key] = sums.get(key, 0.0) + float(value.detach().cpu())
            n_batches += 1
        avg = {key: value / max(1, n_batches) for key, value in sums.items()}
        history["loss"].append(avg.get("total", 0.0))
        history["effective_mask_rate"].append(avg.get("mask_rate", 0.0))
        if epoch == 1 or epoch == args.epochs or epoch % 10 == 0:
            print(
                f"{config.key} epoch={epoch:03d}/{args.epochs} "
                f"loss={avg.get('total', 0.0):.4f} mask={avg.get('mask_rate', 0.0):.4f}",
                flush=True,
            )

    embedding, labels_out = family.extract_embedding(model, eval_loader, device)
    np.save(save_dir / "embedding_final.npy", embedding.astype(np.float32))
    np.save(save_dir / "embeddings_base.npy", embedding.astype(np.float32))
    np.save(save_dir / "labels.npy", labels_out.astype(np.int64))
    np.save(save_dir / "gene_names.npy", bundle.gene_names.astype(str))
    family.save_embedding_h5(save_dir / "embedding.h5", embedding, labels_out)
    save_json(history, str(save_dir / "training_history.json"))
    torch.save(
        {
            "model": model.state_dict(),
            "args": vars(args),
            "variant_config": config.__dict__,
            "gene_names": bundle.gene_names.astype(str),
        },
        save_dir / "model_checkpoint.pth",
    )

    result = None
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
                "variant": config.key,
                "rank": int(config.rank),
                "source_paper": config.source_paper,
                "mask_strategy": config.mask_strategy,
                "preprocessing": "scMAE_family",
            },
        )
        save_json(result["fixed"], str(save_dir / "metrics.json"))

    elapsed = time.time() - start_time
    summary = {
        "dataset": dataset_name,
        "method": args.method_name,
        "variant": config.key,
        "rank": int(config.rank),
        "source_paper": config.source_paper,
        "seed": int(args.seed),
        "n_cells": int(data_np.shape[0]),
        "n_genes": int(data_np.shape[1]),
        "n_clusters": int(n_clusters),
        "runtime_seconds": float(elapsed),
        "fixed_metrics": result["fixed"] if result is not None else {},
    }
    save_json(summary, str(save_dir / "summary.json"))
    print(f"completed {config.key} dataset={dataset_name} seconds={elapsed:.1f}", flush=True)


def list_variants() -> None:
    for key, cfg in sorted(VARIANTS.items(), key=lambda item: item[1].rank):
        print(f"{cfg.rank:03d}\t{key}\t{cfg.title}")


def detect_free_gpu() -> int:
    try:
        result = subprocess.run(
            ["nvidia-smi", "--query-gpu=index,memory.used", "--format=csv,noheader,nounits"],
            text=True,
            capture_output=True,
            check=True,
            timeout=5,
        )
    except (subprocess.SubprocessError, FileNotFoundError):
        return 1
    candidates = []
    for line in result.stdout.splitlines():
        parts = [part.strip() for part in line.split(",")]
        if len(parts) != 2:
            continue
        idx = int(parts[0])
        if idx in {0, 7}:
            continue
        candidates.append((int(parts[1]), idx))
    if not candidates:
        return 1
    return sorted(candidates)[0][1]


def resolve_input_path(data_path: str, dataset_name: str) -> str:
    path = Path(data_path).resolve()
    if path.suffix.lower() != ".h5":
        return str(path)
    out_dir = SCMAES_DIR / "benchmark_data"
    out_dir.mkdir(parents=True, exist_ok=True)
    converted = out_dir / f"{dataset_name}.h5ad"
    if converted.exists():
        return str(converted)
    prepare_script = ROOT / "scripts" / "prepare_dataset.py"
    cmd = [
        sys.executable,
        str(prepare_script),
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
            "Failed to convert input .h5 to benchmark .h5ad.\n"
            f"Command: {' '.join(cmd)}\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"
        )
    return str(converted)
