#!/usr/bin/env python3
"""Gene-program bottleneck ablation on top of the DEC+std-floor winner.

Minimal falsifiable test (see README):
  A0  program_mode=none                       full-z clustering   -> DEC+floor baseline
  A1  program_mode=nmf  split_mode=none        full-z clustering   -> program as full-latent aux regularizer
  A1s program_mode=shuffled split_mode=none    full-z clustering   -> row-permuted target control
  A2f program_mode=nmf  split_mode=fixed       z_type clustering   -> split latent (mechanism reverse-test)
  A2x program_mode=nmf  split_mode=extra       z_type clustering   -> split but z_type keeps full width

The program head is an AUXILIARY loss; it never changes the recon/DEC targets.
NMF program target is built on UNSCALED log-expr (same source as recon target),
then z-scored per program column so lambda is not dominated by a few high-activity
programs. "shuffled" permutes target ROWS (cell<->program pairing) keeping every
program's marginal statistics identical.
"""
from __future__ import annotations

import argparse
import os
import random
import sys
import time
from pathlib import Path

for _var in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_var, "1")
os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")
os.environ.setdefault("NUMBA_CACHE_DIR", "/tmp/numba-cache")

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from sklearn.cluster import KMeans
from sklearn.decomposition import NMF, PCA
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import normalize
from torch.utils.data import DataLoader, Dataset

ROOT = Path(__file__).resolve().parents[2]
ADAPTIVE_SWITCH_DIR = ROOT / "experimental_retired_models" / "Granularity_scMAE_experiments" / "AdaptiveSwitch_scMAE"
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(ADAPTIVE_SWITCH_DIR) not in sys.path:
    sys.path.insert(0, str(ADAPTIVE_SWITCH_DIR))

from clusterability import compute_clusterability
from loss import AdaptiveSwitchLoss, compute_gate
from model import AdaptiveSwitchScMAE
from methods.DeepLearning import scMAE_family as family
from methods.shared_utils import ensure_dir, save_json


def _register_null_h5ad_reader() -> None:
    try:
        import h5py
        from anndata._io.specs.registry import _REGISTRY, IOSpec

        def _read_null(*args, **kwargs):
            return None

        for typ in (h5py.Dataset, h5py.Group):
            try:
                _REGISTRY.register_read(typ, IOSpec("null", "0.1.0"))(_read_null)
            except Exception:
                pass
    except Exception:
        pass


_register_null_h5ad_reader()


class ExprDataset(Dataset):
    def __init__(self, enc: np.ndarray, log: np.ndarray, labels: np.ndarray):
        self.enc = torch.as_tensor(enc, dtype=torch.float32)
        self.log = torch.as_tensor(log, dtype=torch.float32)
        self.labels = torch.as_tensor(labels, dtype=torch.long)

    def __len__(self) -> int:
        return int(self.enc.shape[0])

    def __getitem__(self, idx: int):
        return int(idx), self.enc[idx], self.log[idx], self.labels[idx]


class DecoderBypassScMAE(AdaptiveSwitchScMAE):
    """AdaptiveSwitch 骨干 (scMAE+DEC+std-floor)，唯一变量 = decoder 如何接收 mask。

    scMAE 原始设计里 decoder 输入是 [latent, mask_logits]，即把 G 维的 mask 预测
    整个喂给 decoder。假设：这个 G 维旁路让 encoder 不必在 latent 里编码「哪里被
    corrupt / 该细胞的表达结构」，对重构有利但可能让聚类 latent 变弱（结构捷径）。

    decoder_mode:
      concat : decoder([latent, mask_logits])           —— 原始 scMAE (D0)
      none   : decoder(latent)                          —— mask 只进 BCE loss，不进 decoder (D1)
      lowrank: decoder([latent, pool(mask_logits)->r])  —— mask 只以 r 维摘要进 decoder (D2)

    其它一切（encoder/mask_predictor/DEC centers/std-floor）与赢家完全一致，
    forward 返回的 latent = 全 128 维，KMeans/DEC/floor 照旧作用其上。
    mask_logits 始终返回，BCE mask 判别器 (loss.py) 在三种模式下都保留。
    """

    def __init__(self, num_genes, n_clusters, hidden_size=128, dropout=0.05,
                 decoder_mode="concat", mask_rank=16):
        super().__init__(num_genes, n_clusters, hidden_size, dropout)
        self.decoder_mode = str(decoder_mode)
        self.mask_rank = int(mask_rank)
        if self.decoder_mode == "concat":
            pass  # 继承的 self.decoder = Linear(hidden+genes, genes) 即原始
        elif self.decoder_mode == "none":
            self.decoder = nn.Linear(hidden_size, self.num_genes)
        elif self.decoder_mode == "lowrank":
            self.mask_pool = nn.Linear(self.num_genes, self.mask_rank)
            self.decoder = nn.Linear(hidden_size + self.mask_rank, self.num_genes)
        else:
            raise ValueError(f"unknown decoder_mode {decoder_mode}")

    def forward(self, x: torch.Tensor) -> dict[str, torch.Tensor]:
        latent = self.encoder(x)
        mask_logits = self.mask_predictor(latent)
        if self.decoder_mode == "concat":
            reconstruction = self.decoder(torch.cat([latent, mask_logits], dim=1))
        elif self.decoder_mode == "none":
            reconstruction = self.decoder(latent)
        else:  # lowrank
            u = self.mask_pool(mask_logits)
            reconstruction = self.decoder(torch.cat([latent, u], dim=1))
        q = self.student_q(latent)
        return {"latent": latent, "mask_logits": mask_logits,
                "reconstruction": reconstruction, "cluster_q": q}


def parse_args():
    p = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--data_path", required=True)
    p.add_argument("--save_dir", required=True)
    p.add_argument("--dataset_name", default=None)
    p.add_argument("--label_key", default="auto")
    p.add_argument("--input_mode", default="auto", choices=["auto", "raw", "log1p"])
    p.add_argument("--n_top_genes", type=int, default=1000)
    p.add_argument("--target_sum", type=float, default=10000.0)
    p.add_argument("--scale_input", type=family.str2bool, default=True)
    p.add_argument("--n_clusters", type=int, required=True)
    p.add_argument("--epochs", type=int, default=80)
    p.add_argument("--batch_size", type=int, default=256)
    p.add_argument("--lr", type=float, default=1e-3)
    p.add_argument("--weight_decay", type=float, default=1e-5)
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--gpu", type=int, default=1)
    p.add_argument("--no_cuda", action="store_true")
    p.add_argument("--skip_eval", action="store_true")
    p.add_argument("--hidden_size", type=int, default=128)
    p.add_argument("--dropout", type=float, default=0.05)
    p.add_argument("--mask_prob", type=float, default=0.4)
    p.add_argument("--masked_data_weight", type=float, default=0.75)
    p.add_argument("--mask_weight", type=float, default=0.65)
    p.add_argument("--cluster_weight", type=float, default=0.35)
    p.add_argument("--consistency_weight", type=float, default=0.05)
    p.add_argument("--variance_weight", type=float, default=0.02)
    p.add_argument("--var_mode", default="hinge", choices=["hinge", "cov", "koleo", "both"])
    p.add_argument("--entropy_weight", type=float, default=0.10)
    p.add_argument("--confidence_threshold", type=float, default=0.35)
    p.add_argument("--warmup_epochs", type=int, default=20)
    p.add_argument("--target_update_interval", type=int, default=5)
    p.add_argument("--neighbor_k", type=int, default=15)
    p.add_argument("--knn_pca_dim", type=int, default=50)
    p.add_argument("--gate_kappa", type=float, default=0.15)
    p.add_argument("--force_gate", type=float, default=1.0)
    p.add_argument("--gate_ema", type=float, default=0.5)
    # decoder-bypass ablation knobs
    p.add_argument("--decoder_mode", default="concat", choices=["concat", "none", "lowrank"])
    p.add_argument("--mask_rank", type=int, default=16)
    p.add_argument("--smoke", action="store_true")
    return p.parse_args()


def set_seed(seed: int) -> None:
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False


def get_device(gpu: int, no_cuda: bool) -> torch.device:
    if no_cuda or not torch.cuda.is_available():
        return torch.device("cpu")
    if gpu in {0, 7}:
        raise ValueError("Physical GPU 0 and GPU 7 are intentionally avoided.")
    return torch.device(f"cuda:{gpu}")


def build_neighbor_indices(data: np.ndarray, k: int, pca_dim: int, seed: int):
    n_cells, n_genes = data.shape
    k = max(1, min(int(k), n_cells - 1))
    dim = max(2, min(int(pca_dim), n_cells - 1, n_genes - 1))
    emb = PCA(n_components=dim, random_state=seed).fit_transform(data.astype(np.float64))
    emb = normalize(emb, norm="l2", axis=1)
    nn_ = NearestNeighbors(n_neighbors=k + 1, metric="cosine").fit(emb)
    _, neighbors = nn_.kneighbors(emb, return_distance=True)
    return neighbors[:, 1:].astype(np.int64)


@torch.no_grad()
def extract_all(model, loader, device, want_program=False):
    model.eval()
    emb, q_all, labels, prog_pred = [], [], [], []
    for _, x, _, y in loader:
        out = model(x.to(device))
        emb.append(out["latent"].detach().cpu().numpy())
        q_all.append(out["cluster_q"].detach().cpu().numpy())
        labels.append(y.numpy())
        if want_program and "program_pred" in out:
            prog_pred.append(out["program_pred"].detach().cpu().numpy())
    emb = np.nan_to_num(np.concatenate(emb).astype(np.float32))
    q_all = np.concatenate(q_all).astype(np.float32)
    labels = np.concatenate(labels).astype(np.int64)
    pp = np.concatenate(prog_pred).astype(np.float32) if prog_pred else None
    return emb, q_all, labels, pp


def effective_dimensionality(std: np.ndarray) -> dict:
    var = np.square(std.astype(np.float64))
    pr = float((var.sum() ** 2) / max(float(np.square(var).sum()), 1e-12))
    return {
        "std_min": float(std.min()),
        "std_median": float(np.median(std)),
        "std_max": float(std.max()),
        "effective_dim_pr": pr,
        "dims_std_gt_0p1": int((std > 0.1).sum()),
        "dims_std_gt_1p0": int((std > 1.0).sum()),
    }


def cluster_aligned_eff_dim(emb: np.ndarray, labels: np.ndarray) -> dict:
    """Participation ratio of the between-class scatter eigenspectrum.

    Measures how many latent dims actually carry cell-type-discriminative signal
    (vs eff_dim_pr which just measures how many dims are active at all).
    """
    emb = emb.astype(np.float64)
    grand = emb.mean(axis=0, keepdims=True)
    classes = np.unique(labels)
    d = emb.shape[1]
    Sb = np.zeros((d, d), dtype=np.float64)
    for c in classes:
        m = labels == c
        n_c = int(m.sum())
        if n_c == 0:
            continue
        diff = emb[m].mean(axis=0, keepdims=True) - grand
        Sb += n_c * (diff.T @ diff)
    ev = np.linalg.eigvalsh(Sb).clip(min=0.0)
    pr = float((ev.sum() ** 2) / max(float(np.square(ev).sum()), 1e-12))
    return {"cluster_aligned_eff_dim": pr, "between_class_scatter_trace": float(np.trace(Sb))}


def main() -> int:
    args = parse_args()
    if args.smoke:
        args.epochs = min(args.epochs, 3)
        args.warmup_epochs = min(args.warmup_epochs, 1)
    set_seed(args.seed)
    save_dir = Path(ensure_dir(args.save_dir))
    save_json(vars(args), str(save_dir / "args.json"))
    device = get_device(args.gpu, args.no_cuda)
    dataset_name = args.dataset_name or Path(args.data_path).stem

    target_bundle = family.load_scmae_dataset(
        args.data_path, args.input_mode, args.n_top_genes, args.target_sum, False, args.label_key, args.seed
    )
    if args.scale_input:
        encoder_bundle = family.load_scmae_dataset(
            args.data_path, args.input_mode, args.n_top_genes, args.target_sum, True, args.label_key, args.seed
        )
        encoder_data = np.asarray(encoder_bundle.data, dtype=np.float32)
    else:
        encoder_data = np.asarray(target_bundle.data, dtype=np.float32)
    log_expr = np.asarray(target_bundle.data, dtype=np.float32)
    labels = np.asarray(target_bundle.labels, dtype=np.int64)
    n_clusters = int(args.n_clusters if args.n_clusters > 0 else len(np.unique(labels)))
    save_json(target_bundle.profile, str(save_dir / "dataset_profile.json"))

    nb_indices = build_neighbor_indices(encoder_data, args.neighbor_k, args.knn_pca_dim, args.seed)

    dataset = ExprDataset(encoder_data, log_expr, labels)
    generator = torch.Generator().manual_seed(args.seed)
    train_loader = DataLoader(dataset, batch_size=args.batch_size, shuffle=True, drop_last=True, generator=generator)
    full_loader = DataLoader(dataset, batch_size=max(args.batch_size * 4, 512), shuffle=False, drop_last=False)

    model = DecoderBypassScMAE(
        encoder_data.shape[1], n_clusters, args.hidden_size, args.dropout,
        decoder_mode=args.decoder_mode, mask_rank=args.mask_rank,
    ).to(device)
    criterion = AdaptiveSwitchLoss(
        args.masked_data_weight, args.mask_weight, args.cluster_weight,
        args.consistency_weight, args.variance_weight, args.entropy_weight,
        args.confidence_threshold, args.var_mode,
    )
    optimizer = torch.optim.Adam(model.parameters(), lr=args.lr, weight_decay=args.weight_decay)

    p_targets = None
    clusterab = np.ones(encoder_data.shape[0], dtype=np.float32)
    gate = float(args.force_gate)
    kl_ref = 0.0
    centers_initialized = False
    history = {k: [] for k in ["loss", "base_loss", "sharp_loss", "variance_loss", "gate", "kl_ref"]}
    start = time.time()
    print(
        f"Device={device} dataset={dataset_name} cells={encoder_data.shape[0]} clusters={n_clusters} "
        f"decoder_mode={args.decoder_mode} mask_rank={args.mask_rank} varw={args.variance_weight}"
    )

    for epoch in range(1, max(1, args.epochs) + 1):
        if epoch > args.warmup_epochs and (
            (epoch - args.warmup_epochs - 1) % max(1, args.target_update_interval) == 0 or p_targets is None
        ):
            emb, q_full, _, _ = extract_all(model, full_loader, device)
            if not centers_initialized:
                km = KMeans(n_clusters=n_clusters, n_init=20, random_state=args.seed).fit(emb)
                model.initialize_centers(torch.as_tensor(km.cluster_centers_, dtype=torch.float32, device=device))
                emb, q_full, _, _ = extract_all(model, full_loader, device)
                centers_initialized = True
            sharp_p = AdaptiveSwitchScMAE.sharpen(torch.as_tensor(q_full)).numpy().astype(np.float32)
            p_targets = sharp_p
            clusterab, _ = compute_clusterability(emb, q_full, nb_indices, k=args.neighbor_k)
            g_new, kl_ref = compute_gate(sharp_p, q_full, args.gate_kappa)
            gate = args.gate_ema * gate + (1.0 - args.gate_ema) * g_new if args.force_gate < 0 else float(args.force_gate)

        cluster_scale = 0.0 if epoch <= args.warmup_epochs else min(1.0, (epoch - args.warmup_epochs) / max(1, args.warmup_epochs))
        model.train()
        sums = {k: 0.0 for k in ["loss", "base_loss", "sharp_loss", "variance_loss"]}
        batches = 0
        for idx, x_cpu, log_cpu, _ in train_loader:
            idx_np = idx.numpy().astype(np.int64, copy=False)
            x = x_cpu.to(device)
            target = log_cpu.to(device)
            c_batch = torch.as_tensor(clusterab[idx_np], dtype=torch.float32, device=device)
            strong, mask = model.random_mask(x, args.mask_prob)
            weak, _ = model.random_mask(x, max(0.05, args.mask_prob * 0.5))
            out = model(strong)
            weak_out = model(weak)
            p_batch = None if p_targets is None else torch.as_tensor(p_targets[idx_np], dtype=torch.float32, device=device)
            base_loss, parts = criterion(out, weak_out, target, mask, p_batch, c_batch, cluster_scale, gate)

            loss = base_loss
            if not torch.isfinite(loss):
                raise FloatingPointError(f"Non-finite loss at epoch {epoch}: {parts}")
            optimizer.zero_grad(set_to_none=True)
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 5.0)
            optimizer.step()

            sums["loss"] += float(loss.detach().cpu())
            sums["base_loss"] += float(base_loss.detach().cpu())
            sums["sharp_loss"] += float(parts["sharp_loss"])
            sums["variance_loss"] += float(parts["variance_loss"])
            batches += 1
        for key in sums:
            history[key].append(sums[key] / max(1, batches))
        history["gate"].append(float(gate))
        history["kl_ref"].append(float(kl_ref))
        if epoch == 1 or epoch == args.epochs or epoch % 10 == 0:
            print(
                f"Epoch {epoch:03d}/{args.epochs} loss={history['loss'][-1]:.4f} "
                f"base={history['base_loss'][-1]:.4f} "
                f"sharp={history['sharp_loss'][-1]:.4f} var={history['variance_loss'][-1]:.4f} gate={gate:.3f}"
            )

    embedding, q_out, labels_out, _ = extract_all(model, full_loader, device)
    np.save(save_dir / "embedding_final.npy", embedding.astype(np.float32))
    np.save(save_dir / "labels.npy", labels_out.astype(np.int64))
    save_json(history, str(save_dir / "training_history.json"))
    std_profile = effective_dimensionality(embedding.std(axis=0))
    aligned = cluster_aligned_eff_dim(embedding, labels_out)

    eval_result = None
    preds = None
    if not args.skip_eval:
        eval_result = family.write_kmeans_known_k_outputs(
            save_dir, dataset_name, "decoder-bypass ablation", args.seed,
            embedding, labels_out, n_clusters,
            {"decoder_mode": args.decoder_mode, "mask_rank": int(args.mask_rank),
             "variance_weight": float(args.variance_weight)},
        )
        preds = eval_result["preds"]["kmeans_known_k"]
        save_json(eval_result["fixed"], str(save_dir / "metrics.json"))

    counts = np.bincount(preds.astype(np.int64), minlength=n_clusters).astype(np.float64) if preds is not None else np.zeros(n_clusters)
    frac = counts / max(1.0, counts.sum())
    summary = {
        "dataset": dataset_name,
        "seed": int(args.seed),
        "n_clusters": int(n_clusters),
        "runtime_seconds": float(time.time() - start),
        "decoder_mode": args.decoder_mode,
        "mask_rank": int(args.mask_rank),
        "variance_weight": float(args.variance_weight),
        "fixed_metrics": eval_result["fixed"] if eval_result is not None else {},
        "std_profile": std_profile,
        "cluster_aligned": aligned,
        "final_base_loss": float(history["base_loss"][-1]) if history["base_loss"] else 0.0,
        "cluster_mass_min": float(frac.min()) if frac.size else 0.0,
        "cluster_mass_max": float(frac.max()) if frac.size else 0.0,
        "final_gate": float(gate),
    }
    save_json(summary, str(save_dir / "summary.json"))
    ari = summary["fixed_metrics"].get("kmeans_known_k", {}).get("ari")
    nmi = summary["fixed_metrics"].get("kmeans_known_k", {}).get("nmi")
    print(
        f"[RESULT] {dataset_name} decoder_mode={args.decoder_mode} "
        f"ARI={ari} NMI={nmi} aligned_eff_dim={aligned['cluster_aligned_eff_dim']:.1f}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
