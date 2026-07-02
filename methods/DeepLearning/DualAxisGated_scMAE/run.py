#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os
import sys
import time
from pathlib import Path

# Cap BLAS/OMP threads BEFORE importing numpy/sklearn. On large datasets
# (e.g. Macosko 44k cells) KMeans(n_init=20) + NearestNeighbors otherwise
# spawn enough thread-local buffers to trigger
# "OpenBLAS: too many memory regions" and kill the process during eval.
for _var in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_var, "8")

import numpy as np
import torch

from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors
from torch.utils.data import DataLoader, Dataset

CURRENT_DIR = Path(__file__).resolve().parent
ROOT = next(parent for parent in [CURRENT_DIR, *CURRENT_DIR.parents] if (parent / "methods" / "DeepLearning" / "scMAE_family.py").exists())
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(CURRENT_DIR) not in sys.path:
    sys.path.insert(0, str(CURRENT_DIR))
os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")
os.environ.setdefault("NUMBA_CACHE_DIR", "/tmp/numba-cache")

from loss import GatedFusionLoss
from model import DualAxisGatedScMAE
from reliability import compute_reliability
from methods.DeepLearning import scMAE_family as family
from methods.shared_utils import ensure_dir, save_json


def _register_null_h5ad_reader() -> None:
    """The benchmark_data h5ad files store some uns entries with encoding-type
    'null' (a None value) that installed anndata versions can write but not
    read back. Register a tolerant reader so sc.read_h5ad does not crash. This
    only affects THIS process; it does not modify shared code or data files.
    """
    try:
        import h5py
        from anndata._io.specs.registry import _REGISTRY, IOSpec

        def _read_null(*args, **kwargs):
            return None

        for typ in (h5py.Dataset, h5py.Group):
            for ver in ("0.1.0",):
                try:
                    _REGISTRY.register_read(typ, IOSpec("null", ver))(_read_null)
                except Exception:
                    pass
    except Exception:
        pass


_register_null_h5ad_reader()


METHOD_NAME = "DualAxisGated_scMAE"
DISPLAY_NAME = "dual-axis (gene-module + prototype cell-axis) gated scMAE + DEC"
BASELINES = {
    "Melanoma_5K": {"nmi": 0.735414, "ari": 0.668029},
    "Quake_10x_Spleen": {"nmi": 0.851730, "ari": 0.922275},
    "Macosko": {"nmi": 0.657465, "ari": 0.494268},
}


class ExprDataset(Dataset):
    def __init__(self, encoder_data: np.ndarray, log_expr: np.ndarray, labels: np.ndarray):
        self.encoder_data = torch.as_tensor(encoder_data, dtype=torch.float32)
        self.log_expr = torch.as_tensor(log_expr, dtype=torch.float32)
        self.labels = torch.as_tensor(labels, dtype=torch.long)

    def __len__(self) -> int:
        return int(self.encoder_data.shape[0])

    def __getitem__(self, idx: int):
        return int(idx), self.encoder_data[idx], self.log_expr[idx], self.labels[idx]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Reliability-gated NeighborMix + DEC scMAE candidate.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
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
    p.add_argument("--smoke", action="store_true")
    p.add_argument("--hidden_size", type=int, default=128)
    p.add_argument("--dropout", type=float, default=0.05)
    p.add_argument("--mask_prob", type=float, default=0.4)
    p.add_argument("--masked_data_weight", type=float, default=0.75)
    p.add_argument("--mask_weight", type=float, default=0.65)
    p.add_argument("--cluster_weight", type=float, default=0.35)
    p.add_argument("--confidence_threshold", type=float, default=0.35)
    p.add_argument("--consistency_weight", type=float, default=0.05)
    p.add_argument("--warmup_epochs", type=int, default=20)
    p.add_argument("--target_update_interval", type=int, default=5)
    # NeighborMix + gate hyperparameters
    p.add_argument("--pseudo_weight", type=float, default=0.3)
    p.add_argument("--alpha_min", type=float, default=0.6, help="min self-weight for fully-reliable cells; alpha_i = 1 - (1-alpha_min)*r_i")
    p.add_argument("--mix_neighbors", type=int, default=4)
    p.add_argument("--neighbor_k", type=int, default=15)
    p.add_argument("--knn_pca_dim", type=int, default=50)
    p.add_argument("--reliability_update_interval", type=int, default=5)
    # Dual-axis encoder hyperparameters
    p.add_argument("--n_gene_modules", type=int, default=64, help="gene modules M for the gene axis (shared gene-row weights)")
    p.add_argument("--token_dim", type=int, default=48)
    p.add_argument("--n_prototypes", type=int, default=16, help="K prototypes for the scalable cell axis")
    p.add_argument("--gene_layers", type=int, default=2)
    p.add_argument("--attn_heads", type=int, default=4)
    p.add_argument("--attn_dropout", type=float, default=0.1)
    p.add_argument("--use_mlp_base", type=family.str2bool, default=True, help="fuse a full-gene MLP base with the dual-axis signal (guarantees scMAE-level fidelity)")
    p.add_argument("--mlp_hidden", type=int, default=256)
    p.add_argument("--encoder_type", type=str, default="dualaxis", choices=["dualaxis", "cellaxis_on_mlp", "mlp"], help="cellaxis_on_mlp: drop gene-module axis, run gated prototype cell-axis on MLP features (zero-init gate)")
    # Marker-aware masking (anti-smoothing): preserve/emphasize bimodal marker genes
    p.add_argument("--marker_recon_strength", type=float, default=0.0, help="per-gene reconstruction weight = 1 + strength*(imp-mean); 0 = uniform (scMAE)")
    p.add_argument("--marker_mask_bias", type=float, default=0.0, help="tilt mask probability toward high-marker genes; 0 = uniform masking")
    p.add_argument("--eval_gate", type=family.str2bool, default=False, help="if true, apply the reliability gate at eval too; default false = full encoder at eval")
    return p.parse_args()


def build_gene_marker_weight(data_np: np.ndarray) -> np.ndarray:
    """Label-free per-gene 'marker-ness' via the bimodality coefficient
    BC = (skew^2 + 1) / kurtosis. Marker genes are bimodal across cells (on in
    one cell type, off elsewhere) -> high BC; housekeeping genes are unimodal ->
    low BC. Returned as a mean-normalized weight (~1.0 average) so it can tilt
    masking / reconstruction toward markers without changing overall scale.
    """
    x = np.asarray(data_np, dtype=np.float64)
    n = max(2, x.shape[0])
    mu = x.mean(axis=0)
    xc = x - mu
    var = (xc ** 2).mean(axis=0)
    std = np.sqrt(var) + 1e-8
    skew = (xc ** 3).mean(axis=0) / (std ** 3)
    kurt = (xc ** 4).mean(axis=0) / (std ** 4)          # Pearson (normal=3)
    bc = (skew ** 2 + 1.0) / np.clip(kurt, 1e-6, None)   # Sarle's bimodality coeff
    bc = np.nan_to_num(bc, nan=0.0, posinf=0.0, neginf=0.0)
    # rank to [0,1] percentile for robustness, then mean-normalize to ~1.0
    order = np.argsort(np.argsort(bc)).astype(np.float64)
    imp = order / max(1.0, len(order) - 1.0)             # [0,1]
    imp = imp / max(imp.mean(), 1e-8)                    # mean ~1
    return imp.astype(np.float32)


def build_gene_modules(data_np: np.ndarray, n_modules: int, seed: int) -> np.ndarray:
    """Hard [G, M] gene->module assignment via KMeans on gene profiles.

    Label-free: clusters genes by their expression pattern across cells. Gives
    the dual-axis encoder its shared gene-row weights while keeping the gene-axis
    token count small (M<<G) so attention never blows up.
    """
    g = int(data_np.shape[1])
    M = max(1, min(int(n_modules), g))
    if M >= g:
        ids = np.arange(g, dtype=np.int64)
        M = g
    else:
        ids = KMeans(n_clusters=M, n_init=10, random_state=seed).fit_predict(data_np.T.astype(np.float64)).astype(np.int64)
    assignment = np.zeros((g, M), dtype=np.float32)
    assignment[np.arange(g), ids] = 1.0
    return assignment


def build_knn_graph(data_np: np.ndarray, k: int, pca_dim: int, seed: int):
    """Static PCA-KNN graph (cosine) as in NeighborMix. Returns (indices, probs)."""
    n_cells = int(data_np.shape[0])
    max_k = min(int(k), max(1, n_cells - 1))
    if max_k <= 0:
        return None, None
    dim = min(int(pca_dim), min(data_np.shape) - 1)
    emb = PCA(n_components=max(2, dim), random_state=seed).fit_transform(data_np.astype(np.float64)) if dim >= 2 else data_np.astype(np.float64)
    emb = emb / np.clip(np.linalg.norm(emb, axis=1, keepdims=True), 1e-12, None)
    nn = NearestNeighbors(n_neighbors=max_k + 1, metric="cosine").fit(emb)
    dist, idx = nn.kneighbors(emb)
    indices = idx[:, 1:max_k + 1].astype(np.int64)
    sim = np.clip(1.0 - dist[:, 1:max_k + 1], 1e-6, None).astype(np.float32)
    probs = sim / sim.sum(axis=1, keepdims=True).clip(1e-12)
    return indices, probs.astype(np.float32)


def neighbor_mix_view(data_np, batch_idx, batch_x, indices, probs, r_batch, alpha_min, mix_neighbors, rng):
    """Per-cell reliability-gated NeighborMix. alpha_i = 1 - (1-alpha_min)*r_i.

    r_i≈1 -> alpha≈alpha_min (strong mix). r_i≈0 -> alpha≈1 (mixing OFF).
    """
    if indices is None or mix_neighbors <= 0:
        return batch_x
    bsz = int(batch_idx.shape[0])
    m = min(int(mix_neighbors), indices.shape[1])
    sampled = np.empty((bsz, m), dtype=np.int64)
    weights = np.empty((bsz, m), dtype=np.float32)
    for pos, cell in enumerate(batch_idx):
        pr = probs[cell]
        choice = rng.choice(indices.shape[1], size=m, replace=True, p=pr)
        sampled[pos] = indices[cell, choice]
        w = pr[choice].astype(np.float32)
        weights[pos] = w / max(float(w.sum()), 1e-12)
    neighbor_mean = np.sum(data_np[sampled] * weights[:, :, None], axis=1).astype(np.float32)
    neighbor_t = torch.as_tensor(neighbor_mean, dtype=batch_x.dtype, device=batch_x.device)
    alpha = (1.0 - (1.0 - float(alpha_min)) * r_batch).unsqueeze(1)  # [bsz,1]
    return alpha * batch_x + (1.0 - alpha) * neighbor_t


@torch.no_grad()
def extract_all(model, loader, device, reliability=None):
    """Extract latent/q/labels. If reliability is given (indexed by cell id),
    the gate is applied at inference; otherwise the full encoder is used."""
    model.eval()
    emb, q_all, labels = [], [], []
    for idx, x, _, y in loader:
        r = None
        if reliability is not None:
            r = torch.as_tensor(reliability[idx.numpy()], dtype=torch.float32, device=device)
        out = model(x.to(device), r)
        emb.append(out["latent"].detach().cpu().numpy())
        q_all.append(out["cluster_q"].detach().cpu().numpy())
        labels.append(y.numpy())
    return (np.nan_to_num(np.concatenate(emb).astype(np.float32)),
            np.concatenate(q_all).astype(np.float32),
            np.concatenate(labels).astype(np.int64))


def metric(eval_result, name):
    if not eval_result:
        return None
    v = eval_result.get("fixed", {}).get("kmeans_known_k", {}).get(name)
    return None if v is None else float(v)


def diagnostics(embedding, labels, q, n_clusters, seed, preds, rel_diag, mixed_fraction):
    if preds is None:
        preds = KMeans(n_clusters=n_clusters, n_init=20, random_state=seed).fit_predict(embedding)
    counts = np.bincount(preds.astype(np.int64), minlength=n_clusters).astype(np.float64)
    frac = counts / max(1.0, counts.sum())
    var = float(np.var(embedding, axis=0).mean()) if embedding.size else 0.0
    nn = NearestNeighbors(n_neighbors=min(11, embedding.shape[0])).fit(embedding)
    idx = nn.kneighbors(embedding, return_distance=False)[:, 1:]
    purity = float(np.mean(labels[idx] == labels[:, None]))
    d = {
        "neighbor_purity_proxy": purity,
        "mixed_cell_fraction": float(mixed_fraction),
        "embedding_variance": var,
        "cluster_mass_min": float(frac.min()) if frac.size else 0.0,
        "cluster_mass_max": float(frac.max()) if frac.size else 0.0,
        "collapse_warning": bool((not np.isfinite(var)) or var < 1e-8 or (frac.min() if frac.size else 0) < 0.001 or (frac.max() if frac.size else 1) > 0.95),
    }
    d.update(rel_diag)
    return d


def main() -> int:
    args = parse_args()
    if args.gpu in {0, 7} and not args.no_cuda:
        raise ValueError("GPU 0 and GPU 7 are forbidden. Choose GPU 1-6 or --no_cuda.")
    if args.smoke:
        args.epochs = min(args.epochs, 3)
        args.warmup_epochs = min(args.warmup_epochs, 1)
    family.set_seed(args.seed)
    save_dir = Path(ensure_dir(args.save_dir))
    save_json(vars(args), str(save_dir / "args.json"))
    device = family.get_device(args.gpu, args.no_cuda)
    dataset_name = args.dataset_name or Path(args.data_path).stem
    stage = "smoke" if args.smoke else "screen"
    rng = np.random.default_rng(args.seed)

    # Encoder input = SCALED expression; reconstruction target = UNSCALED log
    # expression. This split is intentional and critical (matches rank13): the
    # unscaled target keeps the reconstruction loss magnitude large (~0.4) so
    # the DEC KL term stays a gentle regularizer instead of dominating and
    # distorting the embedding. Reconstructing scaled->scaled drops scmae loss
    # ~6x and lets DEC take over -> collapse (verified on Quake: 0.47 vs 0.91).
    target_bundle = family.load_scmae_dataset(args.data_path, args.input_mode, args.n_top_genes, args.target_sum, False, args.label_key, args.seed)
    if args.scale_input:
        encoder_bundle = family.load_scmae_dataset(args.data_path, args.input_mode, args.n_top_genes, args.target_sum, True, args.label_key, args.seed)
        encoder_data = encoder_bundle.data
    else:
        encoder_data = target_bundle.data
    log_expr = np.asarray(target_bundle.data, dtype=np.float32)  # unscaled target
    labels = np.asarray(target_bundle.labels, dtype=np.int64)
    n_clusters = int(args.n_clusters if args.n_clusters > 0 else len(np.unique(labels)))
    save_json(target_bundle.profile, str(save_dir / "dataset_profile.json"))

    # NeighborMix graph is built on the ENCODER input (mix happens in encoder-input space)
    print("Building PCA-KNN graph for NeighborMix ...")
    nb_indices, nb_probs = build_knn_graph(encoder_data, args.neighbor_k, args.knn_pca_dim, args.seed)

    # Dual-axis gene modules (shared gene-row weights), built once, label-free.
    print(f"Building {args.n_gene_modules} gene modules for the gene axis ...")
    assignment = build_gene_modules(encoder_data, args.n_gene_modules, args.seed)
    np.save(save_dir / "gene_module_assignment.npy", assignment.astype(np.float32))

    # Label-free per-gene marker weight (bimodality) on the unscaled target space.
    gene_imp = build_gene_marker_weight(log_expr)          # mean ~1
    np.save(save_dir / "gene_marker_weight.npy", gene_imp)
    # Reconstruction weight: 1 + strength*(imp-mean). strength=0 -> all ones.
    recon_w_np = 1.0 + float(args.marker_recon_strength) * (gene_imp - gene_imp.mean())
    recon_w_np = np.clip(recon_w_np, 0.05, None).astype(np.float32)
    recon_weight = torch.as_tensor(recon_w_np, dtype=torch.float32, device=device) if args.marker_recon_strength > 0 else None

    dataset = ExprDataset(encoder_data, log_expr, labels)
    gen = torch.Generator().manual_seed(args.seed)
    train_loader = DataLoader(dataset, batch_size=args.batch_size, shuffle=True, drop_last=True, generator=gen)
    full_loader = DataLoader(dataset, batch_size=max(args.batch_size * 4, 512), shuffle=False, drop_last=False)

    model = DualAxisGatedScMAE(
        encoder_data.shape[1], n_clusters, assignment,
        hidden_size=args.hidden_size, dropout=args.dropout,
        token_dim=args.token_dim, n_prototypes=args.n_prototypes,
        heads=args.attn_heads, gene_layers=args.gene_layers, attn_dropout=args.attn_dropout,
        use_mlp_base=args.use_mlp_base, mlp_hidden=args.mlp_hidden,
        encoder_type=args.encoder_type, gene_weight=gene_imp,
    ).to(device)
    criterion = GatedFusionLoss(args.masked_data_weight, args.mask_weight, args.cluster_weight, args.consistency_weight, args.confidence_threshold)
    optimizer = torch.optim.Adam(model.parameters(), lr=args.lr, weight_decay=args.weight_decay)

    reliability = np.ones(encoder_data.shape[0], dtype=np.float32)  # start fully-reliable (pure NeighborMix) during warmup
    rel_diag = {"reliability_mean": 1.0, "core_fraction": 1.0, "agree_mean": 1.0, "density_mean": 1.0, "confidence_mean": 0.0, "reliability_min": 1.0, "reliability_max": 1.0}
    p_targets = None
    centers_initialized = False
    history = {"loss": [], "scmae_loss": [], "pseudo_loss": [], "cluster_loss": [], "consistency_loss": [], "confidence_fraction": [], "reliability_mean": [], "cluster_scale": [], "stage": stage}
    start = time.time()
    print(f"Device={device} dataset={dataset_name} cells={encoder_data.shape[0]} genes={encoder_data.shape[1]} clusters={n_clusters}")
    print(f"Method={METHOD_NAME} stage={stage} epochs={args.epochs} warmup={args.warmup_epochs}")

    for epoch in range(1, max(1, args.epochs) + 1):
        # After warmup: refresh embedding -> centers -> DEC target p AND reliability field r_i
        if epoch > args.warmup_epochs and ((epoch - args.warmup_epochs - 1) % max(1, args.target_update_interval) == 0 or p_targets is None):
            emb, q_full, _ = extract_all(model, full_loader, device)
            if not centers_initialized:
                km = KMeans(n_clusters=n_clusters, n_init=20, random_state=args.seed).fit(emb)
                model.initialize_centers(torch.as_tensor(km.cluster_centers_, dtype=torch.float32, device=device))
                emb, q_full, _ = extract_all(model, full_loader, device)
                centers_initialized = True
            p_targets = DualAxisGatedScMAE.target_distribution(torch.as_tensor(q_full)).numpy().astype(np.float32)
            reliability, rel_diag = compute_reliability(emb, q_full, nb_indices, k=args.neighbor_k)

        cluster_scale = 0.0 if epoch <= args.warmup_epochs else min(1.0, (epoch - args.warmup_epochs) / max(1, args.warmup_epochs))
        model.train()
        sums = {k: 0.0 for k in ["loss", "scmae_loss", "pseudo_loss", "cluster_loss", "consistency_loss", "confidence_fraction"]}
        n_batches = 0
        for idx, x_cpu, log_cpu, _ in train_loader:
            idx_np = idx.numpy()
            x = x_cpu.to(device)
            target = log_cpu.to(device)
            r_batch = torch.as_tensor(reliability[idx_np], dtype=torch.float32, device=device)

            # rank13 strong + weak zero-masked views (target = unscaled log_expr)
            strong, mask = model.random_mask(x, args.mask_prob, args.marker_mask_bias)
            weak, _ = model.random_mask(x, max(0.05, args.mask_prob * 0.5), args.marker_mask_bias)
            out = model(strong, r_batch)
            weak_out = model(weak, r_batch)

            # NEW: reliability-gated NeighborMix pseudo view — mix toward neighbors
            # in encoder-input space (gated by r_i), zero-mask, target = REAL cell.
            pseudo_out, pseudo_mask = None, None
            if args.pseudo_weight > 0 and nb_indices is not None:
                x_prime = neighbor_mix_view(encoder_data, idx_np, x, nb_indices, nb_probs, r_batch, args.alpha_min, args.mix_neighbors, rng).detach()
                x_prime_masked, pseudo_mask = model.random_mask(x_prime, args.mask_prob, args.marker_mask_bias)
                pseudo_out = model(x_prime_masked, r_batch)

            p_batch = None if p_targets is None else torch.as_tensor(p_targets[idx_np], dtype=torch.float32, device=device)
            loss, parts = criterion(out, weak_out, pseudo_out, target, mask, pseudo_mask, r_batch, args.pseudo_weight, p_batch, cluster_scale, recon_weight)
            if not torch.isfinite(loss):
                raise FloatingPointError(f"Non-finite loss at epoch {epoch}: {parts}")
            optimizer.zero_grad()
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 5.0)
            optimizer.step()
            for k in sums:
                sums[k] += parts[k]
            n_batches += 1
        for k in sums:
            history[k].append(sums[k] / max(1, n_batches))
        history["reliability_mean"].append(float(rel_diag["reliability_mean"]))
        history["cluster_scale"].append(cluster_scale)
        if epoch == 1 or epoch == args.epochs or epoch % 10 == 0:
            print(f"Epoch {epoch:03d}/{args.epochs} loss={history['loss'][-1]:.4f} scmae={history['scmae_loss'][-1]:.4f} pseudo={history['pseudo_loss'][-1]:.4f} cluster={history['cluster_loss'][-1]:.4f} cons={history['consistency_loss'][-1]:.4f} r_mean={history['reliability_mean'][-1]:.3f} conf={history['confidence_fraction'][-1]:.3f} scale={cluster_scale:.2f}")

    # ===EVAL===
    # First pass with full encoder to get an embedding, from which reliability is
    # computed. If --eval_gate, re-extract with the gate applied for the final embedding.
    embedding, q_out, labels_out = extract_all(model, full_loader, device)
    reliability, rel_diag = compute_reliability(embedding, q_out, nb_indices, k=args.neighbor_k)
    if args.eval_gate:
        embedding, q_out, labels_out = extract_all(model, full_loader, device, reliability=reliability)
        reliability, rel_diag = compute_reliability(embedding, q_out, nb_indices, k=args.neighbor_k)
    np.save(save_dir / "embedding_final.npy", embedding.astype(np.float32))
    np.save(save_dir / "reliability.npy", reliability.astype(np.float32))
    save_json(history, str(save_dir / "training_history.json"))

    eval_result, preds = None, None
    if not args.skip_eval:
        eval_result = family.write_kmeans_known_k_outputs(save_dir, dataset_name, DISPLAY_NAME, args.seed, embedding, labels_out, n_clusters, {"variant": METHOD_NAME, "stage": stage})
        preds = eval_result["preds"]["kmeans_known_k"]
        save_json(eval_result["fixed"], str(save_dir / "metrics.json"))

    # mixed_cell_fraction = fraction of cells that actually received neighbor mixing (r>0 and pseudo on)
    mixed_fraction = float((reliability > 1e-3).mean()) if args.pseudo_weight > 0 and nb_indices is not None else 0.0
    diag = diagnostics(embedding, labels_out, q_out, n_clusters, args.seed, preds, rel_diag, mixed_fraction)
    save_json(diag, str(save_dir / "diagnostics.json"))

    baseline = BASELINES.get(dataset_name, {})
    axis_gate = float(model.encoder.axis_gate.detach().cpu()) if hasattr(model.encoder, "axis_gate") else None
    if axis_gate is not None:
        diag["axis_gate"] = axis_gate
        save_json(diag, str(save_dir / "diagnostics.json"))
    nmi, ari, acc, f1 = metric(eval_result, "nmi"), metric(eval_result, "ari"), metric(eval_result, "acc"), metric(eval_result, "f1_macro")
    meets = bool((nmi is not None and nmi >= baseline.get("nmi", np.inf)) or (ari is not None and ari >= baseline.get("ari", np.inf)))
    summary = {
        "dataset": dataset_name, "method": DISPLAY_NAME, "method_dir": METHOD_NAME, "stage": stage,
        "seed": int(args.seed), "n_cells": int(encoder_data.shape[0]), "n_genes": int(encoder_data.shape[1]),
        "n_clusters": int(n_clusters), "runtime_seconds": float(time.time() - start),
        "fixed_metrics": eval_result["fixed"] if eval_result is not None else {},
        "diagnostics": diag, "baseline": baseline, "meets_screen_baseline_any": meets,
        "note": "Independent candidate; NOT appended to 全benchmark结果.csv.",
    }
    save_json(summary, str(save_dir / "summary.json"))
    print(f"[RESULT] {dataset_name} NMI={nmi} ARI={ari} meets_baseline={meets} collapse={diag['collapse_warning']} r_mean={rel_diag['reliability_mean']:.3f} core_frac={rel_diag['core_fraction']:.3f}")
    print(f"Completed {METHOD_NAME}. Results saved to: {save_dir}")
    return 0



if __name__ == "__main__":
    raise SystemExit(main())


