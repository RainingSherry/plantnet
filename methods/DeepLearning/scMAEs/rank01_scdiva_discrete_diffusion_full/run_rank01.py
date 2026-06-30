from __future__ import annotations

import argparse
import os
import subprocess
import sys
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")
os.environ.setdefault("NUMBA_CACHE_DIR", "/tmp/numba-cache")
os.makedirs(os.environ["MPLCONFIGDIR"], exist_ok=True)
os.makedirs(os.environ["NUMBA_CACHE_DIR"], exist_ok=True)

import numpy as np
import torch
import torch.nn.functional as F
from torch.utils.data import DataLoader

CURRENT_DIR = Path(__file__).resolve().parent
ROOT = CURRENT_DIR.parents[3]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from methods.DeepLearning import scMAE_family as family
from methods.shared_utils import ensure_dir, save_json
from model import ScDiVaDiscreteDiffusionMAE
from loss import q_sample_discrete_diffusion, scdiva_loss


def parse_args():
    p = argparse.ArgumentParser("rank01 ScDiVa discrete diffusion")
    p.add_argument("--data_path", required=True)
    p.add_argument("--save_dir", required=True)
    p.add_argument("--dataset_name", default=None)
    p.add_argument("--label_key", default="auto")
    p.add_argument("--input_mode", default="auto", choices=["auto", "raw", "log1p"])
    p.add_argument("--n_top_genes", type=int, default=1000)
    p.add_argument("--target_sum", type=float, default=10000.0)
    p.add_argument("--scale_input", type=family.str2bool, default=True)
    p.add_argument("--n_clusters", type=int, default=0)
    p.add_argument("--hidden_size", type=int, default=128)
    p.add_argument("--depth", type=int, default=4)
    p.add_argument("--heads", type=int, default=4)
    p.add_argument("--n_bins", type=int, default=32)
    p.add_argument("--diffusion_steps", type=int, default=100)
    p.add_argument("--max_mask_ratio", type=float, default=0.75)
    p.add_argument("--epochs", type=int, default=80)
    p.add_argument("--batch_size", type=int, default=128)
    p.add_argument("--lr", type=float, default=1e-3)
    p.add_argument("--weight_decay", type=float, default=1e-4)
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--gpu", type=int, default=1)
    p.add_argument("--no_cuda", action="store_true")
    p.add_argument("--skip_eval", type=family.str2bool, default=False)
    return p.parse_args()


def bins_and_tokens(x, n_bins):
    q = np.linspace(0, 100, n_bins + 1)[1:-1]
    bins = np.unique(np.percentile(x.reshape(-1).astype(np.float32), q).astype(np.float32))
    if bins.size < n_bins - 1:
        bins = np.linspace(float(x.min()), float(x.max()) + 1e-6, n_bins + 1, dtype=np.float32)[1:-1]
    return bins, np.digitize(x.astype(np.float32), bins).astype(np.int64)


def resolve_input_path(data_path: str, dataset_name: str) -> str:
    path = Path(data_path).resolve()
    if path.suffix.lower() != ".h5":
        return str(path)
    out_dir = CURRENT_DIR.parent / "benchmark_data"
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


def main():
    a = parse_args()
    family.set_seed(a.seed)
    out = Path(ensure_dir(a.save_dir))
    save_json(vars(a), str(out / "args.json"))
    (out / "source_manifest.json").write_text(
        (CURRENT_DIR / "source_manifest.json").read_text(encoding="utf-8"),
        encoding="utf-8",
    )
    device = family.get_device(a.gpu, a.no_cuda)
    ds = a.dataset_name or Path(a.data_path).stem
    resolved_data_path = resolve_input_path(a.data_path, ds)
    bundle = family.load_scmae_dataset(resolved_data_path, a.input_mode, a.n_top_genes, a.target_sum, a.scale_input, a.label_key, a.seed)
    x = bundle.data.astype(np.float32, copy=False)
    y = bundle.labels.astype(np.int64)
    n_clusters = int(a.n_clusters) if a.n_clusters > 0 else int(len(np.unique(y)))
    save_json(bundle.profile, str(out / "dataset_profile.json"))
    save_json(bundle.preprocess_config, str(out / "preprocess_config.json"))
    bins, tokens = bins_and_tokens(x, a.n_bins)
    np.save(out / "expression_token_bins.npy", bins)
    token_cpu = torch.as_tensor(tokens, dtype=torch.long)
    loader = DataLoader(family.IndexedExpressionDataset(x, y), batch_size=a.batch_size, shuffle=True, drop_last=True, generator=torch.Generator().manual_seed(a.seed))
    model = ScDiVaDiscreteDiffusionMAE(x.shape[1], a.n_bins, a.hidden_size, a.depth, a.heads, 0.1, a.diffusion_steps).to(device)
    opt = torch.optim.AdamW(model.parameters(), lr=a.lr, weight_decay=a.weight_decay)
    hist = {"loss": [], "mask_rate": []}
    for ep in range(1, max(1, a.epochs) + 1):
        model.train()
        total = 0.0; mr = 0.0; nb = 0
        for idx, xb, _ in loader:
            clean = token_cpu[idx].to(device); xb = xb.to(device)
            t = torch.randint(1, a.diffusion_steps + 1, (clean.shape[0],), device=device)
            corrupt, mask = q_sample_discrete_diffusion(clean, t, a.diffusion_steps, a.max_mask_ratio, model.mask_token_id)
            _, logits, value_pred, mask_logits = model(corrupt, t)
            parts = scdiva_loss(logits, clean, value_pred, xb, mask_logits, mask, a.n_bins)
            opt.zero_grad(set_to_none=True); parts.total.backward(); torch.nn.utils.clip_grad_norm_(model.parameters(), 5.0); opt.step()
            total += float(parts.total.detach().cpu()); mr += float(parts.mask_rate.detach().cpu()); nb += 1
        hist["loss"].append(total / max(1, nb)); hist["mask_rate"].append(mr / max(1, nb))
        if ep == 1 or ep == a.epochs or ep % 10 == 0:
            print(f"rank01 epoch={ep:03d}/{a.epochs} loss={hist['loss'][-1]:.4f} mask={hist['mask_rate'][-1]:.4f}", flush=True)
    model.eval(); embs = []
    bs = max(512, a.batch_size * 4)
    with torch.no_grad():
        for s in range(0, tokens.shape[0], bs):
            embs.append(model.feature(torch.as_tensor(tokens[s:s+bs], dtype=torch.long, device=device)).cpu().numpy())
    emb = np.concatenate(embs, axis=0).astype(np.float32)
    np.save(out / "embedding_final.npy", emb); np.save(out / "labels.npy", y); np.save(out / "gene_names.npy", bundle.gene_names.astype(str))
    family.save_embedding_h5(out / "embedding.h5", emb, y)
    save_json(hist, str(out / "training_history.json"))
    torch.save({"model": model.state_dict(), "args": vars(a)}, out / "model_checkpoint.pth")
    result = None
    if not a.skip_eval:
        result = family.write_kmeans_known_k_outputs(out, ds, "rank01_scdiva_discrete_diffusion_full", a.seed, emb, y, n_clusters, {"rank": 1, "source_paper": "ScDiVa"})
        save_json(result["fixed"], str(out / "metrics.json"))
    save_json({"dataset": ds, "rank": 1, "source_manifest": str((CURRENT_DIR / "source_manifest.json").resolve()), "fixed_metrics": result["fixed"] if result else {}}, str(out / "summary.json"))


if __name__ == "__main__":
    main()
