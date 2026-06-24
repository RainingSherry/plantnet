#!/usr/bin/env python3
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[4]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))


def infer_param_matched_hidden_dim(data_path: str) -> tuple[int, float]:
    import copy

    from methods.DeepLearning.CAAM_scMAE.data.gene_modules import build_gene_modules, normalized_assignment_dense
    from methods.DeepLearning.CAAM_scMAE.data.preprocessing import load_caam_data
    from methods.DeepLearning.CAAM_scMAE.models.caam_model import build_student
    from methods.DeepLearning.CAAM_scMAE.models.common import trainable_parameter_count
    from methods.DeepLearning.CAAM_scMAE.registry import DEFAULT_CONFIG, apply_variant

    base = copy.deepcopy(DEFAULT_CONFIG)
    base["preprocessing"]["input_mode"] = "log1p"
    base["preprocessing"]["n_top_genes"] = 0
    base["preprocessing"]["scale_input"] = False
    base["benchmark_mode"] = True
    base["seed"] = 0
    bundle = load_caam_data(
        data_path,
        input_mode="log1p",
        target_sum=float(base["preprocessing"]["target_sum"]),
        n_top_genes=0,
        scale_input=False,
        benchmark_mode=True,
        seed=0,
    )
    _ids, sparse_assignment = build_gene_modules(
        bundle.x,
        int(base["axial"]["n_gene_modules"]),
        int(base["axial"]["module_svd_dim"]),
        int(base["axial"]["module_seed"]),
        None,
    )
    assignment = normalized_assignment_dense(sparse_assignment)
    axial_cfg = apply_variant(base, "axial")
    target = trainable_parameter_count(build_student(n_genes=bundle.x.shape[1], config=axial_cfg, assignment=assignment))

    def count_for(hidden: int) -> int:
        cfg = apply_variant(copy.deepcopy(base), "control")
        cfg["model"]["mlp_hidden_dim"] = int(hidden)
        return trainable_parameter_count(build_student(n_genes=bundle.x.shape[1], config=cfg))

    low, high = 1, max(4, int(base["model"]["mlp_hidden_dim"]))
    while count_for(high) < target and high < 65536:
        low = high
        high *= 2
    for _ in range(24):
        mid = (low + high) // 2
        if mid <= low:
            break
        if count_for(mid) < target:
            low = mid
        else:
            high = mid

    best_hidden = high
    best_gap = float("inf")
    for hidden in range(max(1, low - 64), high + 65):
        count = count_for(hidden)
        gap = abs(count - target) / max(1, target)
        if gap < best_gap:
            best_gap = gap
            best_hidden = hidden
    return best_hidden, best_gap


def main() -> int:
    parser = argparse.ArgumentParser(description="Run internal CAAM-scMAE ablations.")
    parser.add_argument("--data_path", required=True)
    parser.add_argument("--save_root", default="results/CAAM_scMAE_ablation")
    parser.add_argument("--n_clusters", type=int, required=True)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--epochs", type=int, default=10)
    parser.add_argument("--gpu", type=int, default=1)
    parser.add_argument("--no_cuda", action="store_true")
    args = parser.parse_args()
    root = Path(args.save_root)
    root.mkdir(parents=True, exist_ok=True)
    runner = Path(__file__).resolve().parents[1] / "run.py"
    hidden_dim, gap = infer_param_matched_hidden_dim(args.data_path)
    if gap > 0.05:
        raise SystemExit(
            f"Could not parameter-match MLP within +/-5% of Axial student params; "
            f"best hidden_dim={hidden_dim}, relative_gap={gap:.4f}"
        )
    variants = [
        ("control", "caam_scmae_control", []),
        ("axial", "caam_scmae_axial", []),
        ("advmask", "caam_scmae_advmask", []),
        ("full", "caam_scmae_full", []),
        ("control", "caam_scmae_mlp_parammatched", ["--mlp_hidden_dim", str(hidden_dim)]),
    ]
    for variant, method_name, extra in variants:
        run_name = "mlp_parammatched" if method_name.endswith("parammatched") else variant
        out = root / f"{run_name}__seed{args.seed}"
        cmd = [
            sys.executable,
            str(runner),
            "--data_path",
            args.data_path,
            "--save_dir",
            str(out),
            "--n_clusters",
            str(args.n_clusters),
            "--seed",
            str(args.seed),
            "--epochs",
            str(args.epochs),
            "--variant",
            variant,
            "--method_name",
            method_name,
            "--benchmark_mode",
            "true",
            "--input_mode",
            "log1p",
            "--n_top_genes",
            "0",
            "--scale_input",
            "false",
        ]
        cmd.extend(extra)
        if args.no_cuda:
            cmd.append("--no_cuda")
        else:
            cmd.extend(["--gpu", str(args.gpu)])
        proc = subprocess.run(cmd)
        if proc.returncode != 0:
            return proc.returncode
    print(f"parameter_matched_mlp_hidden_dim={hidden_dim} relative_param_gap={gap:.4f}")
    print("TODO: paired bootstrap aggregation is required before reporting claims.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
