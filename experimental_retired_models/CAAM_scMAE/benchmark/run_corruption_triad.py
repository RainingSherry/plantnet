#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os
import subprocess
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[4]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))


DATASETS = {
    "Quake_Smart-seq2_Lung": Path("processed/Quake_Smart-seq2_Lung.h5ad"),
    "Mouse_Pancreas_1": Path("其他/Mouse_Pancreas_1.h5ad"),
    "Limb_Muscle": Path("processed_scmae/Limb_Muscle.h5ad"),
}
CORRUPTION_TYPES = ("scmae_shuffle", "matched_donor", "nonzero_aware_donor")
SEEDS = (42, 2024, 3407)


def parse_csv(value: str) -> list[str]:
    return [part.strip() for part in value.split(",") if part.strip()]


def resolve_dataset_path(dataset_name: str, data_root: Path | None) -> Path:
    if dataset_name not in DATASETS:
        raise ValueError(f"Unknown Phase 13 dataset {dataset_name!r}")
    root = data_root if data_root is not None else PROJECT_ROOT / "data"
    path = root / DATASETS[dataset_name]
    if not path.exists():
        raise FileNotFoundError(
            f"Dataset {dataset_name!r} not found at {path}. "
            "Pass --data_root explicitly when the repo data/ symlink is not present."
        )
    return path


def main() -> int:
    parser = argparse.ArgumentParser(description="Run Phase 13 CAAM corruption triad.")
    parser.add_argument("--save_root", type=Path, default=PROJECT_ROOT / "results/CAAM_scMAE_correction/corruption_triad")
    parser.add_argument("--data_root", type=Path, default=None)
    parser.add_argument("--datasets", type=str, default=",".join(DATASETS))
    parser.add_argument("--corruption_types", type=str, default=",".join(CORRUPTION_TYPES))
    parser.add_argument("--seeds", type=str, default=",".join(str(seed) for seed in SEEDS))
    parser.add_argument("--epochs", type=int, default=10)
    parser.add_argument("--run_label", choices=("smoke", "formal"), default="formal")
    parser.add_argument("--gpu", type=int, default=1)
    parser.add_argument("--no_cuda", action="store_true")
    args = parser.parse_args()

    if not args.no_cuda and int(args.gpu) not in {1, 2, 3, 4, 5, 6}:
        raise ValueError("Phase 13 may only use physical GPU 1-6.")

    dataset_names = parse_csv(args.datasets)
    corruption_types = parse_csv(args.corruption_types)
    seeds = [int(seed) for seed in parse_csv(args.seeds)]
    unknown_corruptions = sorted(set(corruption_types) - set(CORRUPTION_TYPES))
    if unknown_corruptions:
        raise ValueError(f"Unknown Phase 13 corruption type(s): {unknown_corruptions}")

    runner = PROJECT_ROOT / "methods/DeepLearning/CAAM_scMAE/run.py"
    run_root = args.save_root / args.run_label
    run_root.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    if not args.no_cuda:
        env["CUDA_VISIBLE_DEVICES"] = str(int(args.gpu))

    for dataset_name in dataset_names:
        data_path = resolve_dataset_path(dataset_name, args.data_root)
        for corruption_type in corruption_types:
            for seed in seeds:
                out = run_root / f"{dataset_name}__{corruption_type}__seed{seed}__epochs{int(args.epochs)}"
                cmd = [
                    sys.executable,
                    str(runner),
                    "--data_path",
                    str(data_path),
                    "--save_dir",
                    str(out),
                    "--dataset_name",
                    dataset_name,
                    "--method_name",
                    f"caam_scmae_phase13_{corruption_type}",
                    "--variant",
                    "control",
                    "--n_clusters",
                    "0",
                    "--seed",
                    str(seed),
                    "--epochs",
                    str(int(args.epochs)),
                    "--benchmark_mode",
                    "true",
                    "--input_mode",
                    "log1p",
                    "--n_top_genes",
                    "2000",
                    "--scale_input",
                    "false",
                    "--corruption_type",
                    corruption_type,
                    "--strict_effective_budget",
                    "false",
                ]
                if args.no_cuda:
                    cmd.append("--no_cuda")
                else:
                    cmd.extend(["--gpu", str(int(args.gpu))])
                print("RUN", dataset_name, corruption_type, seed, f"epochs={int(args.epochs)}", flush=True)
                proc = subprocess.run(cmd, cwd=str(PROJECT_ROOT), env=env)
                if proc.returncode != 0:
                    return proc.returncode
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
