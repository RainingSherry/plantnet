#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
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
SEEDS = (42, 2024, 3407)
VARIANTS = ("control", "advmask")
DEFAULT_PHASE13_REPORT = PROJECT_ROOT / "results/CAAM_scMAE_correction/corruption_triad_report.json"


def parse_csv(value: str) -> list[str]:
    return [part.strip() for part in value.split(",") if part.strip()]


def resolve_dataset_path(dataset_name: str, data_root: Path | None) -> Path:
    if dataset_name not in DATASETS:
        raise ValueError(f"Unknown Phase 14 dataset {dataset_name!r}")
    root = data_root if data_root is not None else PROJECT_ROOT / "data"
    path = root / DATASETS[dataset_name]
    if not path.exists():
        raise FileNotFoundError(
            f"Dataset {dataset_name!r} not found at {path}. "
            "Pass --data_root explicitly when the repo data/ symlink is not present."
        )
    return path


def load_recommended_corruptions(report_path: Path) -> list[str]:
    report = json.loads(report_path.read_text(encoding="utf-8"))
    recommendation = report.get("recommendation", {})
    recommended = recommendation.get("recommended_corruption_type")
    if not recommended:
        raise ValueError(f"{report_path}: missing recommended_corruption_type for Phase 14")
    if isinstance(recommended, list):
        out = [str(item) for item in recommended]
    else:
        out = [str(recommended)]
    if not 1 <= len(out) <= 2:
        raise ValueError(f"Phase 14 expects 1-2 recommended corruption types, got {out}")
    return out


def build_run_commands(
    *,
    runner: Path,
    run_root: Path,
    data_root: Path | None,
    dataset_names: list[str],
    corruption_types: list[str],
    seeds: list[int],
    epochs: int,
    gpu: int,
    no_cuda: bool,
) -> list[tuple[str, str, int, str, Path, list[str]]]:
    commands: list[tuple[str, str, int, str, Path, list[str]]] = []
    for dataset_name in dataset_names:
        data_path = resolve_dataset_path(dataset_name, data_root)
        for corruption_type in corruption_types:
            for variant in VARIANTS:
                for seed in seeds:
                    out = run_root / f"{dataset_name}__{corruption_type}__{variant}__seed{seed}__epochs{int(epochs)}"
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
                        f"caam_scmae_phase14_{variant}",
                        "--variant",
                        variant,
                        "--n_clusters",
                        "0",
                        "--seed",
                        str(seed),
                        "--epochs",
                        str(int(epochs)),
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
                    if no_cuda:
                        cmd.append("--no_cuda")
                    else:
                        cmd.extend(["--gpu", str(int(gpu))])
                    commands.append((dataset_name, corruption_type, seed, variant, out, cmd))
    return commands


def main() -> int:
    parser = argparse.ArgumentParser(description="Run Phase 14 CAAM AdvMask triage.")
    parser.add_argument("--save_root", type=Path, default=PROJECT_ROOT / "results/CAAM_scMAE_correction/advmask_triage")
    parser.add_argument("--data_root", type=Path, default=None)
    parser.add_argument("--phase13_report", type=Path, default=DEFAULT_PHASE13_REPORT)
    parser.add_argument("--datasets", type=str, default=",".join(DATASETS))
    parser.add_argument("--corruption_types", type=str, default=None)
    parser.add_argument("--seeds", type=str, default=",".join(str(seed) for seed in SEEDS))
    parser.add_argument("--epochs", type=int, default=3)
    parser.add_argument("--run_label", choices=("smoke", "formal"), default="formal")
    parser.add_argument("--gpu", type=int, default=1)
    parser.add_argument("--no_cuda", action="store_true")
    args = parser.parse_args()

    if not args.no_cuda and int(args.gpu) not in {1, 2, 3, 4, 5, 6}:
        raise ValueError("Phase 14 may only use physical GPU 1-6.")

    corruption_types = (
        parse_csv(args.corruption_types)
        if args.corruption_types is not None
        else load_recommended_corruptions(args.phase13_report)
    )
    if len(corruption_types) > 2:
        raise ValueError("Phase 14 may triage at most two corruption types.")

    runner = PROJECT_ROOT / "methods/DeepLearning/CAAM_scMAE/run.py"
    run_root = args.save_root / args.run_label
    run_root.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    if not args.no_cuda:
        env["CUDA_VISIBLE_DEVICES"] = str(int(args.gpu))

    commands = build_run_commands(
        runner=runner,
        run_root=run_root,
        data_root=args.data_root,
        dataset_names=parse_csv(args.datasets),
        corruption_types=corruption_types,
        seeds=[int(seed) for seed in parse_csv(args.seeds)],
        epochs=int(args.epochs),
        gpu=int(args.gpu),
        no_cuda=bool(args.no_cuda),
    )
    for dataset_name, corruption_type, seed, variant, _out, cmd in commands:
        print("RUN", dataset_name, corruption_type, variant, seed, f"epochs={int(args.epochs)}", flush=True)
        proc = subprocess.run(cmd, cwd=str(PROJECT_ROOT), env=env)
        if proc.returncode != 0:
            return proc.returncode
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
