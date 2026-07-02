#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
FORBIDDEN_GPUS = {0, 7}


@dataclass(frozen=True)
class DatasetSpec:
    name: str
    path: str
    n_clusters: int


DATASETS = {
    "SRP182008": DatasetSpec("SRP182008", "data/SRP182008.h5ad", 15),
    "Melanoma_5K": DatasetSpec("Melanoma_5K", "data/processed_scmae/Melanoma_5K.h5ad", 9),
    "Macosko": DatasetSpec("Macosko", "data/processed_scmae/Macosko.h5ad", 12),
    "Tosches": DatasetSpec("Tosches", "data/processed_scmae/Tosches.h5ad", 15),
    "Wang": DatasetSpec("Wang", "data/processed/Wang.h5ad", 2),
}


def parse_args():
    parser = argparse.ArgumentParser(description="Run RG_NeighborMix_scMAE Phase 2 sensitivity sweeps.")
    parser.add_argument("--sweep", choices=["gate_max", "neighbor_k", "pseudo_weight"], required=True)
    parser.add_argument("--datasets", nargs="+", default=["SRP182008", "Melanoma_5K", "Macosko", "Tosches", "Wang"])
    parser.add_argument("--seeds", nargs="+", type=int, default=[42])
    parser.add_argument("--gpus", nargs="+", type=int, default=[1, 2, 3, 4, 5, 6])
    parser.add_argument("--epochs", type=int, default=80)
    parser.add_argument("--batch_size", type=int, default=256)
    parser.add_argument("--n_top_genes", type=int, default=1000)
    parser.add_argument("--out_dir", default="results/formal/rg_phase2_sensitivity_e80")
    parser.add_argument("--gate_max_values", nargs="+", type=float, default=[0.05, 0.10, 0.15, 0.20])
    parser.add_argument("--neighbor_k_values", nargs="+", type=int, default=[5, 10, 15])
    parser.add_argument("--pseudo_weight_values", nargs="+", type=float, default=[0.10, 0.30])
    parser.add_argument("--base_gate_max", type=float, default=0.15)
    parser.add_argument("--base_neighbor_k", type=int, default=10)
    parser.add_argument("--base_pseudo_weight", type=float, default=0.30)
    parser.add_argument("--dry_run", action="store_true")
    parser.add_argument("--fail_fast", action="store_true")
    return parser.parse_args()


def fmt_value(value: float | int) -> str:
    if isinstance(value, int):
        return str(value)
    return f"{value:.3g}".replace(".", "p")


def validate_gpus(gpus: list[int]) -> None:
    bad = sorted(set(gpus).intersection(FORBIDDEN_GPUS))
    if bad:
        raise SystemExit(f"Forbidden GPU(s) requested: {bad}. Use only 1-6.")


def values_for(args) -> list[float | int]:
    if args.sweep == "gate_max":
        return args.gate_max_values
    if args.sweep == "neighbor_k":
        return args.neighbor_k_values
    return args.pseudo_weight_values


def command_for(args, ds: DatasetSpec, seed: int, value: float | int, save_dir: Path) -> list[str]:
    gate_max = args.base_gate_max
    neighbor_k = args.base_neighbor_k
    pseudo_weight = args.base_pseudo_weight
    if args.sweep == "gate_max":
        gate_max = float(value)
    elif args.sweep == "neighbor_k":
        neighbor_k = int(value)
    elif args.sweep == "pseudo_weight":
        pseudo_weight = float(value)
    return [
        sys.executable,
        str(ROOT / "experimental_retired_models/RG_NeighborMix_scMAE/run.py"),
        "--data_path", str(ROOT / ds.path),
        "--save_dir", str(save_dir),
        "--dataset_name", ds.name,
        "--n_clusters", str(ds.n_clusters),
        "--epochs", str(args.epochs),
        "--batch_size", str(args.batch_size),
        "--n_top_genes", str(args.n_top_genes),
        "--seed", str(seed),
        "--gpu", "{GPU}",
        "--no_save_h5ad",
        "--method_name", "RG_NeighborMix_scMAE",
        "--variant_name", f"rg_nm_phase2_{args.sweep}_{fmt_value(value)}",
        "--mix_mode", "reliability",
        "--gate_mode", "topology",
        "--edge_reliability_mode", "sim_mutual_snn_distance",
        "--neighbor_k", str(neighbor_k),
        "--mix_neighbors", "4",
        "--gate_max", str(gate_max),
        "--pseudo_weight", str(pseudo_weight),
    ]


def main() -> int:
    args = parse_args()
    validate_gpus(args.gpus)
    out_root = ROOT / args.out_dir / args.sweep
    jobs = []
    for ds_name in args.datasets:
        if ds_name not in DATASETS:
            raise SystemExit(f"Unknown dataset {ds_name!r}. Known: {sorted(DATASETS)}")
        ds = DATASETS[ds_name]
        if not (ROOT / ds.path).exists():
            raise SystemExit(f"Missing data file for {ds.name}: {ROOT / ds.path}")
        for seed in args.seeds:
            for value in values_for(args):
                label_value = f"{args.sweep}_{fmt_value(value)}"
                save_dir = out_root / ds.name / label_value / f"seed{seed}"
                jobs.append((ds, seed, value, command_for(args, ds, seed, value, save_dir), save_dir / "run.log"))

    if args.dry_run:
        for _, _, _, cmd_t, _ in jobs:
            print(" ".join(str(args.gpus[0]) if token == "{GPU}" else token for token in cmd_t))
        return 0

    running: list[tuple[subprocess.Popen, str, Path]] = []
    free_gpus = list(args.gpus)
    failures = 0
    next_job = 0

    def launch(job_index: int, gpu: int) -> subprocess.Popen:
        ds, seed, value, cmd_t, log_path = jobs[job_index]
        cmd = [str(gpu) if token == "{GPU}" else token for token in cmd_t]
        env = os.environ.copy()
        env["CUDA_VISIBLE_DEVICES"] = str(gpu)
        env["OMP_NUM_THREADS"] = "4"
        env["OPENBLAS_NUM_THREADS"] = "4"
        env["MKL_NUM_THREADS"] = "4"
        env["NUMEXPR_NUM_THREADS"] = "4"
        env.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")
        env.setdefault("NUMBA_CACHE_DIR", "/tmp/numba-cache")
        log_path.parent.mkdir(parents=True, exist_ok=True)
        handle = open(log_path, "w", encoding="utf-8")
        handle.write("Command: " + " ".join(cmd) + "\n\n")
        handle.flush()
        proc = subprocess.Popen(cmd, cwd=str(ROOT), env=env, stdout=handle, stderr=subprocess.STDOUT, text=True)
        proc._plantnet_log_handle = handle  # type: ignore[attr-defined]
        proc._plantnet_label = f"{ds.name}/{args.sweep}={value}/seed{seed}/gpu{gpu}"  # type: ignore[attr-defined]
        print(f"[start] {proc._plantnet_label} pid={proc.pid}", flush=True)
        return proc

    while next_job < len(jobs) or running:
        while free_gpus and next_job < len(jobs):
            gpu = free_gpus.pop(0)
            proc = launch(next_job, gpu)
            running.append((proc, proc._plantnet_label, jobs[next_job][4]))  # type: ignore[attr-defined]
            next_job += 1
        time.sleep(5)
        still = []
        for proc, label, log_path in running:
            rc = proc.poll()
            if rc is None:
                still.append((proc, label, log_path))
                continue
            proc._plantnet_log_handle.close()  # type: ignore[attr-defined]
            gpu = int(label.rsplit("gpu", 1)[1])
            free_gpus.append(gpu)
            if rc == 0:
                print(f"[ok] {label}", flush=True)
            else:
                failures += 1
                print(f"[fail] {label} rc={rc} log={log_path}", flush=True)
                if args.fail_fast:
                    for p, _, _ in still:
                        p.terminate()
                    return rc
        running = still
        free_gpus.sort()
    print(f"Done. jobs={len(jobs)} failures={failures} out={out_root}", flush=True)
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
