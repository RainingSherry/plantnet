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
    "Wang": DatasetSpec("Wang", "data/processed/Wang.h5ad", 2),
    "Macosko": DatasetSpec("Macosko", "data/processed_scmae/Macosko.h5ad", 12),
    "SRP182008": DatasetSpec("SRP182008", "data/SRP182008.h5ad", 15),
    "Melanoma_5K": DatasetSpec("Melanoma_5K", "data/processed_scmae/Melanoma_5K.h5ad", 9),
    "Quake_Smart-seq2_Lung": DatasetSpec("Quake_Smart-seq2_Lung", "data/processed/Quake_Smart-seq2_Lung.h5ad", 11),
    "Tosches": DatasetSpec("Tosches", "data/processed_scmae/Tosches.h5ad", 15),
    "worm_neuron_cell": DatasetSpec("worm_neuron_cell", "data/processed_scmae/worm_neuron_cell.h5ad", 10),
    "Pollen": DatasetSpec("Pollen", "data/processed/Pollen.h5ad", 11),
}


def parse_args():
    parser = argparse.ArgumentParser(description="Run RA/RG NeighborMix BDD phases.")
    parser.add_argument("--phase", choices=["phase0", "phase1"], default="phase0")
    parser.add_argument("--datasets", nargs="+", default=None)
    parser.add_argument("--seeds", nargs="+", type=int, default=[42])
    parser.add_argument("--gpus", nargs="+", type=int, default=[2, 3, 4, 5, 6])
    parser.add_argument("--epochs", type=int, default=80)
    parser.add_argument("--batch_size", type=int, default=256)
    parser.add_argument("--n_top_genes", type=int, default=1000)
    parser.add_argument("--out_dir", default="results/formal/ra_rg_phase0")
    parser.add_argument("--methods", nargs="+", default=None,
                        choices=[
                            "ra",
                            "ra_nomix",
                            "ra_cell",
                            "ra_intersection",
                            "ra_rec_only",
                            "ra_gene_gate",
                            "rg_none",
                            "rg_fixed",
                            "rg_reliability",
                            "rg_mutual",
                            "rg_random",
                            "rg_far",
                        ])
    parser.add_argument("--dry_run", action="store_true")
    parser.add_argument("--fail_fast", action="store_true")
    return parser.parse_args()


def validate_gpus(gpus: list[int]) -> None:
    bad = sorted(set(gpus).intersection(FORBIDDEN_GPUS))
    if bad:
        raise SystemExit(f"Forbidden GPU(s) requested: {bad}. Use only 1-6.")


def default_datasets(phase: str) -> list[str]:
    if phase == "phase0":
        return ["Wang", "Macosko"]
    return [
        "SRP182008",
        "Melanoma_5K",
        "Quake_Smart-seq2_Lung",
        "Macosko",
        "Tosches",
        "worm_neuron_cell",
        "Wang",
        "Pollen",
    ]


def default_methods(phase: str) -> list[str]:
    if phase == "phase0":
        return ["ra", "rg_none", "rg_fixed", "rg_reliability"]
    return ["rg_none", "rg_fixed", "rg_mutual", "rg_reliability", "rg_random", "rg_far"]


def command_for(method: str, ds: DatasetSpec, seed: int, save_dir: Path, args) -> list[str]:
    base = [
        sys.executable,
        str(ROOT / "methods" / "DeepLearning"),
    ]
    common = [
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
    ]
    if method == "ra":
        return [
            sys.executable,
            str(ROOT / "methods/DeepLearning/RA_NeighborMix_scMAE/run.py"),
            *common,
            "--method_name", "RA_NeighborMix_scMAE",
            "--variant_name", "ra_nm_v1",
            "--warmup_epochs", str(max(0, min(20, args.epochs // 4))),
            "--neighbor_mode", "intersection",
            "--neighbor_k", "10",
            "--mix_neighbors", "4",
            "--adaptive_mix", "true",
            "--mix_strength_max", "0.10",
            "--pseudo_weight", "0.30",
            "--pseudo_mask_loss_weight", "auto",
            "--gene_gate", "true",
            "--gene_gate_strength", "0.50",
            "--gene_adaptive_mask", "true",
        ]
    if method == "ra_gene_gate":
        return [
            sys.executable,
            str(ROOT / "methods/DeepLearning/RA_NeighborMix_scMAE/run.py"),
            *common,
            "--method_name", "RA_NeighborMix_scMAE",
            "--variant_name", "ra_nm_gene_gate_adaptive_mask",
            "--warmup_epochs", str(max(0, min(20, args.epochs // 4))),
            "--neighbor_mode", "intersection",
            "--neighbor_k", "10",
            "--mix_neighbors", "4",
            "--adaptive_mix", "true",
            "--mix_strength_max", "0.10",
            "--pseudo_weight", "0.30",
            "--pseudo_mask_loss_weight", "0.0",
            "--gene_gate", "true",
            "--gene_gate_strength", "0.50",
            "--gene_adaptive_mask", "true",
        ]
    if method == "ra_nomix":
        return [
            sys.executable,
            str(ROOT / "methods/DeepLearning/RA_NeighborMix_scMAE/run.py"),
            *common,
            "--method_name", "RA_NeighborMix_scMAE",
            "--variant_name", "ra_nm_nomix",
            "--warmup_epochs", "0",
            "--neighbor_mode", "none",
            "--use_pseudo", "false",
            "--pseudo_weight", "0",
            "--gene_gate", "false",
        ]
    if method == "ra_cell":
        return [
            sys.executable,
            str(ROOT / "methods/DeepLearning/RA_NeighborMix_scMAE/run.py"),
            *common,
            "--method_name", "RA_NeighborMix_scMAE",
            "--variant_name", "ra_nm_cell_reliability",
            "--warmup_epochs", str(max(0, min(20, args.epochs // 4))),
            "--neighbor_mode", "pca",
            "--neighbor_k", "10",
            "--mix_neighbors", "4",
            "--adaptive_mix", "true",
            "--mix_strength_max", "0.10",
            "--pseudo_weight", "0.30",
            "--pseudo_mask_loss_weight", "1.0",
            "--gene_gate", "false",
            "--gene_adaptive_mask", "false",
        ]
    if method == "ra_intersection":
        return [
            sys.executable,
            str(ROOT / "methods/DeepLearning/RA_NeighborMix_scMAE/run.py"),
            *common,
            "--method_name", "RA_NeighborMix_scMAE",
            "--variant_name", "ra_nm_reliability_intersection_knn",
            "--warmup_epochs", str(max(0, min(20, args.epochs // 4))),
            "--neighbor_mode", "intersection",
            "--neighbor_k", "10",
            "--mix_neighbors", "4",
            "--adaptive_mix", "true",
            "--mix_strength_max", "0.10",
            "--pseudo_weight", "0.30",
            "--pseudo_mask_loss_weight", "1.0",
            "--gene_gate", "false",
            "--gene_adaptive_mask", "false",
        ]
    if method == "ra_rec_only":
        return [
            sys.executable,
            str(ROOT / "methods/DeepLearning/RA_NeighborMix_scMAE/run.py"),
            *common,
            "--method_name", "RA_NeighborMix_scMAE",
            "--variant_name", "ra_nm_pseudo_rec_only",
            "--warmup_epochs", str(max(0, min(20, args.epochs // 4))),
            "--neighbor_mode", "intersection",
            "--neighbor_k", "10",
            "--mix_neighbors", "4",
            "--adaptive_mix", "true",
            "--mix_strength_max", "0.10",
            "--pseudo_weight", "0.30",
            "--pseudo_mask_loss_weight", "0.0",
            "--gene_gate", "false",
            "--gene_adaptive_mask", "false",
        ]

    mode = method.replace("rg_", "")
    gate_mode = "topology"
    edge_mode = "sim_mutual_snn_distance"
    pseudo_weight = "0.30"
    gate_max = "0.15"
    if mode == "none":
        gate_mode = "none"
        edge_mode = "none"
        pseudo_weight = "0.0"
    elif mode in {"fixed", "random", "far"}:
        gate_mode = "constant"
        edge_mode = "none"
        gate_max = "0.10"
    return [
        sys.executable,
        str(ROOT / "methods/DeepLearning/RG_NeighborMix_scMAE/run.py"),
        *common,
        "--method_name", "RG_NeighborMix_scMAE",
        "--variant_name", f"rg_nm_v1_{mode}",
        "--mix_mode", mode,
        "--gate_mode", gate_mode,
        "--edge_reliability_mode", edge_mode,
        "--neighbor_k", "10",
        "--mix_neighbors", "4",
        "--gate_max", gate_max,
        "--pseudo_weight", pseudo_weight,
    ]


def run_one(cmd_template: list[str], gpu: int, log_path: Path, dry_run: bool) -> int:
    cmd = [str(gpu) if token == "{GPU}" else token for token in cmd_template]
    if dry_run:
        print(" ".join(cmd))
        return 0
    env = os.environ.copy()
    env["CUDA_VISIBLE_DEVICES"] = str(gpu)
    env["OMP_NUM_THREADS"] = "4"
    env["OPENBLAS_NUM_THREADS"] = "4"
    env["MKL_NUM_THREADS"] = "4"
    env["NUMEXPR_NUM_THREADS"] = "4"
    env.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")
    env.setdefault("NUMBA_CACHE_DIR", "/tmp/numba-cache")
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with open(log_path, "w", encoding="utf-8") as handle:
        handle.write("Command: " + " ".join(cmd) + "\n\n")
        proc = subprocess.run(cmd, cwd=str(ROOT), env=env, text=True, stdout=handle, stderr=subprocess.STDOUT)
    return int(proc.returncode)


def main() -> int:
    args = parse_args()
    validate_gpus(args.gpus)
    datasets = args.datasets or default_datasets(args.phase)
    methods = args.methods or default_methods(args.phase)
    out_root = ROOT / args.out_dir
    jobs = []
    for ds_name in datasets:
        if ds_name not in DATASETS:
            raise SystemExit(f"Unknown dataset {ds_name!r}. Known: {sorted(DATASETS)}")
        ds = DATASETS[ds_name]
        if not (ROOT / ds.path).exists():
            raise SystemExit(f"Missing data file for {ds.name}: {ROOT / ds.path}")
        for seed in args.seeds:
            for method in methods:
                save_dir = out_root / ds.name / method / f"seed{seed}"
                log_path = save_dir / "run.log"
                jobs.append((method, ds, seed, command_for(method, ds, seed, save_dir, args), log_path))

    if args.dry_run:
        for _, _, _, cmd, _ in jobs:
            run_one(cmd, args.gpus[0], Path("/tmp/dry-run.log"), True)
        return 0

    running: list[tuple[subprocess.Popen, str, Path]] = []
    failures = 0
    next_job = 0

    def launch(job_index: int, gpu: int) -> subprocess.Popen:
        method, ds, seed, cmd_t, log_path = jobs[job_index]
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
        proc._plantnet_label = f"{ds.name}/{method}/seed{seed}/gpu{gpu}"  # type: ignore[attr-defined]
        print(f"[start] {proc._plantnet_label} pid={proc.pid}", flush=True)
        return proc

    free_gpus = list(args.gpus)
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
