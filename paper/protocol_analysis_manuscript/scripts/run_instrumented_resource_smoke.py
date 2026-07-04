#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import os
import subprocess
import sys
import time
from pathlib import Path

import psutil


SCRIPT_DIR = Path(__file__).resolve().parent
MANUSCRIPT_ROOT = SCRIPT_DIR.parent
PROJECT_ROOT = MANUSCRIPT_ROOT.parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

DATA_ROOT = Path("/data/luolie/biopipeline/dimension-reduction/plantnet/data")
DEFAULT_DATASETS = {
    "Quake_Smart-seq2_Lung": DATA_ROOT / "processed/Quake_Smart-seq2_Lung.h5ad",
    "Mouse_Pancreas_1": DATA_ROOT / "其他/Mouse_Pancreas_1.h5ad",
    "Limb_Muscle": DATA_ROOT / "processed_scmae/Limb_Muscle.h5ad",
}
DEFAULT_RUN_ROOT = Path("/tmp/caam_resource_smoke/dev_20260626")
RUNNER = PROJECT_ROOT / "methods/DeepLearning/CAAM_scMAE/run.py"


def load_json(path: Path) -> dict:
    with open(path, "r", encoding="utf-8") as handle:
        return json.load(handle)


def save_json(path: Path, data: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(data, handle, indent=2, sort_keys=True)


def parse_gnu_time(path: Path) -> dict[str, str]:
    out: dict[str, str] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if ":" not in line:
            continue
        key, value = line.split(":", 1)
        out[key.strip()] = value.strip()
    return out


def parse_csv(value: str) -> list[str]:
    return [part.strip() for part in value.split(",") if part.strip()]


def query_gpu_total_memory_mib(gpu: int) -> int:
    proc = subprocess.run(
        [
            "nvidia-smi",
            "-i",
            str(int(gpu)),
            "--query-gpu=memory.used",
            "--format=csv,noheader,nounits",
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    text = proc.stdout.strip()
    if not text:
        raise ValueError(f"nvidia-smi returned no memory value for GPU {gpu}")
    return int(text.splitlines()[0].strip())


def process_tree_pids(proc: psutil.Process) -> set[int]:
    pids = {int(proc.pid)}
    for child in proc.children(recursive=True):
        pids.add(int(child.pid))
    return pids


def process_tree_rss_mib(proc: psutil.Process) -> float:
    processes = [proc, *proc.children(recursive=True)]
    total = 0
    for item in processes:
        try:
            total += int(item.memory_info().rss)
        except psutil.NoSuchProcess:
            continue
    return total / (1024 * 1024)


def infer_param_matched_hidden_dim(data_path: Path) -> int:
    from methods.DeepLearning.CAAM_scMAE.benchmark.run_ablation import infer_param_matched_hidden_dim

    hidden_dim, gap = infer_param_matched_hidden_dim(str(data_path))
    if float(gap) > 0.05:
        raise ValueError(f"Parameter-matched MLP gap too high: hidden_dim={hidden_dim}, relative_gap={gap:.6f}")
    return int(hidden_dim)


def write_smoke_config(path: Path) -> Path:
    save_json(
        path,
        {
            "training": {
                "student_warmup_epochs": 1,
                "generator_update_interval": 5,
            }
        },
    )
    return path


def build_command(
    *,
    config_path: Path,
    data_path: Path,
    save_dir: Path,
    dataset_name: str,
    condition: str,
    variant: str,
    seed: int,
    epochs: int,
    gpu: int,
    mlp_hidden_dim: int | None,
) -> list[str]:
    cmd = [
        "/usr/bin/time",
        "-v",
        "-o",
        str(save_dir / "gnu_time.txt"),
        sys.executable,
        str(RUNNER),
        "--config",
        str(config_path),
        "--data_path",
        str(data_path),
        "--save_dir",
        str(save_dir),
        "--dataset_name",
        dataset_name,
        "--method_name",
        f"caam_scmae_resource_smoke_{condition}",
        "--variant",
        variant,
        "--n_clusters",
        "0",
        "--seed",
        str(int(seed)),
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
        "scmae_shuffle",
        "--strict_effective_budget",
        "false",
        "--gpu",
        str(int(gpu)),
    ]
    if mlp_hidden_dim is not None:
        cmd.extend(["--mlp_hidden_dim", str(int(mlp_hidden_dim))])
    return cmd


def run_instrumented(
    *,
    cmd: list[str],
    save_dir: Path,
    gpu: int,
    env: dict[str, str],
    sample_interval_seconds: float,
) -> dict:
    save_dir.mkdir(parents=True, exist_ok=True)
    baseline_gpu_memory_mib = query_gpu_total_memory_mib(gpu)
    start = time.perf_counter()
    proc = subprocess.Popen(cmd, cwd=str(PROJECT_ROOT), env=env)
    ps_proc = psutil.Process(proc.pid)
    peak_cpu_rss_mib = 0.0
    peak_gpu_total_memory_mib = baseline_gpu_memory_mib
    sample_count = 0
    while proc.poll() is None:
        peak_cpu_rss_mib = max(peak_cpu_rss_mib, process_tree_rss_mib(ps_proc))
        peak_gpu_total_memory_mib = max(peak_gpu_total_memory_mib, query_gpu_total_memory_mib(gpu))
        sample_count += 1
        time.sleep(float(sample_interval_seconds))
    return_code = proc.wait()
    wall_clock_seconds = time.perf_counter() - start
    peak_gpu_delta_mib = max(0, peak_gpu_total_memory_mib - baseline_gpu_memory_mib)
    profile = {
        "wall_clock_seconds": wall_clock_seconds,
        "peak_cpu_rss_mib_process_tree": peak_cpu_rss_mib,
        "baseline_gpu_total_memory_mib": baseline_gpu_memory_mib,
        "peak_gpu_total_memory_mib": peak_gpu_total_memory_mib,
        "peak_gpu_memory_delta_mib": peak_gpu_delta_mib,
        "peak_gpu_memory_mib_process_tree": peak_gpu_delta_mib,
        "gpu_memory_sampling_mode": "total_gpu_memory_delta",
        "sample_count": sample_count,
        "sample_interval_seconds": float(sample_interval_seconds),
        "return_code": int(return_code),
    }
    if (save_dir / "gnu_time.txt").exists():
        profile["gnu_time"] = parse_gnu_time(save_dir / "gnu_time.txt")
    save_json(save_dir / "resource_profile.json", profile)

    runtime_path = save_dir / "runtime.json"
    if runtime_path.exists():
        runtime = load_json(runtime_path)
        runtime.update(
            {
                "wall_clock_seconds": wall_clock_seconds,
                "peak_cpu_rss_mib_process_tree": peak_cpu_rss_mib,
                "baseline_gpu_total_memory_mib": baseline_gpu_memory_mib,
                "peak_gpu_total_memory_mib": peak_gpu_total_memory_mib,
                "peak_gpu_memory_delta_mib": peak_gpu_delta_mib,
                "peak_gpu_memory_mib_process_tree": peak_gpu_delta_mib,
                "gpu_memory_sampling_mode": "total_gpu_memory_delta",
                "resource_profile_path": "resource_profile.json",
            }
        )
        save_json(runtime_path, runtime)
    if return_code != 0:
        raise RuntimeError(f"Run failed with return code {return_code}: {' '.join(cmd)}")
    return profile


def flatten_run(save_dir: Path, *, condition: str, variant: str, dataset_name: str, seed: int, epochs: int) -> dict:
    run_manifest = load_json(save_dir / "run_manifest.json")
    runtime = load_json(save_dir / "runtime.json")
    metrics = load_json(save_dir / "metrics.json")
    return {
        "dataset": dataset_name,
        "condition": condition,
        "variant": variant,
        "seed": int(seed),
        "epochs": int(epochs),
        "student_trainable_params": int(run_manifest["student_trainable_params"]),
        "generator_trainable_params": int(run_manifest["generator_trainable_params"]),
        "total_trainable_params": int(run_manifest["student_trainable_params"]) + int(run_manifest["generator_trainable_params"]),
        "kmeans_known_k.ari": float(metrics["kmeans_known_k"]["ari"]),
        "leiden_fixed.ari": float(metrics["leiden_fixed"]["ari"]),
        "wall_clock_seconds": float(runtime["wall_clock_seconds"]),
        "peak_cpu_rss_mib_process_tree": float(runtime["peak_cpu_rss_mib_process_tree"]),
        "baseline_gpu_total_memory_mib": int(runtime["baseline_gpu_total_memory_mib"]),
        "peak_gpu_total_memory_mib": int(runtime["peak_gpu_total_memory_mib"]),
        "peak_gpu_memory_delta_mib": int(runtime["peak_gpu_memory_delta_mib"]),
        "gpu_memory_sampling_mode": str(runtime["gpu_memory_sampling_mode"]),
        "logical_device": str(runtime["logical_device"]),
        "physical_gpu": "" if runtime["physical_gpu"] is None else str(runtime["physical_gpu"]),
        "run_dir": str(save_dir),
    }


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        raise ValueError(f"No rows to write for {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def latex_escape(value: object) -> str:
    return str(value).replace("_", r"\_")


def write_latex_table(path: Path, rows: list[dict]) -> None:
    lines = [
        r"\begin{table}[ht]",
        r"\centering",
        r"\caption{Instrumented development resource smoke across development datasets, seed 42, three epochs. This table is development evidence only.}",
        r"\label{tab:generated-instrumented-resource-smoke}",
        r"\scriptsize",
        r"\begin{tabular}{llrrrrrr}",
        r"\toprule",
        r"Dataset & Condition & Params & Gen. params & Known-K ARI & Leiden ARI & Wall sec. & GPU delta MiB \\",
        r"\midrule",
    ]
    for row in rows:
        lines.append(
            " & ".join(
                [
                    latex_escape(row["dataset"]),
                    latex_escape(row["condition"]),
                    f"{int(row['total_trainable_params']):,}",
                    f"{int(row['generator_trainable_params']):,}",
                    f"{float(row['kmeans_known_k.ari']):.6f}",
                    f"{float(row['leiden_fixed.ari']):.6f}",
                    f"{float(row['wall_clock_seconds']):.2f}",
                    str(int(row["peak_gpu_memory_delta_mib"])),
                ]
            )
            + r" \\"
        )
    lines.extend([r"\bottomrule", r"\end{tabular}", r"\end{table}", ""])
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")


def write_latex_summary_table(path: Path, rows: list[dict]) -> None:
    grouped: dict[str, list[dict]] = {}
    for row in rows:
        grouped.setdefault(str(row["condition"]), []).append(row)
    lines = [
        r"\begin{table}[ht]",
        r"\centering",
        r"\caption{Condition-level mean of the instrumented development resource smoke. Means are over development datasets only.}",
        r"\label{tab:generated-instrumented-resource-smoke-summary}",
        r"\scriptsize",
        r"\begin{tabular}{lrrrrrr}",
        r"\toprule",
        r"Condition & Datasets & Params & Known-K ARI & Leiden ARI & Wall sec. & GPU delta MiB \\",
        r"\midrule",
    ]
    condition_order = ["control", "advmask", "axial", "mlp_parammatched"]
    for condition in condition_order:
        group_rows = grouped[condition]
        lines.append(
            " & ".join(
                [
                    latex_escape(condition),
                    str(len(group_rows)),
                    f"{int(round(sum(float(row['total_trainable_params']) for row in group_rows) / len(group_rows))):,}",
                    f"{sum(float(row['kmeans_known_k.ari']) for row in group_rows) / len(group_rows):.6f}",
                    f"{sum(float(row['leiden_fixed.ari']) for row in group_rows) / len(group_rows):.6f}",
                    f"{sum(float(row['wall_clock_seconds']) for row in group_rows) / len(group_rows):.2f}",
                    f"{sum(float(row['peak_gpu_memory_delta_mib']) for row in group_rows) / len(group_rows):.1f}",
                ]
            )
            + r" \\"
        )
    lines.extend([r"\bottomrule", r"\end{tabular}", r"\end{table}", ""])
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")


def resolve_dataset_specs(args: argparse.Namespace) -> list[tuple[str, Path]]:
    if args.data_path is not None or args.dataset_name is not None:
        if args.data_path is None or args.dataset_name is None:
            raise ValueError("--data-path and --dataset-name must be provided together for a custom single-dataset run.")
        return [(str(args.dataset_name), Path(args.data_path))]
    specs = []
    for name in parse_csv(args.datasets):
        if name not in DEFAULT_DATASETS:
            raise ValueError(f"Unknown development dataset {name!r}; expected one of {sorted(DEFAULT_DATASETS)}")
        specs.append((name, DEFAULT_DATASETS[name]))
    return specs


def main() -> int:
    parser = argparse.ArgumentParser(description="Run a small instrumented CAAM resource smoke on development data.")
    parser.add_argument("--datasets", default=",".join(DEFAULT_DATASETS))
    parser.add_argument("--data-path", type=Path, default=None)
    parser.add_argument("--dataset-name", default=None)
    parser.add_argument("--run-root", type=Path, default=DEFAULT_RUN_ROOT)
    parser.add_argument("--output-dir", type=Path, default=MANUSCRIPT_ROOT / "generated/instrumented_resource_smoke")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--epochs", type=int, default=3)
    parser.add_argument("--gpu", type=int, default=1)
    parser.add_argument("--sample-interval-seconds", type=float, default=0.5)
    parser.add_argument("--summarize-only", action="store_true")
    args = parser.parse_args()

    dataset_specs = resolve_dataset_specs(args)
    for dataset_name, data_path in dataset_specs:
        if not data_path.exists():
            raise FileNotFoundError(f"Development dataset {dataset_name!r} not found: {data_path}")
    if int(args.gpu) not in {1, 2, 3, 4, 5, 6}:
        raise ValueError("Use physical GPU 1-6 for CAAM development runs.")

    args.run_root.mkdir(parents=True, exist_ok=True)
    config_path = write_smoke_config(args.run_root / "resource_smoke_config.json")
    env = os.environ.copy()
    env["CUDA_VISIBLE_DEVICES"] = str(int(args.gpu))
    rows = []
    hidden_dims = {}
    for dataset_name, data_path in dataset_specs:
        hidden_dim = infer_param_matched_hidden_dim(data_path)
        hidden_dims[dataset_name] = hidden_dim
        specs = [
            ("control", "control", None),
            ("advmask", "advmask", None),
            ("axial", "axial", None),
            ("mlp_parammatched", "control", hidden_dim),
        ]
        for condition, variant, hidden in specs:
            save_dir = args.run_root / f"{dataset_name}__scmae_shuffle__{condition}__seed{int(args.seed)}__epochs{int(args.epochs)}"
            cmd = build_command(
                config_path=config_path,
                data_path=data_path,
                save_dir=save_dir,
                dataset_name=dataset_name,
                condition=condition,
                variant=variant,
                seed=int(args.seed),
                epochs=int(args.epochs),
                gpu=int(args.gpu),
                mlp_hidden_dim=hidden,
            )
            if not args.summarize_only:
                run_instrumented(
                    cmd=cmd,
                    save_dir=save_dir,
                    gpu=int(args.gpu),
                    env=env,
                    sample_interval_seconds=float(args.sample_interval_seconds),
                )
            rows.append(
                flatten_run(
                    save_dir,
                    condition=condition,
                    variant=variant,
                    dataset_name=dataset_name,
                    seed=int(args.seed),
                    epochs=int(args.epochs),
                )
            )

    args.output_dir.mkdir(parents=True, exist_ok=True)
    write_csv(args.output_dir / "instrumented_resource_smoke.csv", rows)
    write_latex_table(args.output_dir / "instrumented_resource_smoke.tex", rows)
    write_latex_summary_table(args.output_dir / "instrumented_resource_smoke_summary.tex", rows)
    save_json(
        args.output_dir / "artifact_manifest.json",
        {
            "run_root": str(args.run_root),
            "datasets": [{"name": name, "data_path": str(path)} for name, path in dataset_specs],
            "seed": int(args.seed),
            "epochs": int(args.epochs),
            "conditions": sorted({row["condition"] for row in rows}),
            "claim_scope": "development-only resource smoke; no validation or sealed test data",
            "parameter_matched_mlp_hidden_dims": hidden_dims,
            "gpu_memory_sampling_mode": "total_gpu_memory_delta",
            "source_runner": str(RUNNER),
        },
    )
    print(f"Wrote instrumented resource smoke to {args.output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
