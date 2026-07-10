from __future__ import annotations

import argparse
import csv
import json
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

from .config import DATASETS, PAPER_ROOT, REMOTE_DATA_ROOT, SEEDS
from .io_utils import require_disk_space, sha256_file, utc_now, write_json
from .remote_store import RemoteStore
from .run_baseline import BASELINES


def run_task(dataset: str, data_path: Path, method: str, seed: int, gpu: int) -> dict:
    command = [
        sys.executable, "-m", "papers.scVICAR.code.run_baseline",
        "--dataset", dataset, "--data-path", str(data_path),
        "--method", method, "--seed", str(seed), "--gpu", str(gpu),
        "--upload", "--cleanup-after-upload",
    ]
    started = utc_now()
    result = subprocess.run(command, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    return {
        "dataset": dataset, "method": method, "seed": seed, "gpu": gpu,
        "started_utc": started, "finished_utc": utc_now(),
        "returncode": result.returncode,
        "status": "complete" if result.returncode == 0 else "failed",
        "message": result.stdout[-4000:],
    }


def write_status(rows: list[dict], output: Path) -> None:
    output.mkdir(parents=True, exist_ok=True)
    write_json(output / "baseline_master.json", rows)
    if rows:
        with (output / "baseline_master.csv").open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description="Run frozen external baselines for scVICAR")
    parser.add_argument("--datasets", nargs="+", default=list(DATASETS), choices=sorted(DATASETS))
    parser.add_argument("--methods", nargs="+", default=list(BASELINES), choices=sorted(BASELINES))
    parser.add_argument("--seeds", nargs="+", type=int, default=list(SEEDS), choices=SEEDS)
    parser.add_argument("--gpus", nargs="+", type=int, default=[2, 3, 4, 5, 6], choices=range(1, 7))
    parser.add_argument("--cache-dir", type=Path, default=PAPER_ROOT / ".staging/data")
    parser.add_argument("--status-dir", type=Path, default=PAPER_ROOT / "experiments/baselines_v1")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    matrix = [
        {"dataset": dataset, "method": method, "seed": seed}
        for dataset in args.datasets for method in args.methods for seed in args.seeds
    ]
    if args.dry_run:
        write_json(args.status_dir / "planned_matrix.json", matrix)
        print(json.dumps({"tasks": len(matrix), "matrix": matrix}, indent=2))
        return 0

    manifest = {
        row["name"]: row for row in json.loads(
            (PAPER_ROOT / "manifests/dataset_upload/dataset_manifest.json").read_text(encoding="utf-8")
        )
    }
    store = RemoteStore()
    store.ensure_layout()
    rows: list[dict] = []
    for dataset in args.datasets:
        require_disk_space(args.cache_dir, 5.0)
        local = args.cache_dir / f"{dataset}.h5ad"
        expected = manifest[dataset]["sha256"]
        if not local.is_file() or sha256_file(local) != expected:
            local.unlink(missing_ok=True)
            store.download_file(f"{REMOTE_DATA_ROOT}/datasets/confirmatory_v1/{dataset}.h5ad", local)
        if sha256_file(local) != expected:
            raise ValueError(f"Dataset checksum failed: {dataset}")

        gpu_tasks = [
            (method, seed) for method in args.methods if BASELINES[method].gpu for seed in args.seeds
        ]
        executors = {gpu: ThreadPoolExecutor(max_workers=1) for gpu in args.gpus}
        try:
            futures = {
                executors[args.gpus[i % len(args.gpus)]].submit(
                    run_task, dataset, local, method, seed, args.gpus[i % len(args.gpus)]
                ): (method, seed) for i, (method, seed) in enumerate(gpu_tasks)
            }
            for future in as_completed(futures):
                row = future.result()
                rows.append(row)
                write_status(rows, args.status_dir)
        finally:
            for executor in executors.values():
                executor.shutdown(wait=True, cancel_futures=True)
        if any(row["status"] == "failed" for row in rows if row["dataset"] == dataset):
            return 1

        cpu_tasks = [
            (method, seed) for method in args.methods if not BASELINES[method].gpu for seed in args.seeds
        ]
        for method, seed in cpu_tasks:
            row = run_task(dataset, local, method, seed, args.gpus[0])
            rows.append(row)
            write_status(rows, args.status_dir)
            if row["status"] == "failed":
                return 1
        local.unlink(missing_ok=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
