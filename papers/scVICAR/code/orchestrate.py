from __future__ import annotations

import argparse
import csv
import json
import os
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

from .config import DATASETS, PAPER_ROOT, REMOTE_DATA_ROOT, SEEDS, VARIANTS
from .io_utils import require_disk_space, sha256_file, utc_now, write_json
from .remote_store import RemoteStore


def load_manifest(path: Path) -> dict[str, dict]:
    rows = json.loads(path.read_text(encoding="utf-8"))
    return {str(row["name"]): row for row in rows}


def run_task(dataset: str, data_path: Path, variant: str, seed: int, gpu: int, staging_root: Path, epochs: int | None) -> dict:
    command = [
        sys.executable, "-m", "papers.scVICAR.code.run_scVICAR",
        "--dataset", dataset,
        "--data-path", str(data_path),
        "--variant", variant,
        "--seed", str(seed),
        "--gpu", str(gpu),
        "--staging-root", str(staging_root),
        "--upload", "--cleanup-after-upload",
    ]
    if epochs is not None:
        command.extend(["--epochs", str(epochs)])
    started = utc_now()
    result = subprocess.run(command, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    return {
        "dataset": dataset,
        "variant": variant,
        "seed": seed,
        "gpu": gpu,
        "started_utc": started,
        "finished_utc": utc_now(),
        "returncode": result.returncode,
        "status": "complete" if result.returncode == 0 else "failed",
        "message": result.stdout[-4000:],
    }


def write_status(rows: list[dict], output_dir: Path) -> None:
    write_json(output_dir / "run_master.json", rows)
    if not rows:
        return
    with (output_dir / "run_master.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description="Run the preregistered scVICAR confirmatory matrix")
    parser.add_argument("--datasets", nargs="+", default=list(DATASETS), choices=sorted(DATASETS))
    parser.add_argument("--variants", nargs="+", default=list(VARIANTS), choices=sorted(VARIANTS))
    parser.add_argument("--seeds", nargs="+", type=int, default=list(SEEDS))
    parser.add_argument("--gpus", nargs="+", type=int, default=[1, 2, 3, 4, 5, 6], choices=range(1, 7))
    parser.add_argument("--epochs", type=int)
    parser.add_argument("--manifest", type=Path, default=PAPER_ROOT / "manifests" / "dataset_upload" / "dataset_manifest.json")
    parser.add_argument("--cache-dir", type=Path, default=PAPER_ROOT / ".staging" / "data")
    parser.add_argument("--staging-root", type=Path, default=PAPER_ROOT / ".staging" / "runs")
    parser.add_argument("--status-dir", type=Path, default=PAPER_ROOT / "experiments" / "protocol_v1")
    parser.add_argument("--keep-data-cache", action="store_true")
    parser.add_argument("--fail-fast", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    invalid_seeds = sorted(set(args.seeds) - set(SEEDS))
    if invalid_seeds:
        raise ValueError(f"Seeds outside preregistration: {invalid_seeds}")
    args.cache_dir.mkdir(parents=True, exist_ok=True)
    args.staging_root.mkdir(parents=True, exist_ok=True)
    args.status_dir.mkdir(parents=True, exist_ok=True)
    if args.dry_run:
        matrix = [
            {"dataset": dataset, "variant": variant, "seed": seed}
            for dataset in args.datasets for variant in args.variants for seed in args.seeds
        ]
        write_json(args.status_dir / "planned_matrix.json", matrix)
        print(json.dumps({"tasks": len(matrix), "matrix": matrix}, indent=2))
        return 0

    store = RemoteStore()
    store.ensure_layout()
    manifest = load_manifest(args.manifest)
    rows: list[dict] = []
    for dataset in args.datasets:
        require_disk_space(args.cache_dir, 5.0)
        if dataset not in manifest:
            raise KeyError(f"Dataset {dataset} missing from {args.manifest}")
        local_data = args.cache_dir / f"{dataset}.h5ad"
        expected = manifest[dataset]["sha256"]
        if not local_data.is_file() or sha256_file(local_data) != expected:
            local_data.unlink(missing_ok=True)
            store.download_file(f"{REMOTE_DATA_ROOT}/datasets/confirmatory_v1/{dataset}.h5ad", local_data)
        if sha256_file(local_data) != expected:
            raise ValueError(f"Downloaded dataset checksum failed: {dataset}")

        tasks = [(variant, seed) for variant in args.variants for seed in args.seeds]
        executors = {gpu: ThreadPoolExecutor(max_workers=1) for gpu in args.gpus}
        try:
            futures = {}
            for index, (variant, seed) in enumerate(tasks):
                gpu = args.gpus[index % len(args.gpus)]
                future = executors[gpu].submit(
                    run_task, dataset, local_data, variant, seed, gpu,
                    args.staging_root, args.epochs,
                )
                futures[future] = (variant, seed, gpu)
            for future in as_completed(futures):
                row = future.result()
                rows.append(row)
                write_status(rows, args.status_dir)
                if row["status"] == "failed" and args.fail_fast:
                    raise RuntimeError(row["message"])
        finally:
            for executor in executors.values():
                executor.shutdown(wait=True, cancel_futures=True)
        if any(row["status"] == "failed" for row in rows if row["dataset"] == dataset):
            # Preserve both the dataset cache and failed run staging, and do
            # not schedule another dataset after an execution/upload failure.
            return 1
        if not args.keep_data_cache:
            local_data.unlink(missing_ok=True)
    write_status(rows, args.status_dir)
    failures = [row for row in rows if row["status"] != "complete"]
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
