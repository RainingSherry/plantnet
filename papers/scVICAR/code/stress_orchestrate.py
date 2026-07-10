from __future__ import annotations

import argparse
import csv
import json
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

from .config import PAPER_ROOT, REMOTE_DATA_ROOT, SEEDS
from .io_utils import require_disk_space, sha256_file, utc_now, write_json
from .remote_store import RemoteStore
from .run_stress import CONTAMINATION_RATIOS, ESTIMATORS, STRESS_DATASETS, STRESS_VARIANTS


def tasks_for_dataset(dataset: str) -> list[tuple[str, int, float, str]]:
    tasks = {
        (variant, seed, ratio, "current")
        for variant in STRESS_VARIANTS for seed in SEEDS for ratio in CONTAMINATION_RATIOS
    }
    tasks.update(
        (variant, seed, 0.0, estimator)
        for variant in STRESS_VARIANTS for seed in SEEDS for estimator in ESTIMATORS
    )
    return sorted(tasks, key=lambda item: (item[2], item[3], item[0], item[1]))


def run_task(dataset: str, data_path: Path, task: tuple[str, int, float, str], gpu: int, epochs: int | None) -> dict:
    variant, seed, ratio, estimator = task
    command = [
        sys.executable, "-m", "papers.scVICAR.code.run_stress",
        "--dataset", dataset, "--data-path", str(data_path), "--variant", variant,
        "--seed", str(seed), "--contamination", str(ratio), "--estimator", estimator,
        "--gpu", str(gpu), "--upload", "--cleanup-after-upload",
    ]
    if epochs is not None:
        command.extend(["--epochs", str(epochs)])
    started = utc_now()
    result = subprocess.run(command, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    return {
        "dataset": dataset, "variant": variant, "seed": seed, "contamination": ratio,
        "estimator": estimator, "gpu": gpu, "started_utc": started,
        "finished_utc": utc_now(), "returncode": result.returncode,
        "status": "complete" if result.returncode == 0 else "failed",
        "message": result.stdout[-4000:],
    }


def write_status(rows: list[dict], target: Path) -> None:
    target.mkdir(parents=True, exist_ok=True)
    write_json(target / "stress_master.json", rows)
    if rows:
        with (target / "stress_master.csv").open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
            writer.writeheader(); writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description="Run the frozen scVICAR contamination and estimator matrix")
    parser.add_argument("--datasets", nargs="+", default=list(STRESS_DATASETS), choices=STRESS_DATASETS)
    parser.add_argument("--gpus", nargs="+", type=int, default=[2, 3, 4, 5, 6], choices=range(1, 7))
    parser.add_argument("--epochs", type=int)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--fail-fast", action="store_true")
    parser.add_argument("--cache-dir", type=Path, default=PAPER_ROOT / ".staging/data")
    parser.add_argument("--status-dir", type=Path, default=PAPER_ROOT / "experiments/stress_v1")
    args = parser.parse_args()
    matrix = [
        {"dataset": dataset, "variant": v, "seed": s, "contamination": r, "estimator": e}
        for dataset in args.datasets for v, s, r, e in tasks_for_dataset(dataset)
    ]
    if args.dry_run:
        write_json(args.status_dir / "planned_matrix.json", matrix)
        print(json.dumps({"tasks": len(matrix)}, indent=2)); return 0
    manifest = {row["name"]: row for row in json.loads((PAPER_ROOT / "manifests/dataset_upload/dataset_manifest.json").read_text())}
    store = RemoteStore(); store.ensure_layout()
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
        tasks = tasks_for_dataset(dataset)
        executors = {gpu: ThreadPoolExecutor(max_workers=1) for gpu in args.gpus}
        try:
            futures = {
                executors[args.gpus[i % len(args.gpus)]].submit(
                    run_task, dataset, local, task, args.gpus[i % len(args.gpus)], args.epochs
                ): task for i, task in enumerate(tasks)
            }
            for future in as_completed(futures):
                row = future.result(); rows.append(row); write_status(rows, args.status_dir)
                if row["status"] == "failed" and args.fail_fast:
                    raise RuntimeError(row["message"])
        finally:
            for executor in executors.values():
                executor.shutdown(wait=True, cancel_futures=True)
        if any(row["status"] == "failed" for row in rows if row["dataset"] == dataset):
            return 1
        local.unlink(missing_ok=True)
    return 1 if any(row["status"] != "complete" for row in rows) else 0


if __name__ == "__main__":
    raise SystemExit(main())
