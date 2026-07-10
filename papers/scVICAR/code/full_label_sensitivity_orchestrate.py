from __future__ import annotations

import argparse
import csv
import json
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

from .config import PAPER_ROOT, SEEDS
from .io_utils import require_disk_space, sha256_file, utc_now, write_json
from .prepare_full_label_sensitivity import VERSION
from .remote_store import RemoteStore
from .run_full_label_sensitivity import SENSITIVITY_VARIANTS, manifest_lookup


def run_task(dataset: str, data_path: Path, variant: str, seed: int, gpu: int) -> dict:
    command = [
        sys.executable, "-m", "papers.scVICAR.code.run_full_label_sensitivity",
        "--dataset", dataset, "--data-path", str(data_path), "--variant", variant,
        "--seed", str(seed), "--gpu", str(gpu), "--upload", "--cleanup-after-upload",
    ]
    started = utc_now()
    result = subprocess.run(command, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    return {
        "dataset": dataset, "variant": variant, "seed": seed, "gpu": gpu,
        "started_utc": started, "finished_utc": utc_now(), "returncode": result.returncode,
        "status": "complete" if result.returncode == 0 else "failed",
        "message": result.stdout[-4000:],
    }


def write_status(rows: list[dict], output: Path) -> None:
    output.mkdir(parents=True, exist_ok=True)
    write_json(output / "scheduler_status.json", {"updated_utc": utc_now(), "results": rows})
    if rows:
        with (output / "scheduler_status.csv").open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
            writer.writeheader(); writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description="Run the frozen 54-task full-label sensitivity matrix")
    parser.add_argument("--datasets", nargs="+")
    parser.add_argument("--gpus", nargs="+", type=int, default=[2, 3, 4, 5, 6], choices=range(1, 7))
    parser.add_argument(
        "--manifest", type=Path,
        default=PAPER_ROOT / f"manifests/{VERSION}/dataset_manifest.json",
    )
    parser.add_argument(
        "--cache-dir", type=Path,
        default=PAPER_ROOT / f".staging/data_{VERSION}",
    )
    parser.add_argument(
        "--status-dir", type=Path,
        default=PAPER_ROOT / f"experiments/{VERSION}",
    )
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--fail-fast", action="store_true")
    args = parser.parse_args()
    manifest = manifest_lookup(args.manifest)
    datasets = args.datasets or list(manifest)
    if not set(datasets).issubset(manifest):
        raise ValueError("Requested dataset is absent from sensitivity manifest")
    matrix = [
        {"dataset": dataset, "variant": variant, "seed": seed}
        for dataset in datasets for variant in SENSITIVITY_VARIANTS for seed in SEEDS
    ]
    if args.dry_run:
        write_json(args.status_dir / "planned_matrix.json", matrix)
        print(json.dumps({"tasks": len(matrix)}, indent=2))
        return 0
    args.cache_dir.mkdir(parents=True, exist_ok=True)
    store = RemoteStore()
    rows = []
    for dataset in datasets:
        require_disk_space(args.cache_dir, 5.0)
        local = args.cache_dir / f"{dataset}.h5ad"
        expected = manifest[dataset]["sha256"]
        if not local.is_file() or sha256_file(local) != expected:
            local.unlink(missing_ok=True)
            store.download_file(manifest[dataset]["remote_path"], local)
        if sha256_file(local) != expected:
            raise ValueError(f"Sensitivity dataset checksum failed: {dataset}")
        tasks = [(variant, seed) for variant in SENSITIVITY_VARIANTS for seed in SEEDS]
        executors = {gpu: ThreadPoolExecutor(max_workers=1) for gpu in args.gpus}
        try:
            futures = {
                executors[args.gpus[index % len(args.gpus)]].submit(
                    run_task, dataset, local, variant, seed, args.gpus[index % len(args.gpus)]
                ): (variant, seed)
                for index, (variant, seed) in enumerate(tasks)
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
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
