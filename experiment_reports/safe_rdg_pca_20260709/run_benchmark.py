#!/usr/bin/env python3
from __future__ import annotations

import argparse
import concurrent.futures as futures
import json
import os
import subprocess
import sys
import time
from pathlib import Path
from typing import Any

import h5py
import numpy as np


PROJECT_ROOT = Path(__file__).resolve().parents[2]
BASE_DIR = Path(__file__).resolve().parent
DEFAULT_SOURCE = PROJECT_ROOT / "result" / "scmae_all_methods_20260705_full"
DEFAULT_DATA_ROOT = PROJECT_ROOT / "data" / "scMAE"
DEFAULT_PROCESSED_ROOT = DEFAULT_SOURCE / "converted_data"
RUNNER = PROJECT_ROOT / "experimental_retired_models" / "Safe_RDG_PCA" / "run.py"
DEFAULT_SEEDS = [42, 3047, 3407]
SMOKE_DATASETS = ["Pollen", "worm_neuron_cell", "Bach", "Macosko", "Wang", "Limb_Muscle"]


def read_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def parse_csv(value: str) -> list[str]:
    return [x.strip() for x in str(value).split(",") if x.strip()]


def n_clusters_from_scmae_h5(path: Path) -> int:
    with h5py.File(path, "r") as handle:
        labels = np.asarray(handle["Y"])
    if labels.dtype.kind in {"S", "O", "U"}:
        labels = np.asarray([v.decode("utf-8") if isinstance(v, bytes) else str(v) for v in labels], dtype=object)
    return int(len(np.unique(labels.astype(str))))


def catalog_by_name(source_root: Path) -> dict[str, dict]:
    path = source_root / "dataset_catalog.json"
    if not path.exists():
        return {}
    return {item["name"]: item for item in read_json(path)}


def discover_raw_scmae(data_root: Path) -> dict[str, Path]:
    rows = {}
    if not data_root.exists():
        return rows
    for path in sorted(data_root.glob("*.h5")):
        rows[path.stem] = path
    return rows


def load_datasets(source_root: Path, data_root: Path, processed_root: Path, names: list[str] | None, data_preference: str) -> list[dict]:
    catalog = catalog_by_name(source_root)
    raw = discover_raw_scmae(data_root)
    wanted = set(names) if names else None
    all_names = sorted(set(catalog) | set(raw))
    rows = []
    for name in all_names:
        if wanted and name not in wanted:
            continue
        item = catalog.get(name, {})
        raw_path = raw.get(name)
        processed_path = processed_root / f"{name}.h5ad"
        catalog_path = Path(item["prepared_path"]) if item.get("prepared_path") else None
        candidates = {
            "processed": processed_path if processed_path.exists() else None,
            "raw": raw_path if raw_path and raw_path.exists() else None,
            "catalog": catalog_path if catalog_path and catalog_path.exists() else None,
        }
        if data_preference == "processed":
            order = ["processed", "catalog", "raw"]
        elif data_preference == "raw":
            order = ["raw", "processed", "catalog"]
        else:
            order = ["catalog", "processed", "raw"]
        selected_kind = next((kind for kind in order if candidates.get(kind) is not None), None)
        if selected_kind is None:
            raise FileNotFoundError(f"No usable data file for {name} under {data_root}, {processed_root}, or {source_root}")
        selected_path = candidates[selected_kind]
        if item.get("n_clusters") is not None:
            n_clusters = int(item["n_clusters"])
        elif raw_path is not None:
            n_clusters = n_clusters_from_scmae_h5(raw_path)
        else:
            raise ValueError(f"n_clusters is absent for {name} and no raw scMAE h5 is available for inference.")
        rows.append(
            {
                "name": name,
                "path": selected_path,
                "n_clusters": n_clusters,
                "source_kind": selected_kind,
                "raw_path": str(raw_path) if raw_path else "",
                "processed_path": str(processed_path) if processed_path.exists() else "",
            }
        )
    if not rows:
        raise SystemExit("No datasets selected.")
    return rows


def build_cmd(args: argparse.Namespace, ds: dict, seed: int, run_dir: Path) -> list[str]:
    return [
        sys.executable,
        str(RUNNER),
        "--data_path",
        str(ds["path"]),
        "--save_dir",
        str(run_dir),
        "--dataset_name",
        ds["name"],
        "--n_clusters",
        str(ds["n_clusters"]),
        "--seed",
        str(seed),
        "--variant_name",
        "stage_a_all",
        "--method_name",
        "Safe-RDG-PCA",
        "--label_key",
        "auto",
        "--input_mode",
        "auto",
        "--n_top_genes",
        str(args.n_top_genes),
        "--target_sum",
        "10000",
        "--scale_input",
        "true",
        "--raw_pca_dim",
        str(args.raw_pca_dim),
        "--gene_bootstrap_B",
        str(args.gene_bootstrap_B),
        "--gene_edge_stability_threshold",
        str(args.gene_edge_stability_threshold),
        "--heuristic_threshold",
        str(args.heuristic_threshold),
        "--include_negative_controls",
        str(args.include_negative_controls).lower(),
    ]


def run_one(task: dict, args: argparse.Namespace) -> dict:
    ds = task["dataset"]
    seed = task["seed"]
    run_dir = Path(task["run_dir"])
    status_path = run_dir / "status.json"
    if args.resume and status_path.exists() and read_json(status_path).get("status") == "success":
        return read_json(status_path)
    run_dir.mkdir(parents=True, exist_ok=True)
    cmd = task["command"]
    (run_dir / "command.txt").write_text(" ".join(cmd) + "\n", encoding="utf-8")
    env = os.environ.copy()
    env.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")
    env.setdefault("NUMBA_CACHE_DIR", "/tmp/numba-cache")
    env.setdefault("OPENBLAS_NUM_THREADS", "1")
    env.setdefault("OMP_NUM_THREADS", "1")
    env.setdefault("MKL_NUM_THREADS", "1")
    env.setdefault("NUMEXPR_NUM_THREADS", "1")
    t0 = time.time()
    proc = subprocess.run(cmd, cwd=PROJECT_ROOT, env=env, text=True, capture_output=True, timeout=args.timeout)
    runtime = time.time() - t0
    (run_dir / "stdout.log").write_text(proc.stdout, encoding="utf-8")
    (run_dir / "stderr.log").write_text(proc.stderr, encoding="utf-8")
    status = {
        "dataset": ds,
        "seed": seed,
        "status": "success" if proc.returncode == 0 else "failed",
        "returncode": proc.returncode,
        "runtime_seconds": runtime,
        "command": " ".join(cmd),
        "error": proc.stderr[-4000:] if proc.returncode else "",
    }
    write_json(status_path, status)
    return status


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--source_root", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--data_root", type=Path, default=DEFAULT_DATA_ROOT)
    parser.add_argument("--processed_root", type=Path, default=DEFAULT_PROCESSED_ROOT)
    parser.add_argument("--data_preference", default="processed", choices=["processed", "raw", "catalog"])
    parser.add_argument("--out_root", type=Path, default=BASE_DIR / "runs")
    parser.add_argument("--datasets", default="all")
    parser.add_argument("--seeds", default="42,3047,3407")
    parser.add_argument("--max_workers", type=int, default=2)
    parser.add_argument("--timeout", type=int, default=7200)
    parser.add_argument("--resume", action="store_true", default=True)
    parser.add_argument("--no_resume", dest="resume", action="store_false")
    parser.add_argument("--smoke", action="store_true")
    parser.add_argument("--dry_run", action="store_true")
    parser.add_argument("--batch_count", type=int, default=1, help="Split selected tasks into this many deterministic modulo batches.")
    parser.add_argument("--batch_index", type=int, default=0, help="Run only tasks whose index modulo batch_count equals this value.")
    parser.add_argument("--max_tasks", type=int, default=0, help="Optionally cap the number of tasks after batching; 0 means no cap.")
    parser.add_argument("--n_top_genes", type=int, default=1000)
    parser.add_argument("--raw_pca_dim", type=int, default=128)
    parser.add_argument("--gene_bootstrap_B", type=int, default=20)
    parser.add_argument("--gene_edge_stability_threshold", type=float, default=0.5)
    parser.add_argument("--heuristic_threshold", type=float, default=0.45)
    parser.add_argument("--include_negative_controls", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    names = None if args.datasets == "all" else parse_csv(args.datasets)
    if args.smoke and names is None:
        names = SMOKE_DATASETS
    datasets = load_datasets(args.source_root, args.data_root, args.processed_root, names, args.data_preference)
    seeds = [int(x) for x in parse_csv(args.seeds)]
    tasks = []
    for ds in datasets:
        for seed in seeds:
            run_dir = args.out_root / ds["name"] / f"seed_{seed}"
            tasks.append({"dataset": ds["name"], "seed": seed, "run_dir": str(run_dir), "command": build_cmd(args, ds, seed, run_dir)})
    total_tasks = len(tasks)
    if args.batch_count < 1:
        raise SystemExit("--batch_count must be >= 1")
    if not 0 <= args.batch_index < args.batch_count:
        raise SystemExit("--batch_index must satisfy 0 <= batch_index < batch_count")
    if args.batch_count > 1:
        tasks = [task for idx, task in enumerate(tasks) if idx % args.batch_count == args.batch_index]
    if args.max_tasks and args.max_tasks > 0:
        tasks = tasks[: args.max_tasks]
    write_json(BASE_DIR / "task_manifest.json", tasks)
    write_json(args.out_root / "task_manifest.json", tasks)
    config = {key: str(value) if isinstance(value, Path) else value for key, value in vars(args).items()}
    config.update({"seeds": seeds, "datasets": [d["name"] for d in datasets], "total_tasks_before_batch": total_tasks, "tasks_after_batch": len(tasks)})
    write_json(BASE_DIR / "suite_config.json", config)
    write_json(args.out_root / "suite_config.json", config)
    print(f"Tasks: {len(tasks)}")
    if args.dry_run:
        for task in tasks:
            print(task["command"])
        return 0
    failures = 0
    with futures.ThreadPoolExecutor(max_workers=args.max_workers) as pool:
        futs = [pool.submit(run_one, task, args) for task in tasks]
        for fut in futures.as_completed(futs):
            status = fut.result()
            print(f"{status['dataset']} seed={status['seed']} {status['status']} {status['runtime_seconds']:.1f}s")
            failures += int(status["status"] != "success")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
