#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]
BASE_DIR = Path(__file__).resolve().parent
RUN_ROOT = PROJECT_ROOT / "results" / "canonical" / "formal_varfloor_scmae_pca_20260704"
REFERENCE_CSV = PROJECT_ROOT / "results" / "260629全benchmark结果.csv"
SEEDS = [42, 2024, 3407]

DATA_ROOT = Path("/data/luolie/biopipeline/scCluBench/data")
DATASETS = [
    ("Bach", DATA_ROOT / "processed_scmae" / "Bach.h5ad", "resolved_label", 8),
    ("Baron", DATA_ROOT / "processed_scmae" / "Baron.h5ad", "resolved_label", 14),
    ("CRA002977_1", DATA_ROOT / "CRA002977_1.h5ad", "Celltype", 7),
    ("CRA007122", DATA_ROOT / "CRA007122.h5ad", "Celltype", 7),
    ("Guo", DATA_ROOT / "processed_scmae" / "Guo.h5ad", "resolved_label", 12),
    ("Limb_Muscle", DATA_ROOT / "processed_scmae" / "Limb_Muscle.h5ad", "resolved_label", 6),
    ("Macosko", DATA_ROOT / "processed_scmae" / "Macosko.h5ad", "resolved_label", 12),
    ("Melanoma_5K", DATA_ROOT / "processed_scmae" / "Melanoma_5K.h5ad", "resolved_label", 9),
    ("Pollen", DATA_ROOT / "processed" / "Pollen.h5ad", "resolved_label", 11),
    ("Quake_10x_Spleen", DATA_ROOT / "processed_scmae" / "Quake_10x_Spleen.h5ad", "resolved_label", 5),
    ("Quake_Smart-seq2_Lung", DATA_ROOT / "processed" / "Quake_Smart-seq2_Lung.h5ad", "resolved_label", 11),
    ("SRP145013", DATA_ROOT / "SRP145013.h5ad", "Celltype", 9),
    ("SRP171040", DATA_ROOT / "SRP171040.h5ad", "Celltype", 12),
    ("SRP182008", DATA_ROOT / "SRP182008.h5ad", "Celltype", 15),
    ("SRP224648", DATA_ROOT / "SRP224648.h5ad", "Celltype", 4),
    ("SRP235541", DATA_ROOT / "SRP235541.h5ad", "Celltype", 18),
    ("SRP309176", DATA_ROOT / "SRP309176.h5ad", "Celltype", 13),
    ("Shekhar", DATA_ROOT / "processed_scmae" / "Shekhar.h5ad", "resolved_label", 18),
    ("Tosches", DATA_ROOT / "processed_scmae" / "Tosches.h5ad", "resolved_label", 15),
    ("Wang", DATA_ROOT / "processed" / "Wang.h5ad", "resolved_label", 2),
    ("Young", DATA_ROOT / "processed_scmae" / "Young.h5ad", "resolved_label", 11),
    ("hrvatin_geo_maintype_counts", DATA_ROOT / "hrvatin_geo_maintype_counts.h5ad", "maintype", 8),
    ("worm_neuron_cell", DATA_ROOT / "processed_scmae" / "worm_neuron_cell.h5ad", "resolved_label", 10),
]

METHODS = {
    "varfloor_scmae": {
        "display": "VarFloor-scMAE",
        "method_raw": "varfloor_scmae",
        "runner": PROJECT_ROOT / "methods" / "DeepLearning" / "scMAE_DEC_StdFloor" / "run.py",
        "uses_gpu": True,
    },
    "pca_kmeans_known_k": {
        "display": "PCA+KMeans known-K",
        "method_raw": "pca_kmeans_known_k",
        "runner": PROJECT_ROOT / "methods" / "Traditional" / "PCA_KMeans" / "run.py",
        "uses_gpu": False,
    },
}


def parse_csv_list(value: str) -> list[str]:
    return [x.strip() for x in str(value).split(",") if x.strip()]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--methods", default="varfloor_scmae,pca_kmeans_known_k")
    parser.add_argument("--datasets", default="all")
    parser.add_argument("--seeds", default="42,2024,3407")
    parser.add_argument("--gpus", default="1,2,3,4,5,6")
    parser.add_argument("--max_parallel", type=int, default=3)
    parser.add_argument("--dry_run", action="store_true")
    parser.add_argument("--merge_only", action="store_true")
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--allow_partial", action="store_true")
    parser.add_argument("--update_reference", action="store_true")
    parser.add_argument("--run_root", default=str(RUN_ROOT))
    parser.add_argument("--reference_csv", default=str(REFERENCE_CSV))
    return parser.parse_args()


def env_for_task() -> dict[str, str]:
    env = os.environ.copy()
    env.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")
    env.setdefault("NUMBA_CACHE_DIR", "/tmp/numba-cache")
    env.setdefault("OPENBLAS_NUM_THREADS", "1")
    env.setdefault("OMP_NUM_THREADS", "1")
    env.setdefault("MKL_NUM_THREADS", "1")
    env.setdefault("NUMEXPR_NUM_THREADS", "1")
    return env


def dataset_rows(selected: list[str]) -> list[dict]:
    rows = []
    wanted = None if selected == ["all"] else set(selected)
    for dataset, path, label_key, n_clusters in DATASETS:
        if wanted is not None and dataset not in wanted:
            continue
        if not path.exists():
            raise FileNotFoundError(f"{dataset}: missing data file {path}")
        rows.append(
            {
                "dataset": dataset,
                "data_path": str(path),
                "label_key": label_key,
                "n_clusters": int(n_clusters),
            }
        )
    if not rows:
        raise SystemExit(f"No datasets selected: {selected}")
    return rows


def build_command(method_key: str, ds: dict, seed: int, gpu: int | None, run_dir: Path) -> list[str]:
    method = METHODS[method_key]
    cmd = [
        sys.executable,
        str(method["runner"]),
        "--data_path",
        ds["data_path"],
        "--save_dir",
        str(run_dir),
        "--dataset_name",
        ds["dataset"],
        "--label_key",
        ds["label_key"],
        "--input_mode",
        "auto",
        "--n_top_genes",
        "1000",
        "--target_sum",
        "10000",
        "--scale_input",
        "true",
        "--n_clusters",
        str(ds["n_clusters"]),
        "--seed",
        str(seed),
        "--method_name",
        method["display"],
        "--variant_name",
        method["method_raw"],
    ]
    if method_key == "varfloor_scmae":
        cmd.extend(["--epochs", "80", "--batch_size", "256"])
        if gpu is not None:
            cmd.extend(["--gpu", str(gpu)])
    elif method_key == "pca_kmeans_known_k":
        cmd.extend(["--pca_dim", "128"])
    return cmd


def task_status_path(run_dir: Path) -> Path:
    return run_dir / "status.json"


def is_success(run_dir: Path) -> bool:
    status_path = task_status_path(run_dir)
    eval_path = run_dir / "eval_fixed.csv"
    if not status_path.exists() or not eval_path.exists():
        return False
    try:
        return json.loads(status_path.read_text()).get("status") == "success"
    except Exception:
        return False


def write_json(path: Path, payload: dict | list) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def build_tasks(args: argparse.Namespace) -> list[dict]:
    method_keys = parse_csv_list(args.methods)
    unknown = sorted(set(method_keys) - set(METHODS))
    if unknown:
        raise SystemExit(f"Unknown methods: {unknown}")
    datasets = dataset_rows(parse_csv_list(args.datasets))
    seeds = [int(x) for x in parse_csv_list(args.seeds)]
    gpus = [int(x) for x in parse_csv_list(args.gpus)]
    if any(g in {0, 7} for g in gpus):
        raise SystemExit("GPU 0 and GPU 7 are forbidden.")

    run_root = Path(args.run_root)
    tasks = []
    gpu_idx = 0
    for ds in datasets:
        for method_key in method_keys:
            for seed in seeds:
                method = METHODS[method_key]
                gpu = None
                if method["uses_gpu"]:
                    gpu = gpus[gpu_idx % len(gpus)]
                    gpu_idx += 1
                run_dir = run_root / ds["dataset"] / f"{method_key}__seed{seed}"
                cmd = build_command(method_key, ds, seed, gpu, run_dir)
                tasks.append(
                    {
                        "dataset": ds["dataset"],
                        "data_path": ds["data_path"],
                        "label_key": ds["label_key"],
                        "n_clusters": ds["n_clusters"],
                        "method_key": method_key,
                        "method": method["display"],
                        "method_raw": method["method_raw"],
                        "seed": int(seed),
                        "gpu": gpu,
                        "run_dir": str(run_dir),
                        "command": cmd,
                    }
                )
    return tasks


def run_tasks(tasks: list[dict], args: argparse.Namespace) -> None:
    run_root = Path(args.run_root)
    write_json(BASE_DIR / "task_manifest.json", tasks)
    if args.dry_run:
        print(f"Dry run: {len(tasks)} tasks")
        for task in tasks[:20]:
            print(" ".join(task["command"]))
        if len(tasks) > 20:
            print(f"... {len(tasks) - 20} more")
        return

    pending = []
    for task in tasks:
        run_dir = Path(task["run_dir"])
        if not args.force and is_success(run_dir):
            continue
        pending.append(task)

    print(f"Tasks total={len(tasks)} pending={len(pending)} run_root={run_root}")
    active: list[dict] = []
    completed = 0
    failures = 0
    env = env_for_task()

    def launch(task: dict) -> dict:
        run_dir = Path(task["run_dir"])
        run_dir.mkdir(parents=True, exist_ok=True)
        log_handle = (run_dir / "run.log").open("w", encoding="utf-8")
        (run_dir / "command.txt").write_text(" ".join(task["command"]) + "\n", encoding="utf-8")
        write_json(run_dir / "run_config.json", {k: v for k, v in task.items() if k != "command"})
        start = time.time()
        status = {**task, "status": "running", "start_time": datetime.now().isoformat()}
        write_json(task_status_path(run_dir), {k: v for k, v in status.items() if k != "command"})
        proc = subprocess.Popen(task["command"], cwd=PROJECT_ROOT, env=env, stdout=log_handle, stderr=subprocess.STDOUT)
        return {"task": task, "proc": proc, "log_handle": log_handle, "start": start}

    i = 0
    max_parallel = max(1, int(args.max_parallel))
    while i < len(pending) or active:
        while i < len(pending) and len(active) < max_parallel:
            task = pending[i]
            i += 1
            print(f"[launch] {task['dataset']} {task['method_raw']} seed={task['seed']} gpu={task['gpu']}")
            active.append(launch(task))

        time.sleep(5)
        still_active = []
        for item in active:
            proc = item["proc"]
            rc = proc.poll()
            if rc is None:
                still_active.append(item)
                continue
            item["log_handle"].close()
            task = item["task"]
            elapsed = time.time() - item["start"]
            success = rc == 0 and (Path(task["run_dir"]) / "eval_fixed.csv").exists()
            status = {
                **task,
                "status": "success" if success else "failed",
                "return_code": int(rc),
                "elapsed_seconds": float(elapsed),
                "end_time": datetime.now().isoformat(),
            }
            write_json(task_status_path(Path(task["run_dir"])), {k: v for k, v in status.items() if k != "command"})
            completed += 1
            failures += 0 if success else 1
            print(f"[done] {task['dataset']} {task['method_raw']} seed={task['seed']} status={status['status']} elapsed={elapsed:.1f}s")
        active = still_active

    print(f"Completed pending tasks={completed}; failures={failures}")
    if failures and not args.allow_partial:
        raise SystemExit("At least one run failed; not aggregating as complete.")


def read_run_row(task: dict) -> dict | None:
    run_dir = Path(task["run_dir"])
    eval_path = run_dir / "eval_fixed.csv"
    if not eval_path.exists():
        return None
    row = pd.read_csv(eval_path).iloc[0].to_dict()
    summary_path = run_dir / "summary.json"
    runtime = None
    if summary_path.exists():
        try:
            runtime = json.loads(summary_path.read_text()).get("runtime_seconds")
        except Exception:
            runtime = None
    return {
        "dataset": task["dataset"],
        "method": task["method"],
        "method_raw": task["method_raw"],
        "seed": int(task["seed"]),
        "run_dir": str(run_dir),
        "acc": row.get("acc"),
        "nmi": row.get("nmi"),
        "ari": row.get("ari"),
        "f1_macro": row.get("f1_macro"),
        "runtime_seconds": runtime,
    }


def aggregate(tasks: list[dict]) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows = [row for task in tasks if (row := read_run_row(task)) is not None]
    runs = pd.DataFrame(rows)
    if runs.empty:
        raise SystemExit("No eval_fixed.csv files found to aggregate.")
    runs.to_csv(BASE_DIR / "formal_completion_runs.csv", index=False)

    out_rows = []
    for (dataset, method, method_raw), sub in runs.groupby(["dataset", "method", "method_raw"], sort=True):
        out = {"dataset": dataset, "method": method}
        for metric in ["acc", "nmi", "ari", "f1_macro"]:
            vals = pd.to_numeric(sub[metric], errors="coerce").dropna()
            out[f"{metric}_mean"] = vals.mean() if len(vals) else pd.NA
            out[f"{metric}_std"] = vals.std(ddof=1) if len(vals) > 1 else pd.NA
            out[f"{metric}_count"] = int(len(vals))
        vals = pd.to_numeric(sub["runtime_seconds"], errors="coerce").dropna()
        out["runtime_seconds_mean"] = vals.mean() if len(vals) else pd.NA
        out["runtime_seconds_std"] = vals.std(ddof=1) if len(vals) > 1 else pd.NA
        seeds = sorted({int(x) for x in sub["seed"].dropna().tolist()})
        out["seed_<lambda>"] = ",".join(str(x) for x in seeds)
        out["n_seeds"] = int(len(seeds))
        out["method_raw"] = method_raw
        out_rows.append(out)
    summary = pd.DataFrame(out_rows)
    summary.to_csv(BASE_DIR / "formal_completion_summary.csv", index=False)
    return runs, summary


def merge_reference(summary: pd.DataFrame, reference_csv: Path, update_reference: bool) -> Path:
    ref = pd.read_csv(reference_csv)
    method_raws = set(summary["method_raw"].astype(str))
    displays = set(summary["method"].astype(str))
    kept = ref[~(ref["method_raw"].astype(str).isin(method_raws) | ref["method"].astype(str).isin(displays))]
    merged = pd.concat([kept, summary[ref.columns]], ignore_index=True)
    out_path = BASE_DIR / "260629全benchmark结果.with_varfloor_pca.csv"
    merged.to_csv(out_path, index=False)
    if update_reference:
        stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        backup = reference_csv.with_name(reference_csv.name + f".bak_{stamp}")
        shutil.copy2(reference_csv, backup)
        merged.to_csv(reference_csv, index=False)
        print(f"Updated reference CSV: {reference_csv}")
        print(f"Backup: {backup}")
    return out_path


def validate_completeness(tasks: list[dict], summary: pd.DataFrame, allow_partial: bool) -> None:
    expected = {(t["dataset"], t["method_raw"], int(t["seed"])) for t in tasks}
    actual = set()
    for task in tasks:
        if is_success(Path(task["run_dir"])):
            actual.add((task["dataset"], task["method_raw"], int(task["seed"])))
    missing = sorted(expected - actual)
    write_json(BASE_DIR / "formal_completion_missing.json", missing)
    if missing and not allow_partial:
        raise SystemExit(f"Missing successful runs: {len(missing)}. See formal_completion_missing.json")
    print("Summary rows:")
    print(summary[["dataset", "method_raw", "ari_mean", "ari_std", "nmi_mean", "nmi_std", "n_seeds"]].to_string(index=False))


def main() -> int:
    args = parse_args()
    tasks = build_tasks(args)
    if not args.merge_only:
        run_tasks(tasks, args)
        if args.dry_run:
            return 0
    runs, summary = aggregate(tasks)
    validate_completeness(tasks, summary, args.allow_partial)
    out_path = merge_reference(summary, Path(args.reference_csv), args.update_reference)
    print(f"Run rows: {BASE_DIR / 'formal_completion_runs.csv'} ({len(runs)})")
    print(f"Summary: {BASE_DIR / 'formal_completion_summary.csv'} ({len(summary)})")
    print(f"Merged CSV: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

