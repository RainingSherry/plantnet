#!/usr/bin/env python
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time
from pathlib import Path

import yaml

SCRIPT_DIR = Path(__file__).resolve().parent
PKG_DIR = SCRIPT_DIR.parent
ROOT = SCRIPT_DIR.parents[3]


def parse_args():
    parser = argparse.ArgumentParser(description="Run the 8-plant PlantSPADE-LGCL experiment suite.")
    parser.add_argument("--datasets_config", default=str(PKG_DIR / "configs" / "datasets_8plant.yaml"))
    parser.add_argument("--main_config", default=str(PKG_DIR / "configs" / "main_lgcl.yaml"))
    parser.add_argument("--baselines_config", default=str(PKG_DIR / "configs" / "baselines.yaml"))
    parser.add_argument("--datasets", default="all", help="Comma-separated dataset names or all.")
    parser.add_argument("--groups", default=None, help="Comma-separated group names.")
    parser.add_argument("--methods", default=None, help="Comma-separated methods; defaults to suite_methods in main config.")
    parser.add_argument("--seeds", default=None, help="Comma-separated seeds; defaults to main config seeds.")
    parser.add_argument("--output_dir", default=None)
    parser.add_argument("--gpu", type=int, default=None)
    parser.add_argument("--gpus", default=None, help="Comma-separated GPU ids for parallel scheduling, e.g. 1,2,3,4,5,6.")
    parser.add_argument("--jobs", type=int, default=None, help="Total concurrent runs. Defaults to main config jobs, normally 1.")
    parser.add_argument("--jobs_per_gpu", type=int, default=None, help="Concurrent runs per GPU. Default comes from config or 1.")
    parser.add_argument("--max_parallel_datasets", type=int, default=None, help="Maximum distinct datasets processed at once. Default 1.")
    parser.add_argument("--no_cuda", action="store_true")
    parser.add_argument("--save_h5ad", action="store_true")
    parser.add_argument("--dry_run", action="store_true")
    parser.add_argument("--keep_going", action="store_true")
    return parser.parse_args()


def load_yaml(path: str) -> dict:
    with open(path, "r", encoding="utf-8") as handle:
        return yaml.safe_load(handle) or {}


def split_csv(value):
    if value is None:
        return None
    return [item.strip() for item in str(value).split(",") if item.strip()]


def parse_int_csv(value):
    parsed = split_csv(value)
    if parsed is None:
        return None
    return [int(item) for item in parsed]


def select_datasets(cfg: dict, datasets_arg: str, groups_arg: str | None):
    datasets = cfg.get("datasets", [])
    wanted_names = None if datasets_arg == "all" else set(split_csv(datasets_arg))
    wanted_groups = set(split_csv(groups_arg) or [])
    out = []
    for entry in datasets:
        if wanted_names is not None and entry["dataset_name"] not in wanted_names:
            continue
        if wanted_groups and entry.get("group") not in wanted_groups:
            continue
        out.append(entry)
    return out


def resolve_gpus(args, main_cfg: dict):
    if args.no_cuda:
        return []
    if args.gpus:
        gpus = parse_int_csv(args.gpus)
    elif args.gpu is not None:
        gpus = [int(args.gpu)]
    elif "gpus" in main_cfg:
        gpus = [int(gpu) for gpu in main_cfg.get("gpus", [])]
    else:
        gpus = [int(main_cfg.get("gpu", 1))]
    if not gpus:
        raise ValueError("No GPU ids configured. Use --gpus, --gpu, or --no_cuda.")
    forbidden = sorted(set(gpus).intersection({0, 7}))
    if forbidden:
        print(f"Warning: GPUs {forbidden} are forbidden (system reserved); removing them from the pool.")
        gpus = [g for g in gpus if g not in {0, 7}]
        if not gpus:
            raise ValueError("All requested GPUs are forbidden. Use --gpus 1,2,3,4,5,6 or --no_cuda.")
    return gpus


def build_tasks(args, datasets, methods, seeds):
    tasks = []
    for entry in datasets:
        for method in methods:
            for seed in seeds:
                tasks.append(
                    {
                        "dataset": entry["dataset_name"],
                        "method": method,
                        "seed": int(seed),
                    }
                )
    return tasks


def output_root(args, main_cfg: dict) -> Path:
    return Path(args.output_dir or main_cfg.get("output_dir", ROOT / "results" / "PlantSPADE_LGCL_protocol"))


def task_output_dir(base_output: Path, task: dict) -> Path:
    return base_output / task["dataset"] / task["method"] / f"seed_{task['seed']}"


def build_command(args, task: dict, gpu: int | None):
    cmd = [
        sys.executable,
        str(SCRIPT_DIR / "run_single.py"),
        "--datasets_config",
        args.datasets_config,
        "--main_config",
        args.main_config,
        "--baselines_config",
        args.baselines_config,
        "--dataset",
        task["dataset"],
        "--method",
        task["method"],
        "--seed",
        str(task["seed"]),
    ]
    if args.output_dir:
        cmd.extend(["--output_dir", args.output_dir])
    if args.no_cuda:
        cmd.append("--no_cuda")
    elif gpu is not None:
        cmd.extend(["--gpu", "0"])
    if args.save_h5ad:
        cmd.append("--save_h5ad")
    if args.dry_run:
        cmd.append("--dry_run")
    return cmd


def command_for_log(cmd: list[str], env: dict | None = None) -> str:
    if not env:
        return " ".join(cmd)
    keys = [
        "CUDA_VISIBLE_DEVICES",
        "OMP_NUM_THREADS",
        "MKL_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "NUMBA_NUM_THREADS",
        "VECLIB_MAXIMUM_THREADS",
    ]
    prefix = [f"{key}={env[key]}" for key in keys if key in env]
    return " ".join(prefix + cmd)


def make_slots(args, main_cfg: dict, gpus: list[int]):
    if args.no_cuda:
        jobs = int(args.jobs or main_cfg.get("jobs", 1))
        return [None] * max(1, jobs)

    if args.jobs is not None:
        jobs = int(args.jobs)
    elif args.jobs_per_gpu is not None:
        jobs_per_gpu = int(args.jobs_per_gpu)
        if jobs_per_gpu < 1:
            raise ValueError("--jobs_per_gpu must be >= 1")
        jobs = len(gpus) * jobs_per_gpu
    else:
        jobs = int(main_cfg.get("jobs", 1))
    if jobs < 1:
        raise ValueError("--jobs must be >= 1")
    return [gpus[idx % len(gpus)] for idx in range(jobs)]


def thread_limited_env(slot_gpu: int | None) -> dict:
    env = os.environ.copy()
    env.update(
        {
            "OMP_NUM_THREADS": "1",
            "MKL_NUM_THREADS": "1",
            "OPENBLAS_NUM_THREADS": "1",
            "NUMBA_NUM_THREADS": "1",
            "VECLIB_MAXIMUM_THREADS": "1",
            "PYTHONUNBUFFERED": "1",
        }
    )
    if slot_gpu is not None:
        env["CUDA_VISIBLE_DEVICES"] = str(slot_gpu)
        env["PLANTSPADE_ASSIGNED_GPU"] = str(slot_gpu)
    return env


def acquire_gpu_lock(locks_dir: Path, gpu: int | None, task: dict) -> Path | None:
    if gpu is None:
        return None
    locks_dir.mkdir(parents=True, exist_ok=True)
    lock_path = locks_dir / f"gpu_{gpu}.lock"
    payload = {
        "pid": os.getpid(),
        "gpu": int(gpu),
        "task": task,
        "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
    }
    try:
        fd = os.open(str(lock_path), os.O_CREAT | os.O_EXCL | os.O_WRONLY)
    except FileExistsError:
        try:
            with open(lock_path, "r", encoding="utf-8") as handle:
                existing = json.load(handle)
            pid = int(existing.get("pid", -1))
            if pid > 0:
                os.kill(pid, 0)
        except ProcessLookupError:
            lock_path.unlink(missing_ok=True)
            fd = os.open(str(lock_path), os.O_CREAT | os.O_EXCL | os.O_WRONLY)
        except Exception:
            raise FileExistsError(str(lock_path))
        else:
            raise FileExistsError(str(lock_path))
    with os.fdopen(fd, "w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2)
    return lock_path


def release_gpu_lock(lock_path: Path | None) -> None:
    if lock_path is None:
        return
    try:
        lock_path.unlink()
    except FileNotFoundError:
        pass


def pop_launchable_task(pending: list[dict], running: list[dict], max_parallel_datasets: int) -> dict | None:
    if not pending:
        return None
    active = {item["task"]["dataset"] for item in running}
    if not active or len(active) < max_parallel_datasets:
        return pending.pop(0)
    for idx, task in enumerate(pending):
        if task["dataset"] in active:
            return pending.pop(idx)
    return None


def write_failure_report(base_output: Path, task: dict, gpu: int | None, returncode: int, log_path: Path) -> Path:
    run_dir = task_output_dir(base_output, task)
    run_dir.mkdir(parents=True, exist_ok=True)
    report = {
        "dataset": task["dataset"],
        "method": task["method"],
        "seed": task["seed"],
        "gpu": gpu,
        "returncode": int(returncode),
        "log": str(log_path),
        "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
        "recovery_attempt": bool(task.get("recovery_attempt", False)),
    }
    path = run_dir / "failure_report.json"
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(report, handle, indent=2)
    return path


def has_training_artifacts(base_output: Path, task: dict) -> bool:
    run_dir = task_output_dir(base_output, task)
    if task["method"].startswith("plantspade_lgcl"):
        required = [
            "embedding_baseline.npy",
            "training_history.json",
            "labels.npy",
            "support_matrix.npz",
            "amplitude_matrix.npz",
            "gene_embedding.npy",
            "global_embedding_svd_projected.npy",
        ]
        return all((run_dir / name).exists() for name in required)
    return (run_dir / "embedding_final.npy").exists()


def run_parallel(args, tasks: list[dict], slots: list[int | None], main_cfg: dict):
    base_output = output_root(args, main_cfg)
    base_output.mkdir(parents=True, exist_ok=True)
    locks_dir = base_output / ".locks"
    max_parallel_datasets = int(
        args.max_parallel_datasets
        if args.max_parallel_datasets is not None
        else main_cfg.get("max_parallel_datasets", 1)
    )
    if max_parallel_datasets < 1:
        raise ValueError("--max_parallel_datasets must be >= 1")
    failures = []
    running = []
    pending = list(tasks)
    free_slots = list(enumerate(slots))

    print(f"Parallel slots: {len(slots)} ({'CPU' if args.no_cuda else 'GPUs ' + ','.join(map(str, slots))})")
    print(f"Max parallel datasets: {max_parallel_datasets}")
    while pending or running:
        while pending and free_slots:
            task = pop_launchable_task(pending, running, max_parallel_datasets)
            if task is None:
                break
            slot_id, slot_gpu = free_slots.pop(0)
            lock_path = None
            try:
                lock_path = acquire_gpu_lock(locks_dir, slot_gpu, task)
            except FileExistsError:
                pending.insert(0, task)
                free_slots.append((slot_id, slot_gpu))
                break
            cmd = build_command(args, task, slot_gpu)
            env = thread_limited_env(slot_gpu)
            run_dir = task_output_dir(base_output, task)
            log_dir = run_dir / "logs"
            log_dir.mkdir(parents=True, exist_ok=True)
            log_path = log_dir / "run_suite.log"
            command_text = command_for_log(cmd, env)
            print(
                f"[launch] slot={slot_id} gpu={slot_gpu if slot_gpu is not None else 'cpu'} "
                f"{task['dataset']} {task['method']} seed={task['seed']}"
            )
            print(command_text)
            if args.dry_run:
                release_gpu_lock(lock_path)
                free_slots.append((slot_id, slot_gpu))
                continue

            log_handle = open(log_path, "a", encoding="utf-8")
            log_handle.write(f"\n\n===== run_suite launch {time.strftime('%Y-%m-%d %H:%M:%S')} =====\n")
            log_handle.write(command_text + "\n")
            log_handle.flush()
            proc = subprocess.Popen(
                cmd,
                cwd=str(ROOT),
                stdout=log_handle,
                stderr=subprocess.STDOUT,
                env=env,
            )
            running.append(
                {
                    "task": task,
                    "slot_id": slot_id,
                    "gpu": slot_gpu,
                    "proc": proc,
                    "log_handle": log_handle,
                    "log_path": log_path,
                    "lock_path": lock_path,
                }
            )

        if args.dry_run:
            break

        time.sleep(2.0)
        still_running = []
        for item in running:
            returncode = item["proc"].poll()
            if returncode is None:
                still_running.append(item)
                continue

            item["log_handle"].close()
            task = item["task"]
            release_gpu_lock(item.get("lock_path"))
            free_slots.append((item["slot_id"], item["gpu"]))
            status = "ok" if returncode == 0 else f"failed rc={returncode}"
            print(
                f"[done] {status} slot={item['slot_id']} gpu={item['gpu'] if item['gpu'] is not None else 'cpu'} "
                f"{task['dataset']} {task['method']} seed={task['seed']} log={item['log_path']}"
            )
            if returncode != 0:
                report_path = write_failure_report(base_output, task, item["gpu"], returncode, item["log_path"])
                failures.append(
                    {
                        "dataset": task["dataset"],
                        "method": task["method"],
                        "seed": task["seed"],
                        "gpu": item["gpu"],
                        "returncode": returncode,
                        "log": str(item["log_path"]),
                        "failure_report": str(report_path),
                    }
                )
                if returncode in {-11, 139}:
                    print("[SIGSEGV] failure report written; waiting 30 seconds before recovery scheduling")
                    time.sleep(30)
                if returncode in {-11, 139} and not task.get("recovery_attempt") and has_training_artifacts(base_output, task):
                    recovery_task = dict(task)
                    recovery_task["recovery_attempt"] = True
                    pending.insert(0, recovery_task)
                    failures.pop()
                    print(
                        f"[recovery queued] embeddings exist; run_single will skip training and run CPU eval "
                        f"for {task['dataset']} {task['method']} seed={task['seed']}"
                    )
                    continue
                if not args.keep_going:
                    for other in still_running:
                        other["proc"].terminate()
                        other["log_handle"].close()
                        release_gpu_lock(other.get("lock_path"))
                    raise SystemExit(f"Stopping after failed run: {failures[-1]}")
        running = still_running

    return failures


def main():
    args = parse_args()
    datasets_cfg = load_yaml(args.datasets_config)
    main_cfg = load_yaml(args.main_config)
    datasets = select_datasets(datasets_cfg, args.datasets, args.groups)
    if args.datasets == "all" and args.groups is None:
        excluded = set(main_cfg.get("exclude_from_full_suite", []))
        if excluded:
            before = len(datasets)
            datasets = [entry for entry in datasets if entry.get("dataset_name") not in excluded]
            removed = before - len(datasets)
            if removed:
                print(f"Excluded stress-test datasets from full suite: {sorted(excluded)}")
    methods = split_csv(args.methods)
    if not methods or (len(methods) == 1 and methods[0].lower() == "all"):
        methods = main_cfg.get("suite_methods", [])
    seeds = [int(seed) for seed in (split_csv(args.seeds) or main_cfg.get("seeds", [1, 2, 3, 4, 5]))]
    gpus = resolve_gpus(args, main_cfg)
    slots = make_slots(args, main_cfg, gpus)
    tasks = build_tasks(args, datasets, methods, seeds)

    total = len(tasks)
    print(f"Scheduled {total} runs: {len(datasets)} datasets x {len(methods)} methods x {len(seeds)} seeds")
    failures = run_parallel(args, tasks, slots, main_cfg)
    if failures:
        print("Failures:")
        for failure in failures:
            print(failure)
        raise SystemExit(1)


if __name__ == "__main__":
    main()
