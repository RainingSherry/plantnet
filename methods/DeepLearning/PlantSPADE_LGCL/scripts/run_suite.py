#!/usr/bin/env python
from __future__ import annotations

import argparse
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
    parser.add_argument("--jobs", type=int, default=None, help="Total concurrent runs. Defaults to len(gpus) * jobs_per_gpu.")
    parser.add_argument("--jobs_per_gpu", type=int, default=None, help="Concurrent runs per GPU. Default comes from config or 1.")
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
        cmd.extend(["--gpu", str(gpu)])
    if args.save_h5ad:
        cmd.append("--save_h5ad")
    if args.dry_run:
        cmd.append("--dry_run")
    return cmd


def make_slots(args, main_cfg: dict, gpus: list[int]):
    if args.no_cuda:
        jobs = int(args.jobs or main_cfg.get("jobs", 1))
        return [None] * max(1, jobs)

    jobs_per_gpu = int(args.jobs_per_gpu or main_cfg.get("jobs_per_gpu", 1))
    if jobs_per_gpu < 1:
        raise ValueError("--jobs_per_gpu must be >= 1")
    if args.jobs is None:
        slots = []
        for _ in range(jobs_per_gpu):
            slots.extend(gpus)
        return slots

    jobs = int(args.jobs)
    if jobs < 1:
        raise ValueError("--jobs must be >= 1")
    return [gpus[idx % len(gpus)] for idx in range(jobs)]


def run_parallel(args, tasks: list[dict], slots: list[int | None], main_cfg: dict):
    base_output = output_root(args, main_cfg)
    base_output.mkdir(parents=True, exist_ok=True)
    failures = []
    running = []
    pending = list(tasks)
    free_slots = list(enumerate(slots))

    print(f"Parallel slots: {len(slots)} ({'CPU' if args.no_cuda else 'GPUs ' + ','.join(map(str, slots))})")
    while pending or running:
        while pending and free_slots:
            task = pending.pop(0)
            slot_id, slot_gpu = free_slots.pop(0)
            cmd = build_command(args, task, slot_gpu)
            run_dir = task_output_dir(base_output, task)
            log_dir = run_dir / "logs"
            log_dir.mkdir(parents=True, exist_ok=True)
            log_path = log_dir / "run_suite.log"
            command_text = " ".join(cmd)
            print(
                f"[launch] slot={slot_id} gpu={slot_gpu if slot_gpu is not None else 'cpu'} "
                f"{task['dataset']} {task['method']} seed={task['seed']}"
            )
            print(command_text)
            if args.dry_run:
                free_slots.append((slot_id, slot_gpu))
                continue

            log_handle = open(log_path, "a", encoding="utf-8")
            log_handle.write(f"\n\n===== run_suite launch {time.strftime('%Y-%m-%d %H:%M:%S')} =====\n")
            log_handle.write(command_text + "\n")
            log_handle.flush()
            env = os.environ.copy()
            if slot_gpu is not None:
                env["PLANTSPADE_ASSIGNED_GPU"] = str(slot_gpu)
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
            free_slots.append((item["slot_id"], item["gpu"]))
            status = "ok" if returncode == 0 else f"failed rc={returncode}"
            print(
                f"[done] {status} slot={item['slot_id']} gpu={item['gpu'] if item['gpu'] is not None else 'cpu'} "
                f"{task['dataset']} {task['method']} seed={task['seed']} log={item['log_path']}"
            )
            if returncode != 0:
                failures.append(
                    {
                        "dataset": task["dataset"],
                        "method": task["method"],
                        "seed": task["seed"],
                        "gpu": item["gpu"],
                        "returncode": returncode,
                        "log": str(item["log_path"]),
                    }
                )
                if not args.keep_going:
                    for other in still_running:
                        other["proc"].terminate()
                        other["log_handle"].close()
                    raise SystemExit(f"Stopping after failed run: {failures[-1]}")
        running = still_running

    return failures


def main():
    args = parse_args()
    datasets_cfg = load_yaml(args.datasets_config)
    main_cfg = load_yaml(args.main_config)
    datasets = select_datasets(datasets_cfg, args.datasets, args.groups)
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
