#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_SOURCE = PROJECT_ROOT / "result" / "scmae_all_methods_20260705_full"


def read_json(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


def parse_csv(value: str) -> list[str]:
    return [x.strip() for x in str(value).split(",") if x.strip()]


def load_dataset_names(source_root: Path, datasets: str) -> list[str]:
    if datasets != "all":
        return parse_csv(datasets)
    catalog = read_json(source_root / "dataset_catalog.json")
    return [item["name"] for item in catalog]


def status_for(run_root: Path, dataset: str, seed: int) -> dict:
    run_dir = run_root / dataset / f"seed_{seed}"
    status_path = run_dir / "status.json"
    if not status_path.exists():
        return {"dataset": dataset, "seed": seed, "status": "missing", "runtime_seconds": 0.0, "run_dir": str(run_dir), "error": ""}
    status = read_json(status_path)
    return {
        "dataset": dataset,
        "seed": seed,
        "status": status.get("status", "unknown"),
        "runtime_seconds": float(status.get("runtime_seconds", 0.0)),
        "run_dir": str(run_dir),
        "error": str(status.get("error", ""))[:240],
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--run_root", type=Path, default=Path(__file__).resolve().parent / "runs")
    parser.add_argument("--source_root", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--datasets", default="all")
    parser.add_argument("--seeds", default="42,3047,3407")
    parser.add_argument("--csv", type=Path, default=None)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    datasets = load_dataset_names(args.source_root, args.datasets)
    seeds = [int(x) for x in parse_csv(args.seeds)]
    rows = [status_for(args.run_root, dataset, seed) for dataset in datasets for seed in seeds]
    df = pd.DataFrame(rows)
    if args.csv:
        args.csv.parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(args.csv, index=False)
    counts = df["status"].value_counts(dropna=False).to_dict()
    total_runtime = float(df.loc[df["status"] == "success", "runtime_seconds"].sum())
    print(f"Run root: {args.run_root}")
    print(f"Expected tasks: {len(df)}")
    print(f"Status counts: {counts}")
    print(f"Completed runtime seconds: {total_runtime:.1f}")
    print()
    print(
        df.sort_values(["status", "dataset", "seed"])[
            ["dataset", "seed", "status", "runtime_seconds", "error"]
        ].to_string(index=False)
    )
    return 1 if (df["status"] == "failed").any() else 0


if __name__ == "__main__":
    raise SystemExit(main())
