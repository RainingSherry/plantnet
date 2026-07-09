#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]
BASE_DIR = Path(__file__).resolve().parent
DEFAULT_OLD_ROOT = PROJECT_ROOT / "result" / "scmae_all_methods_20260705_full"
DEFAULT_NEW_RUN_ROOT = BASE_DIR / "runs"
DEFAULT_SEEDS = "42,3047,3407"
METRICS = ["ari", "nmi", "acc"]


def read_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False, allow_nan=True) + "\n", encoding="utf-8")


def parse_csv(value: str | None) -> list[str] | None:
    if value is None or value in {"", "all", "auto"}:
        return None
    return [x.strip() for x in str(value).split(",") if x.strip()]


def discover_datasets(old_root: Path) -> list[str]:
    catalog = old_root / "dataset_catalog.json"
    if catalog.exists():
        return [str(item["name"]) for item in read_json(catalog)]
    return sorted(p.name for p in old_root.iterdir() if (p / "pca_kmeans_known_k").exists())


def load_old_master(old_root: Path) -> pd.DataFrame | None:
    path = old_root / "all_runs_master.csv"
    if not path.exists():
        return None
    df = pd.read_csv(path)
    return df[df["model"].astype(str) == "pca_kmeans_known_k"].copy()


def old_metrics(old_root: Path, old_master: pd.DataFrame | None, dataset: str, seed: int) -> dict[str, float] | None:
    metrics_path = old_root / dataset / "pca_kmeans_known_k" / f"seed_{seed}" / "metrics.json"
    if metrics_path.exists():
        payload = read_json(metrics_path)
        fixed = payload.get("kmeans_known_k", payload)
        return {key: float(fixed[key]) for key in METRICS}
    if old_master is not None:
        row = old_master[(old_master["dataset"].astype(str) == dataset) & (old_master["seed"].astype(int) == int(seed))]
        if not row.empty:
            rec = row.iloc[0]
            return {key: float(rec[key]) for key in METRICS}
    return None


def new_metrics(new_run_root: Path, dataset: str, seed: int) -> dict[str, float] | None:
    metrics_path = new_run_root / dataset / f"seed_{seed}" / "pca_kmeans" / "metrics.json"
    if not metrics_path.exists():
        return None
    payload = read_json(metrics_path)
    return {key: float(payload[key]) for key in METRICS}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--old_root", type=Path, default=DEFAULT_OLD_ROOT)
    parser.add_argument("--new_run_root", type=Path, default=DEFAULT_NEW_RUN_ROOT)
    parser.add_argument("--out_dir", type=Path, default=BASE_DIR / "analysis_full")
    parser.add_argument("--datasets", default="all")
    parser.add_argument("--seeds", default=DEFAULT_SEEDS)
    parser.add_argument("--tolerance", type=float, default=1e-10)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    datasets = parse_csv(args.datasets) or discover_datasets(args.old_root)
    seeds = [int(x) for x in parse_csv(args.seeds) or parse_csv(DEFAULT_SEEDS)]
    old_master = load_old_master(args.old_root)
    rows: list[dict[str, Any]] = []
    failures: list[dict[str, Any]] = []

    for dataset in datasets:
        for seed in seeds:
            old = old_metrics(args.old_root, old_master, dataset, seed)
            new = new_metrics(args.new_run_root, dataset, seed)
            row: dict[str, Any] = {"dataset": dataset, "seed": int(seed)}
            if old is None:
                row["status"] = "missing_old"
                failures.append(row.copy())
            elif new is None:
                row["status"] = "missing_new"
                failures.append(row.copy())
            else:
                row["status"] = "ok"
                for key in METRICS:
                    diff = float(new[key] - old[key])
                    row[f"old_{key}"] = old[key]
                    row[f"new_{key}"] = new[key]
                    row[f"{key}_diff"] = diff
                    if not np.isfinite(diff) or abs(diff) > args.tolerance:
                        row["status"] = "diff_exceeds_tolerance"
                if row["status"] != "ok":
                    failures.append(row.copy())
            rows.append(row)

    out_dir = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame(rows)
    csv_path = out_dir / "pca_compatibility_report.csv"
    json_path = out_dir / "pca_compatibility_report.json"
    df.to_csv(csv_path, index=False)
    diff_cols = [f"{key}_diff" for key in METRICS if f"{key}_diff" in df.columns]
    summary = {
        "old_root": str(args.old_root),
        "new_run_root": str(args.new_run_root),
        "n_expected": int(len(datasets) * len(seeds)),
        "n_checked": int(len(df)),
        "tolerance": float(args.tolerance),
        "passed": not failures,
        "n_failures": int(len(failures)),
        "max_abs_diff": {col: float(df[col].abs().max()) for col in diff_cols},
        "mean_diff": {col: float(df[col].mean()) for col in diff_cols},
        "mean_abs_diff": {col: float(df[col].abs().mean()) for col in diff_cols},
        "failures": failures[:50],
        "csv": str(csv_path),
    }
    write_json(json_path, summary)
    print(json.dumps(summary, indent=2, ensure_ascii=False, allow_nan=True))
    return 0 if summary["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
