#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path
from statistics import mean


ROOT = Path(__file__).resolve().parents[1]
SWEEP_ORDER = {"gate_max": 0, "neighbor_k": 1, "pseudo_weight": 2}
RAW_COLUMNS = [
    "sweep",
    "dataset",
    "label",
    "value",
    "seed",
    "gate_max",
    "neighbor_k",
    "pseudo_weight",
    "ari",
    "acc",
    "nmi",
    "mean_node_gate",
    "effective_neighbor_count",
    "max_edge_weight_p95",
    "path",
]
AGG_COLUMNS = ["sweep", "value", "n_datasets", "mean_ari", "min_ari", "max_ari", "datasets"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Summarize RG_NeighborMix_scMAE Phase 2 sweep outputs.")
    parser.add_argument("--out_dir", default="results/formal/rg_phase2_sensitivity_e80")
    parser.add_argument("--raw_name", default="rg_phase2_all_sweeps_raw.csv")
    parser.add_argument("--aggregate_name", default="rg_phase2_all_sweeps_aggregate.csv")
    return parser.parse_args()


def read_json(path: Path) -> dict:
    if not path.exists():
        return {}
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def first_csv_row(path: Path) -> dict[str, str]:
    with path.open("r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        try:
            return next(reader)
        except StopIteration as exc:
            raise ValueError(f"Empty CSV: {path}") from exc


def to_float(value) -> float | None:
    if value is None or value == "":
        return None
    return float(value)


def rel_path(path: Path) -> str:
    try:
        return str(path.relative_to(ROOT))
    except ValueError:
        return str(path)


def collect_rows(out_root: Path) -> list[dict]:
    rows = []
    for eval_path in sorted(out_root.glob("*/*/*/seed*/eval_fixed.csv")):
        run_dir = eval_path.parent
        sweep, dataset, label, seed_label = eval_path.relative_to(out_root).parts[:4]
        if sweep not in SWEEP_ORDER:
            continue
        seed = int(seed_label.removeprefix("seed"))
        eval_row = first_csv_row(eval_path)
        args = read_json(run_dir / "args.json")
        gate_summary = read_json(run_dir / "gate_summary.json")
        edge_summary = read_json(run_dir / "edge_weight_summary.json")

        gate_max = to_float(args.get("gate_max", eval_row.get("gate_max")))
        neighbor_k = args.get("neighbor_k")
        if neighbor_k is not None:
            neighbor_k = int(neighbor_k)
        pseudo_weight = to_float(args.get("pseudo_weight", eval_row.get("pseudo_weight")))

        if sweep == "gate_max":
            value = gate_max
        elif sweep == "neighbor_k":
            value = neighbor_k
        else:
            value = pseudo_weight

        rows.append(
            {
                "sweep": sweep,
                "dataset": dataset,
                "label": label,
                "value": value,
                "seed": seed,
                "gate_max": gate_max,
                "neighbor_k": neighbor_k,
                "pseudo_weight": pseudo_weight,
                "ari": to_float(eval_row.get("ari")),
                "acc": to_float(eval_row.get("acc")),
                "nmi": to_float(eval_row.get("nmi")),
                "mean_node_gate": to_float(gate_summary.get("mean_node_gate")),
                "effective_neighbor_count": to_float(edge_summary.get("effective_neighbor_count")),
                "max_edge_weight_p95": to_float(edge_summary.get("max_edge_weight_p95")),
                "path": rel_path(run_dir),
            }
        )

    rows.sort(key=lambda r: (SWEEP_ORDER[r["sweep"]], float(r["value"]), r["dataset"], r["seed"]))
    return rows


def aggregate_rows(rows: list[dict]) -> list[dict]:
    by_dataset: dict[tuple[str, float, str], list[float]] = defaultdict(list)
    for row in rows:
        ari = row["ari"]
        if ari is None:
            continue
        by_dataset[(row["sweep"], float(row["value"]), row["dataset"])].append(float(ari))

    grouped: dict[tuple[str, float], list[tuple[str, float]]] = defaultdict(list)
    for (sweep, value, dataset), values in by_dataset.items():
        grouped[(sweep, value)].append((dataset, mean(values)))

    aggregate = []
    for (sweep, value), dataset_scores in sorted(grouped.items(), key=lambda x: (SWEEP_ORDER[x[0][0]], x[0][1])):
        dataset_scores.sort(key=lambda item: item[0])
        scores = [score for _, score in dataset_scores]
        aggregate.append(
            {
                "sweep": sweep,
                "value": value,
                "n_datasets": len(dataset_scores),
                "mean_ari": mean(scores),
                "min_ari": min(scores),
                "max_ari": max(scores),
                "datasets": ";".join(f"{dataset}:{score:.4f}" for dataset, score in dataset_scores),
            }
        )
    return aggregate


def write_csv(path: Path, columns: list[str], rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    args = parse_args()
    out_root = Path(args.out_dir)
    if not out_root.is_absolute():
        out_root = ROOT / out_root
    rows = collect_rows(out_root)
    if not rows:
        raise SystemExit(f"No eval_fixed.csv files found under {out_root}")
    aggregate = aggregate_rows(rows)
    raw_path = out_root / args.raw_name
    aggregate_path = out_root / args.aggregate_name
    write_csv(raw_path, RAW_COLUMNS, rows)
    write_csv(aggregate_path, AGG_COLUMNS, aggregate)
    datasets = sorted({row["dataset"] for row in rows})
    print(f"Wrote {len(rows)} raw rows across {len(datasets)} datasets: {', '.join(datasets)}")
    print(f"Wrote {len(aggregate)} aggregate rows")
    print(f"raw={rel_path(raw_path)}")
    print(f"aggregate={rel_path(aggregate_path)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
