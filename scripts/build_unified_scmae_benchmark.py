#!/usr/bin/env python3
"""Build a disk-efficient, comparable view of the scMAE benchmark results."""

from __future__ import annotations

import argparse
import csv
import json
import os
import shutil
import statistics
from collections import defaultdict
from pathlib import Path
from typing import Any, Iterable


METRICS = ("acc", "ari", "nmi", "f1_macro", "runtime_seconds")
COMMON_FIELDS = (
    "benchmark_source",
    "record_scope",
    "source_run_dir",
    "dataset",
    "model",
    "seed",
    "status",
    "acc",
    "ari",
    "nmi",
    "f1_macro",
    "fmi",
    "v_measure",
    "homogeneity",
    "completeness",
    "runtime_seconds",
)


def hardlink_copy(source: Path, destination: Path, excluded: set[str] | None = None) -> None:
    """Copy a static result tree without duplicating file contents on this filesystem."""
    excluded = excluded or set()
    destination.mkdir(parents=True, exist_ok=True)
    for item in source.iterdir():
        if item.name in excluded:
            continue
        target = destination / item.name
        if item.is_dir():
            hardlink_copy(item, target, excluded)
        elif item.is_symlink():
            target.symlink_to(os.readlink(item))
        else:
            os.link(item, target)


def read_rows(path: Path, source: str, scope: str, stage7_only: bool = False) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            if stage7_only and row.get("stage") != "stage7":
                continue
            row["benchmark_source"] = source
            row["record_scope"] = scope
            row["source_run_dir"] = str(path.parent)
            rows.append(row)
    return rows


def numeric(row: dict[str, str], field: str) -> float | None:
    try:
        value = float(row.get(field, ""))
    except (TypeError, ValueError):
        return None
    return value if value == value else None


def summarise(rows: Iterable[dict[str, str]], fields: tuple[str, ...]) -> list[dict[str, Any]]:
    groups: dict[tuple[str, ...], list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        if row.get("status") == "success":
            groups[tuple(row.get(field, "") for field in fields)].append(row)

    summaries: list[dict[str, Any]] = []
    for key, group in sorted(groups.items()):
        record: dict[str, Any] = dict(zip(fields, key))
        record["n_success"] = len(group)
        for metric in METRICS:
            values = [value for row in group if (value := numeric(row, metric)) is not None]
            record[f"{metric}_mean"] = statistics.fmean(values) if values else ""
            record[f"{metric}_std"] = statistics.stdev(values) if len(values) > 1 else 0.0 if values else ""
        summaries.append(record)
    return summaries


def write_csv(path: Path, rows: list[dict[str, Any]], preferred: tuple[str, ...] = ()) -> None:
    all_fields = {field for row in rows for field in row}
    fieldnames = [field for field in preferred if field in all_fields]
    fieldnames.extend(sorted(all_fields - set(fieldnames)))
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--active", type=Path, required=True)
    parser.add_argument("--retired", type=Path, required=True)
    parser.add_argument("--tuned", type=Path, required=True)
    parser.add_argument("--history", type=Path, required=True)
    args = parser.parse_args()

    output = args.output
    if output.exists():
        raise SystemExit(f"Refusing to overwrite existing output directory: {output}")
    output.mkdir(parents=True)

    # Keep all static artifacts useful for comparison while excluding duplicated input data.
    hardlink_copy(args.active, output / "active_benchmark", {"converted_data"})
    hardlink_copy(args.retired, output / "retired_benchmark", {"converted_data"})
    hardlink_copy(args.history, output / "tuning_history_8datasets", {"converted_data"})
    for item in args.tuned.iterdir():
        if item.name in {"runs", "converted_data"}:
            continue
        target = output / "tuned_full16_final" / item.name
        target.parent.mkdir(parents=True, exist_ok=True)
        if item.is_dir():
            hardlink_copy(item, target)
        else:
            os.link(item, target)
    hardlink_copy(args.tuned / "runs" / "stage7", output / "tuned_full16_final" / "runs" / "stage7")

    combined_rows = []
    combined_rows.extend(read_rows(args.active / "all_runs_master.csv", "active_methods", "full_benchmark"))
    combined_rows.extend(read_rows(args.retired / "retired_all_runs_master.csv", "retired_methods", "full_benchmark"))
    combined_rows.extend(read_rows(args.tuned / "run_master.csv", "tuned_nm_rg_contrast", "final_best", stage7_only=True))
    write_csv(output / "combined_run_master.csv", combined_rows, COMMON_FIELDS)
    write_csv(output / "summary_by_source_model.csv", summarise(combined_rows, ("benchmark_source", "model")))
    write_csv(
        output / "summary_by_source_model_dataset.csv",
        summarise(combined_rows, ("benchmark_source", "model", "dataset")),
    )

    manifest = {
        "description": "Unified comparison view. Source trees are hard-link copies; deleting a source tree does not delete these artifacts.",
        "sources": {
            "active_methods": str(args.active),
            "retired_methods": str(args.retired),
            "tuned_nm_rg_contrast": str(args.tuned),
            "tuning_history_8datasets": str(args.history),
        },
        "combined_master_includes": {
            "active_methods": "all rows",
            "retired_methods": "all rows",
            "tuned_nm_rg_contrast": "stage7 final-best rows only",
        },
        "excluded_from_artifact_copies": ["converted_data"],
    }
    (output / "comparison_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    (output / "README.md").write_text(
        "# Unified scMAE Comparison Benchmark\n\n"
        "`combined_run_master.csv` is the comparison table. `tuned_nm_rg_contrast` rows are restricted to stage7 final-best trials; active and retired sources retain their full benchmark rows. "
        "The source artifact trees are hard-link copies, excluding converted input data.\n"
    )


if __name__ == "__main__":
    main()
