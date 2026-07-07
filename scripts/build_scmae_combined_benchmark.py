#!/usr/bin/env python3
"""Build combined scMAE benchmark tables from active and retired result suites."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


METRICS = ("ari", "nmi", "acc")


def load_summary(path: Path, source: str) -> pd.DataFrame:
    df = pd.read_csv(path)
    required = ["model", "dataset", "n_success"] + [f"{m}_mean" for m in METRICS]
    missing = [col for col in required if col not in df.columns]
    if missing:
        raise ValueError(f"{path} is missing required columns: {missing}")

    out = df[required].copy()
    out.insert(0, "source", source)
    out = out.rename(columns={f"{m}_mean": m.upper() for m in METRICS})
    return out


def format_value(value: float, is_best: bool) -> str:
    if pd.isna(value):
        return "-"
    text = f"{value:.4f}"
    return f"**{text}**" if is_best else text


def to_markdown_table(df: pd.DataFrame, best_lookup: dict[tuple[str, str], float]) -> str:
    display = df.copy()
    for metric in ("ARI", "NMI", "ACC"):
        display[metric] = [
            format_value(
                value,
                pd.notna(value)
                and np.isclose(value, best_lookup.get((dataset, metric), np.nan), rtol=0, atol=1e-12),
            )
            for dataset, value in zip(display["dataset"], display[metric])
        ]
    display["n_success"] = display["n_success"].apply(lambda v: "-" if pd.isna(v) else str(int(v)))

    headers = ["source", "model", "dataset", "n_success", "ARI", "NMI", "ACC"]
    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join(["---"] * len(headers)) + " |",
    ]
    for row in display[headers].itertuples(index=False, name=None):
        lines.append("| " + " | ".join(str(value) for value in row) + " |")
    return "\n".join(lines) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--all-summary",
        type=Path,
        default=Path("result/scmae_all_methods_20260705_full/summary_by_model_dataset.csv"),
    )
    parser.add_argument(
        "--retired-summary",
        type=Path,
        default=Path("result/scmae_retired_methods_20260706_full/summary_by_model_dataset.csv"),
    )
    parser.add_argument(
        "--out-prefix",
        type=Path,
        default=Path("result/scmae_combined_benchmark_ari_nmi_acc_20260707"),
    )
    args = parser.parse_args()

    combined = pd.concat(
        [
            load_summary(args.all_summary, "scmae_all_methods_20260705_full"),
            load_summary(args.retired_summary, "scmae_retired_methods_20260706_full"),
        ],
        ignore_index=True,
    )
    duplicates = combined[combined.duplicated(["model", "dataset"], keep=False)]
    if not duplicates.empty:
        dup_pairs = duplicates[["model", "dataset"]].drop_duplicates().to_dict("records")
        raise ValueError(f"Duplicate model/dataset rows found: {dup_pairs[:20]}")

    datasets = sorted(combined["dataset"].unique())
    model_order = (
        combined[["source", "model"]]
        .drop_duplicates()
        .sort_values(["source", "model"], kind="stable")
        .reset_index(drop=True)
    )
    grid = (
        model_order.assign(key=1)
        .merge(pd.DataFrame({"dataset": datasets, "key": 1}), on="key")
        .drop(columns="key")
    )
    complete = grid.merge(combined, on=["source", "model", "dataset"], how="left")
    complete = complete[["source", "model", "dataset", "n_success", "ARI", "NMI", "ACC"]]

    best_lookup: dict[tuple[str, str], float] = {}
    for dataset, sub in combined.groupby("dataset"):
        for metric in ("ARI", "NMI", "ACC"):
            best_lookup[(dataset, metric)] = sub[metric].max(skipna=True)

    csv_path = args.out_prefix.with_suffix(".csv")
    md_path = args.out_prefix.with_suffix(".md")
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    complete.to_csv(csv_path, index=False)

    note = (
        "# Combined scMAE Benchmark (ARI/NMI/ACC)\n\n"
        "Values are mean scores from `summary_by_model_dataset.csv`. "
        "Bold values mark the best method for each dataset and metric. "
        "Missing model/dataset combinations are shown as `-`.\n\n"
    )
    md_path.write_text(note + to_markdown_table(complete, best_lookup), encoding="utf-8")

    observed_rows = int(combined.shape[0])
    total_rows = int(complete.shape[0])
    missing_rows = int(complete["ARI"].isna().sum())
    print(f"Wrote {csv_path} ({total_rows} rows, {missing_rows} missing combinations)")
    print(f"Wrote {md_path} ({observed_rows} observed model/dataset rows)")


if __name__ == "__main__":
    main()
