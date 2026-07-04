#!/usr/bin/env python3
from __future__ import annotations

import csv
import json
import statistics
from collections import defaultdict
from pathlib import Path


def load_summary(path: Path) -> dict:
    with open(path, "r", encoding="utf-8") as handle:
        return json.load(handle)


def mean_sd(values: list[float]) -> tuple[float, float]:
    if not values:
        return float("nan"), float("nan")
    if len(values) == 1:
        return float(values[0]), 0.0
    return float(statistics.mean(values)), float(statistics.stdev(values))


def fmt_mean_sd(values: list[float], digits: int = 4) -> str:
    mean, sd = mean_sd(values)
    return f"{mean:.{digits}f} +/- {sd:.{digits}f}"


def main() -> int:
    root = Path(__file__).resolve().parent
    rows = []
    for summary_path in sorted((root / "runs").glob("*/summary.json")):
        if summary_path.parent.name.startswith("smoke"):
            continue
        summary = load_summary(summary_path)
        metrics = summary.get("fixed_metrics", {}).get("kmeans_known_k", {})
        std = summary.get("std_profile", {})
        rows.append(
            {
                "run": summary_path.parent.name,
                "dataset": summary.get("dataset"),
                "seed": summary.get("seed"),
                "mix_mode": summary.get("mix_mode"),
                "variance_weight": summary.get("variance_weight"),
                "denoise_weight": summary.get("denoise_weight"),
                "ari": metrics.get("ari"),
                "nmi": metrics.get("nmi"),
                "acc": metrics.get("acc"),
                "f1_macro": metrics.get("f1_macro"),
                "effective_dim_pr": std.get("effective_dim_pr"),
                "std_min": std.get("std_min"),
                "std_median": std.get("std_median"),
                "std_max": std.get("std_max"),
                "dims_std_gt_1p0": std.get("dims_std_gt_1p0"),
                "cluster_mass_min": summary.get("cluster_mass_min"),
                "cluster_mass_max": summary.get("cluster_mass_max"),
                "runtime_seconds": summary.get("runtime_seconds"),
            }
        )

    out_csv = root / "summary.csv"
    if rows:
        with open(out_csv, "w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
            writer.writeheader()
            writer.writerows(rows)
    out_md = root / "SUMMARY.md"
    with open(out_md, "w", encoding="utf-8") as handle:
        handle.write("# NeighborMix x std-floor ablation summary\n\n")
        if not rows:
            handle.write("No completed runs found under `runs/*/summary.json`.\n")
        else:
            groups = defaultdict(list)
            for row in rows:
                groups[(row["mix_mode"], row["variance_weight"])].append(row)

            handle.write("## Arm-level aggregate\n\n")
            handle.write(
                "| mix | varw | n | ARI mean+/-sd | NMI mean+/-sd | "
                "eff_dim mean+/-sd | std_min mean+/-sd | dims_std>1 mean+/-sd |\n"
            )
            handle.write("|---|---:|---:|---:|---:|---:|---:|---:|\n")
            for key in sorted(groups, key=lambda item: (str(item[0]), float(item[1]))):
                group_rows = groups[key]
                mix_mode, variance_weight = key
                handle.write(
                    f"| {mix_mode} | {variance_weight} | {len(group_rows)} | "
                    f"{fmt_mean_sd([float(row['ari']) for row in group_rows])} | "
                    f"{fmt_mean_sd([float(row['nmi']) for row in group_rows])} | "
                    f"{fmt_mean_sd([float(row['effective_dim_pr']) for row in group_rows], 1)} | "
                    f"{fmt_mean_sd([float(row['std_min']) for row in group_rows], 3)} | "
                    f"{fmt_mean_sd([float(row['dims_std_gt_1p0']) for row in group_rows], 1)} |\n"
                )

            handle.write("\n## Per-run details\n\n")
            handle.write("| run | mix | varw | ARI | NMI | eff_dim | std_min | std_med | dims_std>1 |\n")
            handle.write("|---|---:|---:|---:|---:|---:|---:|---:|---:|\n")
            for row in rows:
                handle.write(
                    f"| {row['run']} | {row['mix_mode']} | {row['variance_weight']} | "
                    f"{float(row['ari']):.4f} | {float(row['nmi']):.4f} | "
                    f"{float(row['effective_dim_pr']):.1f} | {float(row['std_min']):.3f} | "
                    f"{float(row['std_median']):.3f} | {row['dims_std_gt_1p0']} |\n"
                )
            handle.write("\nInterpretation guide:\n\n")
            handle.write("- `none + varw=0`: DEC-only control.\n")
            handle.write("- `none + varw=0.02`: std-floor intervention.\n")
            handle.write("- `neighbor + varw=0`: NeighborMix auxiliary denoising without std-floor.\n")
            handle.write("- `neighbor + varw=0.02`: test whether NeighborMix and std-floor are complementary.\n")
    print(f"Wrote {out_csv}")
    print(f"Wrote {out_md}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
