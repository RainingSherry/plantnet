#!/usr/bin/env python3
from __future__ import annotations

import csv
import json
import math
import statistics
from collections import defaultdict
from pathlib import Path


def load_json(path: Path) -> dict:
    with open(path, "r", encoding="utf-8") as handle:
        return json.load(handle)


def mean_sd(values: list[float]) -> tuple[float, float]:
    values = [float(value) for value in values if value is not None and math.isfinite(float(value))]
    if not values:
        return float("nan"), float("nan")
    if len(values) == 1:
        return float(values[0]), 0.0
    return float(statistics.mean(values)), float(statistics.stdev(values))


def fmt(values: list[float], digits: int = 4) -> str:
    mean, sd = mean_sd(values)
    if not math.isfinite(mean):
        return "NA"
    return f"{mean:.{digits}f} +/- {sd:.{digits}f}"


def metric(summary: dict, name: str, key: str) -> float:
    return float(summary.get("fixed_metrics", {}).get(name, {}).get(key, float("nan")))


def main() -> int:
    root = Path(__file__).resolve().parent
    rows = []
    for summary_path in sorted((root / "runs").glob("*/summary.json")):
        if summary_path.parent.name.startswith("smoke"):
            continue
        summary = load_json(summary_path)
        std = summary.get("std_profile", {})
        rows.append(
            {
                "run": summary_path.parent.name,
                "dataset": summary.get("dataset"),
                "seed": summary.get("seed"),
                "method": summary.get("method"),
                "assignment_mode": summary.get("assignment_mode"),
                "variance_weight": summary.get("variance_weight"),
                "latent_dim": summary.get("latent_dim"),
                "kmeans_ari": metric(summary, "kmeans_known_k", "ari"),
                "kmeans_nmi": metric(summary, "kmeans_known_k", "nmi"),
                "direct_ari": metric(summary, "direct_prototype_argmax", "ari"),
                "direct_nmi": metric(summary, "direct_prototype_argmax", "nmi"),
                "pca_ari": metric(summary, "pca_kmeans_known_k", "ari"),
                "pca_nmi": metric(summary, "pca_kmeans_known_k", "nmi"),
                "effective_dim_pr": std.get("effective_dim_pr"),
                "std_min": std.get("std_min"),
                "std_median": std.get("std_median"),
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
        handle.write("# Objective-swap scMAE experiment summary\n\n")
        if not rows:
            handle.write("No completed non-smoke runs found under `runs/*/summary.json`.\n")
        else:
            groups = defaultdict(list)
            for row in rows:
                groups[(row["dataset"], row["assignment_mode"], row["variance_weight"], row["latent_dim"])].append(row)
            handle.write("## Aggregate\n\n")
            handle.write(
                "| dataset | assignment | varw | latent | n | KMeans ARI | direct ARI | PCA ARI | eff_dim | std_min |\n"
            )
            handle.write("|---|---|---:|---:|---:|---:|---:|---:|---:|---:|\n")
            for key in sorted(groups, key=lambda item: (str(item[0]), str(item[1]), float(item[2]), int(item[3]))):
                group = groups[key]
                dataset, assignment, varw, latent = key
                handle.write(
                    f"| {dataset} | {assignment} | {varw} | {latent} | {len(group)} | "
                    f"{fmt([float(row['kmeans_ari']) for row in group])} | "
                    f"{fmt([float(row['direct_ari']) for row in group])} | "
                    f"{fmt([float(row['pca_ari']) for row in group])} | "
                    f"{fmt([float(row['effective_dim_pr']) for row in group], 1)} | "
                    f"{fmt([float(row['std_min']) for row in group], 3)} |\n"
                )
            handle.write("\n## Per-run details\n\n")
            handle.write(
                "| run | assignment | KMeans ARI | direct ARI | PCA ARI | eff_dim | std_min | dims_std>1 |\n"
            )
            handle.write("|---|---:|---:|---:|---:|---:|---:|---:|\n")
            for row in rows:
                handle.write(
                    f"| {row['run']} | {row['assignment_mode']} | {float(row['kmeans_ari']):.4f} | "
                    f"{float(row['direct_ari']):.4f} | {float(row['pca_ari']):.4f} | "
                    f"{float(row['effective_dim_pr']):.1f} | {float(row['std_min']):.3f} | "
                    f"{row['dims_std_gt_1p0']} |\n"
                )

    print(f"Wrote {out_csv}")
    print(f"Wrote {out_md}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
