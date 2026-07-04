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
    values = [v for v in values if v is not None]
    if not values:
        return float("nan"), float("nan")
    if len(values) == 1:
        return float(values[0]), 0.0
    return float(statistics.mean(values)), float(statistics.stdev(values))


def fmt(values: list[float], digits: int = 4) -> str:
    mean, sd = mean_sd(values)
    return f"{mean:.{digits}f} +/- {sd:.{digits}f}"


def arm_name(run: str) -> str:
    # strip trailing _seedNN
    return run.rsplit("_seed", 1)[0]


def main() -> int:
    root = Path(__file__).resolve().parent
    rows = []
    for summary_path in sorted((root / "runs").glob("*/summary.json")):
        if summary_path.parent.name.startswith("smoke"):
            continue
        s = load_summary(summary_path)
        m = s.get("fixed_metrics", {}).get("kmeans_known_k", {})
        std = s.get("std_profile", {})
        aligned = s.get("cluster_aligned", {})
        rows.append({
            "run": summary_path.parent.name,
            "arm": arm_name(summary_path.parent.name),
            "seed": s.get("seed"),
            "program_mode": s.get("program_mode"),
            "split_mode": s.get("split_mode"),
            "program_weight": s.get("program_weight"),
            "type_dim": s.get("type_dim"),
            "ari": m.get("ari"),
            "nmi": m.get("nmi"),
            "acc": m.get("acc"),
            "effective_dim_pr": std.get("effective_dim_pr"),
            "cluster_aligned_eff_dim": aligned.get("cluster_aligned_eff_dim"),
            "std_min": std.get("std_min"),
            "dims_std_gt_1p0": std.get("dims_std_gt_1p0"),
            "program_r2": s.get("program_r2"),
            "final_base_loss": s.get("final_base_loss"),
            "final_program_loss": s.get("final_program_loss"),
            "runtime_seconds": s.get("runtime_seconds"),
        })

    out_csv = root / "summary.csv"
    if rows:
        with open(out_csv, "w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
            writer.writeheader()
            writer.writerows(rows)

    out_md = root / "SUMMARY.md"
    with open(out_md, "w", encoding="utf-8") as handle:
        handle.write("# Gene-program bottleneck ablation summary (Macosko)\n\n")
        if not rows:
            handle.write("No completed runs found under `runs/*/summary.json`.\n")
            print(f"Wrote {out_md} (empty)")
            return 0

        groups = defaultdict(list)
        for row in rows:
            groups[row["arm"]].append(row)

        # arm ordering for readability
        order = ["a0_baseline", "a1_prog_w02", "a1_prog_w10", "a1_shuffle_w10", "a2_fixed_w10", "a2_extra_w10"]
        keys = [k for k in order if k in groups] + [k for k in sorted(groups) if k not in order]

        handle.write("## Arm-level aggregate\n\n")
        handle.write(
            "| arm | program | split | pw | n | ARI mean+/-sd | NMI mean+/-sd | "
            "eff_dim | aligned_eff_dim | prog_R2 |\n"
        )
        handle.write("|---|---|---|---:|---:|---:|---:|---:|---:|---:|\n")
        for k in keys:
            g = groups[k]
            r0 = g[0]
            handle.write(
                f"| {k} | {r0['program_mode']} | {r0['split_mode']} | {r0['program_weight']} | {len(g)} | "
                f"{fmt([r['ari'] for r in g])} | {fmt([r['nmi'] for r in g])} | "
                f"{fmt([r['effective_dim_pr'] for r in g], 1)} | "
                f"{fmt([r['cluster_aligned_eff_dim'] for r in g], 2)} | "
                f"{fmt([r['program_r2'] for r in g], 3)} |\n"
            )

        # decision helper: deltas vs a0_baseline
        if "a0_baseline" in groups:
            a0_ari, _ = mean_sd([r["ari"] for r in groups["a0_baseline"]])
            handle.write("\n## Delta vs a0_baseline (ARI)\n\n")
            handle.write("| arm | ARI mean | delta | verdict (|delta|>=0.02) |\n|---|---:|---:|---|\n")
            for k in keys:
                mean, _ = mean_sd([r["ari"] for r in groups[k]])
                d = mean - a0_ari
                verdict = "up" if d >= 0.02 else ("down" if d <= -0.02 else "tie")
                handle.write(f"| {k} | {mean:.4f} | {d:+.4f} | {verdict} |\n")

        handle.write("\n## Per-run details\n\n")
        handle.write("| run | ARI | NMI | eff_dim | aligned_eff_dim | std_min | prog_R2 | base_loss | prog_loss |\n")
        handle.write("|---|---:|---:|---:|---:|---:|---:|---:|---:|\n")
        for row in rows:
            def g(x, d=4):
                return f"{float(x):.{d}f}" if x is not None else "n/a"
            handle.write(
                f"| {row['run']} | {g(row['ari'])} | {g(row['nmi'])} | {g(row['effective_dim_pr'],1)} | "
                f"{g(row['cluster_aligned_eff_dim'],2)} | {g(row['std_min'],3)} | {g(row['program_r2'],3)} | "
                f"{g(row['final_base_loss'],4)} | {g(row['final_program_loss'],4)} |\n"
            )

        handle.write("\n## Decision rule\n\n")
        handle.write("- **a1 > a0** (ARI +>=0.02, stable): program aux-regularizer helps -> proceed to posterior fusion.\n")
        handle.write("- **a1 ~= a0**: program is neutral decoration -> stop.\n")
        handle.write("- **a1 < a0**: expression-derived program semantics hurt clustering -> stop.\n")
        handle.write("- **a1_shuffle ~= a1**: gain (if any) is generic aux-loss regularization, NOT program semantics.\n")
        handle.write("- **a2_fixed < a1**: shrinking the clustering subspace conflicts with the std-floor mechanism.\n")
        handle.write("- **a2_extra vs a2_fixed**: isolates 'split' (bad) from 'fewer clustering dims' (confound).\n")
        handle.write("- **aligned_eff_dim**: participation ratio of between-class scatter; the key semantic-alignment metric.\n")
    print(f"Wrote {out_csv}")
    print(f"Wrote {out_md}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
