#!/usr/bin/env python3
"""聚合 conditional-shuffle 消融结果：按 (corruption, dataset) 汇总 ARI/NMI mean±sd，

并计算两组 delta：
  - vs zero      : 相对「复现赢家」零填充基线
  - vs swap_global(S0): 相对原始 scMAE 全局 swap 基线（Phase 2 的真正判据）
同时输出每个 arm 的平均 donor-pool 大小与实际改动率（监控 swap 是否退化）。

判据：S1/S2/S3 相对 S0 的 ARI +>=0.02 且多种子稳定，才算「从 scMAE 内部改对了」。
"""
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
        pool = s.get("donor_pool", {})
        rows.append({
            "run": summary_path.parent.name,
            "dataset": s.get("dataset"),
            "corruption": s.get("corruption"),
            "seed": s.get("seed"),
            "ari": m.get("ari"),
            "nmi": m.get("nmi"),
            "acc": m.get("acc"),
            "effective_dim_pr": std.get("effective_dim_pr"),
            "cluster_aligned_eff_dim": aligned.get("cluster_aligned_eff_dim"),
            "mean_effective_change_rate": s.get("mean_effective_change_rate"),
            "donor_pool_mean": pool.get("pool_mean"),
            "donor_pool_min": pool.get("pool_min"),
            "n_bins": pool.get("n_bins"),
            "final_scmae_loss": s.get("final_scmae_loss"),
            "runtime_seconds": s.get("runtime_seconds"),
        })

    out_csv = root / "summary.csv"
    if rows:
        with open(out_csv, "w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
            writer.writeheader()
            writer.writerows(rows)

    out_md = root / "SUMMARY.md"
    _write_md(out_md, rows)
    print(f"Wrote {out_csv}")
    print(f"Wrote {out_md}")
    return 0
ARM_ORDER = ["zero", "swap_global", "swap_lib", "swap_ndet", "swap_zerolib"]


def _write_md(out_md: Path, rows: list[dict]) -> None:
    with open(out_md, "w", encoding="utf-8") as handle:
        handle.write("# Conditional / nuisance-matched shuffle ablation summary\n\n")
        if not rows:
            handle.write("No completed runs found under `runs/*/summary.json`.\n")
            print(f"Wrote {out_md} (empty)")
            return

        # group by (dataset, corruption)
        groups: dict[tuple[str, str], list[dict]] = defaultdict(list)
        datasets: list[str] = []
        for row in rows:
            groups[(row["dataset"], row["corruption"])].append(row)
            if row["dataset"] not in datasets:
                datasets.append(row["dataset"])

        handle.write("## Arm-level aggregate (per dataset)\n\n")
        handle.write(
            "| dataset | corruption | n | ARI mean+/-sd | NMI mean+/-sd | eff_dim | "
            "aligned_eff_dim | eff_change | donor_pool_mean | n_bins |\n"
        )
        handle.write("|---|---|---:|---:|---:|---:|---:|---:|---:|---:|\n")
        for ds in datasets:
            arms = [a for a in ARM_ORDER if (ds, a) in groups]
            arms += [a for (d, a) in groups if d == ds and a not in ARM_ORDER]
            for arm in arms:
                g = groups[(ds, arm)]
                r0 = g[0]
                handle.write(
                    f"| {ds} | {arm} | {len(g)} | {fmt([r['ari'] for r in g])} | "
                    f"{fmt([r['nmi'] for r in g])} | {fmt([r['effective_dim_pr'] for r in g], 1)} | "
                    f"{fmt([r['cluster_aligned_eff_dim'] for r in g], 2)} | "
                    f"{fmt([r['mean_effective_change_rate'] for r in g], 3)} | "
                    f"{r0['donor_pool_mean']:.0f} | {r0['n_bins']} |\n"
                )

        # delta tables: vs zero and vs swap_global (S0)
        for ref in ("zero", "swap_global"):
            handle.write(f"\n## Delta vs {ref} (ARI)\n\n")
            handle.write("| dataset | corruption | ARI mean | delta | verdict (delta>=+0.02) |\n")
            handle.write("|---|---|---:|---:|---|\n")
            for ds in datasets:
                if (ds, ref) not in groups:
                    continue
                ref_ari, _ = mean_sd([r["ari"] for r in groups[(ds, ref)]])
                arms = [a for a in ARM_ORDER if (ds, a) in groups]
                for arm in arms:
                    mean, _ = mean_sd([r["ari"] for r in groups[(ds, arm)]])
                    d = mean - ref_ari
                    verdict = "up" if d >= 0.02 else ("down" if d <= -0.02 else "tie")
                    handle.write(f"| {ds} | {arm} | {mean:.4f} | {d:+.4f} | {verdict} |\n")

        handle.write("\n## Per-run details\n\n")
        handle.write(
            "| run | dataset | corruption | seed | ARI | NMI | eff_dim | aligned | "
            "eff_change | pool_mean | pool_min |\n"
        )
        handle.write("|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|\n")
        for row in rows:
            def g(x, d=4):
                return f"{float(x):.{d}f}" if x is not None else "n/a"
            handle.write(
                f"| {row['run']} | {row['dataset']} | {row['corruption']} | {row['seed']} | "
                f"{g(row['ari'])} | {g(row['nmi'])} | {g(row['effective_dim_pr'],1)} | "
                f"{g(row['cluster_aligned_eff_dim'],2)} | {g(row['mean_effective_change_rate'],3)} | "
                f"{g(row['donor_pool_mean'],0)} | {g(row['donor_pool_min'],0)} |\n"
            )

        handle.write("\n## Decision rule\n\n")
        handle.write("- **S1/S2/S3 > S0 (swap_global)**：ARI +>=0.02 且多种子稳定 -> "
                     "限制供体池确实逼模型学更细结构，是从 scMAE 内部改对了。\n")
        handle.write("- **S1/S2/S3 ~= S0**：nuisance 匹配无效 -> swap 捷径不是瓶颈。\n")
        handle.write("- **swap_global(S0) vs zero**：先确认 swap-noise 本身相对零填充的效应方向。\n")
        handle.write("- **donor_pool_mean / pool_min 过小（个位数）**：swap 退化成近似恒等，"
                     "该 arm 结果不可信 —— 需调大 n_nuisance_bins 对应的箱宽（减少箱数）。\n")
        handle.write("- **eff_change 远低于 mask_prob**：供体与目标高度相同（池太窄或过稀疏），同样是退化信号。\n")


if __name__ == "__main__":
    raise SystemExit(main())
