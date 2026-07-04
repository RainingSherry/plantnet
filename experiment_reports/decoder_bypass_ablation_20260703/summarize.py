#!/usr/bin/env python3
"""汇总 decoder-bypass 消融：按 (dataset, decoder_mode) 聚合，并算相对 concat(D0) 的 delta。"""
from __future__ import annotations

import csv
import json
import statistics
from collections import defaultdict
from pathlib import Path


def mean_sd(vals):
    vals = [v for v in vals if v is not None]
    if not vals:
        return float("nan"), float("nan")
    if len(vals) == 1:
        return float(vals[0]), 0.0
    return float(statistics.mean(vals)), float(statistics.stdev(vals))


def fmt(vals, d=4):
    m, s = mean_sd(vals)
    return f"{m:.{d}f} +/- {s:.{d}f}"


def main() -> int:
    root = Path(__file__).resolve().parent
    rows = []
    for sp in sorted((root / "runs").glob("*/summary.json")):
        s = json.load(open(sp))
        m = s.get("fixed_metrics", {}).get("kmeans_known_k", {})
        std = s.get("std_profile", {})
        al = s.get("cluster_aligned", {})
        rows.append({
            "run": sp.parent.name,
            "dataset": s.get("dataset"),
            "decoder_mode": s.get("decoder_mode"),
            "seed": s.get("seed"),
            "ari": m.get("ari"), "nmi": m.get("nmi"), "acc": m.get("acc"),
            "effective_dim_pr": std.get("effective_dim_pr"),
            "cluster_aligned_eff_dim": al.get("cluster_aligned_eff_dim"),
            "std_min": std.get("std_min"),
            "final_base_loss": s.get("final_base_loss"),
            "runtime_seconds": s.get("runtime_seconds"),
        })

    csv_path = root / "summary.csv"
    if rows:
        with open(csv_path, "w", newline="") as h:
            w = csv.DictWriter(h, fieldnames=list(rows[0].keys()))
            w.writeheader()
            w.writerows(rows)

    md = root / "SUMMARY.md"
    with open(md, "w") as h:
        h.write("# Decoder-bypass 消融汇总\n\n")
        if not rows:
            h.write("runs/*/summary.json 尚无完成的 run。\n")
            print(f"Wrote {md} (empty)")
            return 0
        by = defaultdict(list)
        for r in rows:
            by[(r["dataset"], r["decoder_mode"])].append(r)

        h.write("## 按 (数据集 × decoder_mode) 聚合\n\n")
        h.write("| 数据集 | decoder_mode | n | ARI mean±sd | NMI mean±sd | eff_dim | aligned_eff_dim |\n")
        h.write("|---|---|---:|---:|---:|---:|---:|\n")
        order = {"concat": 0, "none": 1, "lowrank": 2}
        for ds in sorted({r["dataset"] for r in rows}):
            for mode in sorted({m for (d, m) in by if d == ds}, key=lambda x: order.get(x, 9)):
                g = by[(ds, mode)]
                tag = " (D0=赢家对照)" if mode == "concat" else ""
                h.write(
                    f"| {ds} | {mode}{tag} | {len(g)} | {fmt([r['ari'] for r in g])} | "
                    f"{fmt([r['nmi'] for r in g])} | {fmt([r['effective_dim_pr'] for r in g],1)} | "
                    f"{fmt([r['cluster_aligned_eff_dim'] for r in g],2)} |\n"
                )

        h.write("\n## 相对 concat (D0 原始 scMAE decoder) 的 ARI delta\n\n")
        h.write("| 数据集 | decoder_mode | ARI mean | delta vs D0 | 判定(|Δ|≥0.02) |\n|---|---|---:|---:|---|\n")
        for ds in sorted({r["dataset"] for r in rows}):
            base = by.get((ds, "concat"))
            if not base:
                continue
            b, _ = mean_sd([r["ari"] for r in base])
            for mode in sorted({m for (d, m) in by if d == ds}, key=lambda x: order.get(x, 9)):
                mn, _ = mean_sd([r["ari"] for r in by[(ds, mode)]])
                d = mn - b
                verdict = "↑改善" if d >= 0.02 else ("↓变差" if d <= -0.02 else "≈持平")
                h.write(f"| {ds} | {mode} | {mn:.4f} | {d:+.4f} | {verdict} |\n")

        h.write("\n## 每-run 明细\n\n")
        h.write("| run | ARI | NMI | eff_dim | aligned | std_min | base_loss |\n|---|---:|---:|---:|---:|---:|---:|\n")
        for r in rows:
            def g(x, dd=4):
                return f"{float(x):.{dd}f}" if x is not None else "n/a"
            h.write(f"| {r['run']} | {g(r['ari'])} | {g(r['nmi'])} | {g(r['effective_dim_pr'],1)} | "
                    f"{g(r['cluster_aligned_eff_dim'],2)} | {g(r['std_min'],3)} | {g(r['final_base_loss'],4)} |\n")

        h.write("\n## 判据\n\n")
        h.write("- **D1(none)/D2(lowrank) 相对 D0(concat) ARI +≥0.02 且多种子多数据集稳定** → 原始 decoder 的 G 维 mask 旁路是结构捷径，去掉/压缩它能得到更适合聚类的 embedding。\n")
        h.write("- **D1/D2 ≈ D0** → mask 旁路无害无益，原 scMAE decoder 结构对聚类不是瓶颈。\n")
        h.write("- **D1/D2 < D0** → mask 旁路对重构/聚类是有益的，不应移除。\n")
        h.write("- concat 应复现 DEC+std-floor 赢家（Macosko ARI≈0.70）作为 sanity check。\n")
    print(f"Wrote {csv_path}")
    print(f"Wrote {md}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
