#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
MANUSCRIPT_ROOT = SCRIPT_DIR.parent
GENERATED_ROOT = MANUSCRIPT_ROOT / "generated"
DEFAULT_RUN_ROOT = Path("/tmp/caam_feature_space_smoke/dev_20260626_gpu")

RUNS = {
    "hvg2000": "Quake_Smart-seq2_Lung__hvg2000__control__seed42__epochs3",
    "full_gene_stress": "Quake_Smart-seq2_Lung__full_gene_stress__control__seed42__epochs3",
}


def load_json(path: Path) -> dict:
    with open(path, "r", encoding="utf-8") as handle:
        return json.load(handle)


def write_json(path: Path, data: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, sort_keys=True), encoding="utf-8")


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        raise ValueError(f"No rows to write for {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(rows[0].keys())
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def fmt(value: float) -> str:
    return f"{float(value):.6f}"


def latex_escape(value: object) -> str:
    return str(value).replace("_", r"\_")


def flatten_run(run_root: Path, role: str, run_name: str) -> dict:
    run_dir = run_root / run_name
    required = ("metrics.json", "dataset_profile.json", "run_manifest.json", "preprocess_config.json", "artifact_manifest.json")
    missing = [name for name in required if not (run_dir / name).exists()]
    if missing:
        raise FileNotFoundError(f"Missing required feature-space smoke files for {role}: {missing} in {run_dir}")

    metrics = load_json(run_dir / "metrics.json")
    profile = load_json(run_dir / "dataset_profile.json")
    run_manifest = load_json(run_dir / "run_manifest.json")
    preprocess = load_json(run_dir / "preprocess_config.json")
    artifact = load_json(run_dir / "artifact_manifest.json")
    return {
        "role": role,
        "dataset": str(artifact["dataset"]),
        "seed": int(artifact["seed"]),
        "epochs": int(load_json(run_dir / "args.json")["epochs"]),
        "feature_space_source": str(profile["feature_space_source"]),
        "n_top_genes": int(preprocess["n_top_genes"]),
        "n_cells": int(profile["n_cells"]),
        "n_genes": int(profile["n_genes"]),
        "n_genes_original": int(profile["n_genes_original"]),
        "student_trainable_params": int(run_manifest["student_trainable_params"]),
        "generator_trainable_params": int(run_manifest["generator_trainable_params"]),
        "kmeans_known_k_ari": float(metrics["kmeans_known_k"]["ari"]),
        "kmeans_known_k_nmi": float(metrics["kmeans_known_k"]["nmi"]),
        "kmeans_known_k_acc": float(metrics["kmeans_known_k"]["acc"]),
        "kmeans_known_k_f1_macro": float(metrics["kmeans_known_k"]["f1_macro"]),
        "leiden_fixed_ari": float(metrics["leiden_fixed"]["ari"]),
        "leiden_fixed_nmi": float(metrics["leiden_fixed"]["nmi"]),
        "leiden_fixed_acc": float(metrics["leiden_fixed"]["acc"]),
        "leiden_fixed_f1_macro": float(metrics["leiden_fixed"]["f1_macro"]),
        "run_dir": str(run_dir),
    }


def paired_deltas(rows: list[dict]) -> dict:
    by_role = {row["role"]: row for row in rows}
    hvg = by_role["hvg2000"]
    full = by_role["full_gene_stress"]
    return {
        "full_minus_hvg_known_k_ari": full["kmeans_known_k_ari"] - hvg["kmeans_known_k_ari"],
        "full_minus_hvg_known_k_nmi": full["kmeans_known_k_nmi"] - hvg["kmeans_known_k_nmi"],
        "full_minus_hvg_known_k_f1_macro": full["kmeans_known_k_f1_macro"] - hvg["kmeans_known_k_f1_macro"],
        "full_minus_hvg_leiden_fixed_ari": full["leiden_fixed_ari"] - hvg["leiden_fixed_ari"],
        "full_minus_hvg_leiden_fixed_nmi": full["leiden_fixed_nmi"] - hvg["leiden_fixed_nmi"],
        "full_minus_hvg_leiden_fixed_f1_macro": full["leiden_fixed_f1_macro"] - hvg["leiden_fixed_f1_macro"],
        "full_over_hvg_parameter_ratio": full["student_trainable_params"] / hvg["student_trainable_params"],
        "full_over_hvg_gene_ratio": full["n_genes"] / hvg["n_genes"],
    }


def write_latex(path: Path, rows: list[dict], deltas: dict) -> None:
    lines = [
        r"\begin{table}[ht]",
        r"\centering",
        r"\caption{Feature-space development smoke on Quake\_Smart-seq2\_Lung. This is a single-dataset development diagnostic, not validation evidence.}",
        r"\label{tab:feature-space-smoke}",
        r"\scriptsize",
        r"\begin{tabular}{lrrrrr}",
        r"\toprule",
        r"Feature space & Genes & Student params & known-\(K\) ARI & fixed-Leiden ARI & fixed-Leiden F1 \\",
        r"\midrule",
    ]
    labels = {"hvg2000": "HVG 2000", "full_gene_stress": "Full-gene stress"}
    for row in rows:
        lines.append(
            f"{labels[row['role']]} & {row['n_genes']} & {row['student_trainable_params']} & "
            f"{fmt(row['kmeans_known_k_ari'])} & {fmt(row['leiden_fixed_ari'])} & {fmt(row['leiden_fixed_f1_macro'])} \\\\"
        )
    lines.extend(
        [
            r"\midrule",
            rf"Full minus HVG & -- & {deltas['full_over_hvg_parameter_ratio']:.2f}$\times$ & "
            rf"{fmt(deltas['full_minus_hvg_known_k_ari'])} & {fmt(deltas['full_minus_hvg_leiden_fixed_ari'])} & "
            rf"{fmt(deltas['full_minus_hvg_leiden_fixed_f1_macro'])} \\",
            r"\bottomrule",
            r"\end{tabular}",
            r"\end{table}",
            "",
        ]
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")


def write_markdown(path: Path, rows: list[dict], deltas: dict) -> None:
    labels = {"hvg2000": "HVG 2000", "full_gene_stress": "Full-gene stress"}
    lines = [
        "# Feature-space Smoke Summary",
        "",
        "Status: generated development-only feature-space smoke. This is not validation evidence.",
        "",
        "| feature space | genes | student params | known-K ARI | known-K NMI | fixed-Leiden ARI | fixed-Leiden F1 |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    for row in rows:
        lines.append(
            f"| {labels[row['role']]} | {row['n_genes']} | {row['student_trainable_params']} | "
            f"{fmt(row['kmeans_known_k_ari'])} | {fmt(row['kmeans_known_k_nmi'])} | "
            f"{fmt(row['leiden_fixed_ari'])} | {fmt(row['leiden_fixed_f1_macro'])} |"
        )
    lines.extend(
        [
            "",
            "## Paired Deltas",
            "",
            f"- Full-gene minus HVG known-K ARI: {fmt(deltas['full_minus_hvg_known_k_ari'])}.",
            f"- Full-gene minus HVG fixed-Leiden ARI: {fmt(deltas['full_minus_hvg_leiden_fixed_ari'])}.",
            f"- Full-gene minus HVG fixed-Leiden F1: {fmt(deltas['full_minus_hvg_leiden_fixed_f1_macro'])}.",
            f"- Full-gene/HVG student parameter ratio: {deltas['full_over_hvg_parameter_ratio']:.2f}x.",
            f"- Full-gene/HVG gene-count ratio: {deltas['full_over_hvg_gene_ratio']:.2f}x.",
            "",
            "## Claim Boundary",
            "",
            "This supports keeping HVG 2000 as the current dense-MLP protocol default under development evidence. It does not validate HVG 2000, reject all full-gene approaches, or evaluate sparse/gene-token full-gene architectures.",
            "",
        ]
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")


def build_outputs(args: argparse.Namespace) -> dict:
    rows = [flatten_run(args.run_root, role, run_name) for role, run_name in RUNS.items()]
    deltas = paired_deltas(rows)
    payload = {
        "status": "development_only",
        "run_root": str(args.run_root),
        "n_runs": len(rows),
        "rows": rows,
        "deltas": deltas,
        "claim_boundary": (
            "single-dataset development smoke; not validation; not a universal claim against full-gene or sparse/gene-token models"
        ),
    }
    args.output_dir.mkdir(parents=True, exist_ok=True)
    write_csv(args.output_dir / "feature_space_smoke.csv", rows)
    write_json(args.output_dir / "feature_space_smoke.json", payload)
    write_markdown(args.output_dir / "feature_space_smoke.md", rows, deltas)
    write_latex(args.output_dir / "feature_space_smoke.tex", rows, deltas)
    return payload


def main() -> int:
    parser = argparse.ArgumentParser(description="Summarize the CAAM/scMAE feature-space development smoke.")
    parser.add_argument("--run-root", type=Path, default=DEFAULT_RUN_ROOT)
    parser.add_argument("--output-dir", type=Path, default=GENERATED_ROOT / "feature_space_smoke")
    args = parser.parse_args()
    payload = build_outputs(args)
    print(f"feature_space_smoke_runs={payload['n_runs']}")
    print(f"feature_space_smoke_output={args.output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
