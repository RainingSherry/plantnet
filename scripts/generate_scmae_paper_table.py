#!/usr/bin/env python3
"""Generate manuscript-ready overall tables from the unified scMAE benchmark."""

from __future__ import annotations

import argparse
import csv
import statistics
from collections import defaultdict
from pathlib import Path
from typing import Any


METRICS = ("acc", "nmi", "ari", "f1_macro")
SOURCE_LABELS = {
    "active_methods": "Methods benchmark",
    "tuned_nm_rg_contrast": "Ours",
}
OUR_FINAL_MODELS = {
    "neighbormix_scmae",
    "rg_neighbormix_scmae",
    "rg_neighbormix_scmae_contrast_safe",
}


def number(value: str | None) -> float | None:
    try:
        result = float(value or "")
    except ValueError:
        return None
    return result if result == result else None


def true(value: str | None) -> bool:
    return (value or "").strip().lower() in {"1", "true", "yes"}


def display_name(source: str, model: str) -> str:
    names = {
        "dec": "DEC",
        "attentionae_sc": "AttentionAE-sc",
        "desc": "DESC",
        "genecompass": "GeneCompass",
        "geneformer": "Geneformer",
        "leiden": "Leiden",
        "louvain": "Louvain",
        "phytocluster": "PhytoCluster",
        "sc3": "SC3",
        "sccdcg": "scCDCG",
        "scdcc": "scDCC",
        "scdeepcluster": "scDeepCluster",
        "scdsc": "scDSC",
        "scdsc_gse70256": "scDSC-GSE70256",
        "scgnn": "scGNN",
        "scgpt": "scGPT",
        "scmae": "scMAE",
        "scname": "scNAME",
        "scrcl": "scRCL",
        "scvi": "scVI",
        "sczidesk": "scziDesk",
        "neighbormix_scmae": "scVICAR-F",
        "rg_neighbormix_scmae": "scVICAR-T",
        "rg_neighbormix_scmae_contrast_safe": "scVICAR-T + CL",
        "canm_cut_reweighted_mix": "CutAware NeighborMix",
        "canm_cut_reweighted_mix_contrast": "CutAware NeighborMix + CL",
        "scmae_dec_stdfloor": "scMAE-DEC-StdFloor",
        "pca_kmeans_known_k": "PCA + KMeans (known K)",
    }
    name = names.get(model, model)
    return name


def metric_stats(rows: list[dict[str, str]]) -> dict[str, float]:
    result: dict[str, float] = {}
    for metric in METRICS:
        values = [value for row in rows if (value := number(row.get(metric))) is not None]
        result[f"{metric}_mean"] = statistics.fmean(values)
        result[f"{metric}_std"] = statistics.stdev(values) if len(values) > 1 else 0.0
    return result


def fmt(mean: float, std: float) -> str:
    return f"{mean:.4f} +/- {std:.4f}"


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    fields = [
        "rank_by_ari",
        "method",
        "source",
        "n_success",
        "acc_mean",
        "acc_std",
        "nmi_mean",
        "nmi_std",
        "ari_mean",
        "ari_std",
        "f1_macro_mean",
        "f1_macro_std",
        "note",
    ]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def write_markdown(path: Path, rows: list[dict[str, Any]], n_datasets: int, n_runs: int) -> None:
    lines = [
        "| Rank (ARI) | Method | Source | Runs | ACC | NMI | ARI | Macro-F1 |",
        "|---:|---|---|---:|---:|---:|---:|---:|",
    ]
    for row in rows:
        rank = row["rank_by_ari"] if row["rank_by_ari"] else "-"
        lines.append(
            "| {rank} | {method} | {source} | {runs} | {acc} | {nmi} | {ari} | {f1} |".format(
                rank=rank,
                method=row["method"],
                source=row["source"],
                runs=row["n_success"],
                acc=fmt(row["acc_mean"], row["acc_std"]),
                nmi=fmt(row["nmi_mean"], row["nmi_std"]),
                ari=fmt(row["ari_mean"], row["ari_std"]),
                f1=fmt(row["f1_macro_mean"], row["f1_macro_std"]),
            )
        )
    lines.extend(
        [
            "",
            f"Values are mean +/- standard deviation over {n_datasets} datasets x 3 seeds ({n_runs} runs).",
            f"Only methods with {n_runs}/{n_runs} successful non-fallback runs are ranked.",
        ]
    )
    path.write_text("\n".join(lines) + "\n")


def latex_escape(text: str) -> str:
    return text.replace("_", "\\_").replace("+", "+")


def write_latex(path: Path, rows: list[dict[str, Any]], n_datasets: int, n_runs: int) -> None:
    ranked = [row for row in rows if row["rank_by_ari"]]
    best = {metric: max(row[f"{metric}_mean"] for row in ranked) for metric in METRICS}
    lines = [
        "% Requires booktabs.",
        "\\begin{table}[p]",
        "\\centering",
        f"\\caption{{Development-stage clustering performance on {n_datasets} scMAE datasets across three seeds. Values are mean $\\pm$ standard deviation over {n_runs} runs. Methods with incomplete or fallback runs are excluded.}}",
        "\\label{tab:scmae-overall}",
        "\\resizebox{\\textwidth}{!}{%",
        "\\begin{tabular}{llrrrr}",
        "\\toprule",
        "Method & Source & ACC & NMI & ARI & Macro-F1 \\\\",
        "\\midrule",
    ]
    for row in rows:
        values = []
        for metric in METRICS:
            cell = f"{row[f'{metric}_mean']:.4f} $\\pm$ {row[f'{metric}_std']:.4f}"
            if row["rank_by_ari"] and row[f"{metric}_mean"] == best[metric]:
                cell = f"\\textbf{{{cell}}}"
            values.append(cell)
        lines.append(f"{latex_escape(row['method'])} & {latex_escape(row['source'])} & " + " & ".join(values) + " \\\\")
    lines.extend(["\\bottomrule", "\\end{tabular}%", "}", "\\end{table}"])
    path.write_text("\n".join(lines) + "\n")


def write_main_latex(path: Path, rows: list[dict[str, Any]], n_datasets: int, n_runs: int) -> None:
    """Write the full ranked benchmark in a two-column manuscript layout."""
    rows = [dict(row) for row in rows if row["method"] != "scVICAR-T + CL"]
    for rank, row in enumerate(rows, start=1):
        row["rank_by_ari"] = rank
    ranked = [row for row in rows if row["rank_by_ari"]]
    best = {metric: max(row[f"{metric}_mean"] for row in ranked) for metric in METRICS}
    second = {
        metric: sorted({row[f"{metric}_mean"] for row in ranked}, reverse=True)[1]
        for metric in METRICS
    }
    lines = [
        "% Requires booktabs.",
        "\\begin{table*}[t]",
        "\\centering",
        "\\scriptsize",
        "\\setlength{\\tabcolsep}{4pt}",
        (
            f"\\caption{{Primary development comparison on {n_datasets} datasets and {len(rows)} method configurations. "
            "Each row summarizes 45 runs (three seeds per dataset). "
            "Bold and underlined values denote the best and second-best mean, respectively. "
            "The exploratory contrastive extension is reported in the supplement.}"
        ),
        "\\label{tab:development-benchmark}",
        "\\begin{tabular}{rlrrrr}",
        "\\toprule",
        "Rank & Method & ACC & NMI & ARI & Macro-F1 \\\\",
        "\\midrule",
    ]
    for row in rows:
        values = []
        for metric in METRICS:
            cell = f"{row[f'{metric}_mean']:.4f} $\\pm$ {row[f'{metric}_std']:.4f}"
            if row[f"{metric}_mean"] == best[metric]:
                cell = f"\\textbf{{{cell}}}"
            elif row[f"{metric}_mean"] == second[metric]:
                cell = f"\\underline{{{cell}}}"
            values.append(cell)
        lines.append(
            f"{row['rank_by_ari']} & {latex_escape(row['method'])} & "
            + " & ".join(values)
            + " \\\\"
        )
    lines.extend(["\\bottomrule", "\\end{tabular}", "\\end{table*}"])
    path.write_text("\n".join(lines) + "\n")


def write_per_dataset(path: Path, rows: list[dict[str, Any]]) -> None:
    fields = ["source", "model", "method", "dataset", "n_seeds", *METRICS]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def write_effects(path: Path, per_dataset: list[dict[str, Any]]) -> None:
    by_method_dataset = {
        (row["model"], row["dataset"]): row for row in per_dataset
    }
    fields = [
        "model", "method", "reference", "metric", "n_datasets", "mean_delta",
        "median_delta", "wins", "ties", "losses", "min_delta", "max_delta",
    ]
    output = []
    for model in ("neighbormix_scmae", "rg_neighbormix_scmae"):
        datasets = sorted(
            dataset for candidate, dataset in by_method_dataset
            if candidate == model and ("scmae", dataset) in by_method_dataset
        )
        for metric in METRICS:
            deltas = [
                float(by_method_dataset[(model, dataset)][metric])
                - float(by_method_dataset[("scmae", dataset)][metric])
                for dataset in datasets
            ]
            output.append({
                "model": model,
                "method": display_name("tuned_nm_rg_contrast", model),
                "reference": "scMAE",
                "metric": metric,
                "n_datasets": len(datasets),
                "mean_delta": statistics.fmean(deltas),
                "median_delta": statistics.median(deltas),
                "wins": sum(delta > 1e-12 for delta in deltas),
                "ties": sum(abs(delta) <= 1e-12 for delta in deltas),
                "losses": sum(delta < -1e-12 for delta in deltas),
                "min_delta": min(deltas),
                "max_delta": max(deltas),
            })
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(output)


def write_principal_dataset_latex(path: Path, per_dataset: list[dict[str, Any]]) -> None:
    indexed = {
        (row["model"], row["dataset"]): row for row in per_dataset
    }
    datasets = sorted(
        dataset for model, dataset in indexed
        if model == "scmae"
        and ("neighbormix_scmae", dataset) in indexed
        and ("rg_neighbormix_scmae", dataset) in indexed
    )
    lines = [
        r"\begin{table}[p]",
        r"\centering",
        r"\small",
        r"\caption{Dataset-level ARI in the 15-dataset development benchmark. Values average three seeds.}",
        r"\label{tab:development-dataset-ari}",
        r"\begin{tabular}{lrrrrr}",
        r"\toprule",
        "Dataset & scMAE & scVICAR-F & $\\Delta$F & scVICAR-T & $\\Delta$T \\\\",
        r"\midrule",
    ]
    for dataset in datasets:
        base = float(indexed[("scmae", dataset)]["ari"])
        fixed = float(indexed[("neighbormix_scmae", dataset)]["ari"])
        topo = float(indexed[("rg_neighbormix_scmae", dataset)]["ari"])
        label = latex_escape(dataset).replace("Quake\\_", "Q. ")
        lines.append(
            f"{label} & {base:.4f} & {fixed:.4f} & {fixed-base:+.4f} & "
            f"{topo:.4f} & {topo-base:+.4f} \\\\"
        )
    lines.extend([r"\bottomrule", r"\end{tabular}", r"\end{table}"])
    path.write_text("\n".join(lines) + "\n")


def write_attempted_method_status(
    lightweight_root: Path, output_root: Path, ranked_models: set[str]
) -> None:
    master = lightweight_root / "all_runs_master.csv"
    if not master.is_file():
        raise FileNotFoundError(master)
    with master.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    counts: dict[str, dict[str, int]] = defaultdict(lambda: defaultdict(int))
    for row in rows:
        counts[row["model"]][row.get("status", "unknown")] += 1
    incomplete = []
    for model, status in sorted(counts.items()):
        if model in ranked_models:
            continue
        if status.get("success", 0) == 48 and sum(status.values()) == 48:
            continue
        incomplete.append({
            "model": display_name("active_methods", model),
            "success": status.get("success", 0),
            "failed": status.get("failed", 0),
            "fallback": status.get("fallback", 0),
            "pending": status.get("pending", 0),
            "records": sum(status.values()),
        })
    csv_path = output_root / "attempted_incomplete_methods.csv"
    with csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["model", "records", "success", "failed", "fallback", "pending"],
        )
        writer.writeheader()
        writer.writerows(incomplete)
    lines = [
        r"\begin{table}[p]",
        r"\centering",
        r"\small",
        r"\caption{Additional development attempts excluded from the ranked benchmark because coverage was incomplete or used fallback outputs.}",
        r"\label{tab:development-incomplete}",
        r"\begin{tabular}{lrrrrr}",
        r"\toprule",
        "Method & Records & Success & Failed & Fallback & Pending \\\\",
        r"\midrule",
    ]
    for row in incomplete:
        lines.append(
            f"{latex_escape(row['model'])} & {row['records']} & {row['success']} & "
            f"{row['failed']} & {row['fallback']} & {row['pending']} \\\\"
        )
    lines.extend([r"\bottomrule", r"\end{tabular}", r"\end{table}"])
    (output_root / "attempted_incomplete_methods.tex").write_text(
        "\n".join(lines) + "\n"
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--benchmark-root", type=Path, required=True)
    parser.add_argument("--output-root", type=Path)
    parser.add_argument("--exclude-datasets", nargs="*", default=[])
    parser.add_argument("--lightweight-root", type=Path)
    args = parser.parse_args()
    root = args.benchmark_root
    output_root = args.output_root or root
    output_root.mkdir(parents=True, exist_ok=True)
    with (root / "combined_run_master.csv").open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    excluded = set(args.exclude_datasets)
    rows = [row for row in rows if row.get("dataset") not in excluded]
    datasets = sorted({row["dataset"] for row in rows if row.get("dataset")})
    expected_runs = len(datasets) * 3

    groups: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        groups[(row["benchmark_source"], row["model"])].append(row)

    table: list[dict[str, Any]] = []
    per_dataset: list[dict[str, Any]] = []
    for (source, model), group in groups.items():
        if source == "retired_methods":
            continue
        if source == "tuned_nm_rg_contrast" and model not in OUR_FINAL_MODELS:
            continue
        success = [row for row in group if row.get("status") == "success"]
        if len(group) != expected_runs or len(success) != expected_runs or any(true(row.get("fallback")) or true(row.get("substitute_model_used")) for row in group):
            continue
        record: dict[str, Any] = {
            "rank_by_ari": "",
            "method": display_name(source, model),
            "source": SOURCE_LABELS[source],
            "n_success": expected_runs,
            "note": "",
        }
        record.update(metric_stats(success))
        table.append(record)
        by_dataset: dict[str, list[dict[str, str]]] = defaultdict(list)
        for row in success:
            by_dataset[row["dataset"]].append(row)
        for dataset, dataset_rows in sorted(by_dataset.items()):
            dataset_record: dict[str, Any] = {
                "source": source,
                "model": model,
                "method": display_name(source, model),
                "dataset": dataset,
                "n_seeds": len(dataset_rows),
            }
            for metric in METRICS:
                values = [number(row.get(metric)) for row in dataset_rows]
                dataset_record[metric] = statistics.fmean(
                    value for value in values if value is not None
                )
            per_dataset.append(dataset_record)

    table.sort(key=lambda row: row["ari_mean"], reverse=True)
    for rank, row in enumerate(table, start=1):
        row["rank_by_ari"] = rank

    write_csv(output_root / "paper_table_overall.csv", table)
    write_markdown(output_root / "paper_table_overall.md", table, len(datasets), expected_runs)
    write_latex(output_root / "paper_table_overall.tex", table, len(datasets), expected_runs)
    write_main_latex(output_root / "development_benchmark_main.tex", table, len(datasets), expected_runs)
    write_per_dataset(output_root / "per_dataset_method_metrics.csv", per_dataset)
    write_effects(output_root / "development_effects_vs_scmae.csv", per_dataset)
    write_principal_dataset_latex(output_root / "development_principal_dataset_ari.tex", per_dataset)
    if args.lightweight_root:
        write_attempted_method_status(
            args.lightweight_root,
            output_root,
            {row["model"] for row in per_dataset},
        )


if __name__ == "__main__":
    main()
