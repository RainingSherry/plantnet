#!/usr/bin/env python3
"""
Create paper-ready benchmark tables from all_results_master_table.csv.

Outputs:
  - paper_benchmark_summary_numeric.csv
  - paper_benchmark_table_formatted.csv
  - paper_benchmark_main.tex
  - paper_benchmark_compact_primary_metric.tex
  - paper_benchmark.xlsx (if openpyxl is available)

Example:
  python make_paper_benchmark_tables.py \
    --input experiment_reports/results_master_20260619/all_results_master_table.csv \
    --outdir experiment_reports/results_master_20260619/paper_tables/formal \
    --experiment-group canonical/formal_benchmark_20260607_08 \
    --metrics acc nmi ari f1_macro \
    --primary-metric acc
"""

from __future__ import annotations

import argparse
import math
import re
import sys
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd


DEFAULT_METRICS = ["acc", "nmi", "ari", "f1_macro"]

METRIC_LABELS = {
    "acc": "ACC",
    "nmi": "NMI",
    "ari": "ARI",
    "f1_macro": "Macro-F1",
    "macro_f1": "Macro-F1",
    "fmi": "FMI",
    "v_measure": "V-measure",
    "homogeneity": "Homogeneity",
    "completeness": "Completeness",
    "silhouette": "Silhouette",
}

# Edit these labels to match the naming convention in your paper.
METHOD_LABELS = {
    "nm_scmae_nomix": "scMAE w/o NeighBorMix",
    "snn_neighbormix": "scMAE + SNN-NeighBorMix",
    "random_beta_uniform_0.1": "Random beta U(0, 0.1)",
    "random_beta_uniform_0.05": "Random beta U(0, 0.05)",
    "fixed_beta_0.1": "Fixed beta = 0.1",
    "fixed_beta_0.2": "Fixed beta = 0.2",
    "global_random_neighbor_control": "Global-random neighbor control",
    "gaussian_noise_matched_anchor_target": "Gaussian-noise control",
    "dec": "DEC",
    "desc": "DESC",
    "kmeans": "K-means",
    "leiden": "Leiden",
    "louvain": "Louvain",
    "sc3": "SC3",
    "scdcc": "scDCC",
    "scdsc": "scDSC",
    "scname": "scNAME",
    "sczidesk": "scziDesk",
    "attentionae_sc": "AttentionAE-SC",
    "scgpt_plantnet": "scGPT",
    "apa_scmae": "APA-scMAE",
    "scmae": "scMAE",
    "scmae_raw": "scMAE",
    "scmae_nomix": "scMAE w/o Mix",
    "neighbormix_scmae": "scMAE + NeighborMix",
}

STATUS_OK = {"success", "completed_no_status_json", "finished", "completed"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Aggregate experiment-run CSV into paper-ready benchmark tables."
    )
    parser.add_argument("--input", required=True, help="CSV path, raw GitHub URL, or GitHub blob URL.")
    parser.add_argument("--outdir", required=True, help="Directory where output tables will be written.")
    parser.add_argument(
        "--experiment-group",
        default=None,
        help="Exact experiment_group to keep. Strongly recommended to avoid mixing unrelated experiments.",
    )
    parser.add_argument(
        "--experiment-group-contains",
        default=None,
        help="Regex/string filter for experiment_group when exact matching is inconvenient.",
    )
    parser.add_argument(
        "--dataset-group",
        default=None,
        help="Optional exact dataset_group to keep.",
    )
    parser.add_argument(
        "--stage",
        default=None,
        help="Optional exact stage to keep, e.g. full, stage1, stage2. Recommended for ablations.",
    )
    parser.add_argument(
        "--datasets",
        nargs="*",
        default=None,
        help="Optional whitelist of dataset names.",
    )
    parser.add_argument(
        "--include-method-regex",
        default=None,
        help="Only keep methods matching this regex.",
    )
    parser.add_argument(
        "--exclude-method-regex",
        default=None,
        help="Drop methods matching this regex.",
    )
    parser.add_argument(
        "--method-col",
        default="method",
        choices=["method", "variant", "model_key"],
        help="Column used as the method identifier.",
    )
    parser.add_argument(
        "--metrics",
        nargs="+",
        default=DEFAULT_METRICS,
        help="Metrics to include. Default: acc nmi ari f1_macro.",
    )
    parser.add_argument(
        "--primary-metric",
        default="acc",
        help="Metric used for compact table, sorting, and best/second-best highlighting.",
    )
    parser.add_argument(
        "--decimals",
        type=int,
        default=3,
        help="Decimal places for paper-formatted tables.",
    )
    parser.add_argument(
        "--valid-status",
        nargs="*",
        default=sorted(STATUS_OK),
        help="Statuses treated as usable runs.",
    )
    parser.add_argument(
        "--score-source",
        default=None,
        help="Optional exact score_source filter, e.g. eval_fixed.csv or metrics.json:top_level.",
    )
    parser.add_argument(
        "--keep-secondary",
        action="store_true",
        help="Do not filter is_primary == true.",
    )
    parser.add_argument(
        "--keep-non-run-rows",
        action="store_true",
        help="Do not filter row_granularity == run.",
    )
    parser.add_argument(
        "--no-deduplicate",
        action="store_true",
        help="Do not drop duplicate dataset/method/seed rows after filtering.",
    )
    parser.add_argument(
        "--split-stage",
        action="store_true",
        help="Always include stage in the grouping, even if --stage is provided.",
    )
    parser.add_argument(
        "--no-auto-split-stage",
        action="store_true",
        help="Do not automatically split by stage when multiple stages are present.",
    )
    parser.add_argument(
        "--sort-datasets",
        action="store_true",
        help="Sort dataset names alphabetically instead of preserving file order.",
    )
    parser.add_argument(
        "--sort-methods-by-primary",
        action="store_true",
        help="Within each dataset, sort methods by primary metric descending. Default preserves file order.",
    )
    return parser.parse_args()


def normalize_github_url(path_or_url: str) -> str:
    """Convert a GitHub blob URL into a raw URL so pandas can read it."""
    if path_or_url.startswith("https://github.com/") and "/blob/" in path_or_url:
        # https://github.com/OWNER/REPO/blob/BRANCH/path/to/file.csv
        tail = path_or_url.removeprefix("https://github.com/")
        owner, repo, _, rest = tail.split("/", 3)
        branch, file_path = rest.split("/", 1)
        return f"https://raw.githubusercontent.com/{owner}/{repo}/{branch}/{file_path}"
    return path_or_url


def as_bool(series: pd.Series) -> pd.Series:
    return series.astype(str).str.strip().str.lower().isin({"true", "1", "yes", "y"})


def natural_key(value: object) -> tuple[object, ...]:
    text = "" if pd.isna(value) else str(value)
    return tuple(int(chunk) if chunk.isdigit() else chunk.lower() for chunk in re.split(r"(\d+)", text))


def latex_escape(text: object) -> str:
    text = "" if pd.isna(text) else str(text)
    repl = {
        "\\": r"\textbackslash{}",
        "&": r"\&",
        "%": r"\%",
        "$": r"\$",
        "#": r"\#",
        "_": r"\_",
        "{": r"\{",
        "}": r"\}",
        "~": r"\textasciitilde{}",
        "^": r"\textasciicircum{}",
    }
    return "".join(repl.get(ch, ch) for ch in text)


def pretty_method(name: object) -> str:
    if pd.isna(name):
        return "Unknown"
    raw = str(name)
    if raw in METHOD_LABELS:
        return METHOD_LABELS[raw]
    # Useful fallback for sweep names such as random_beta_uniform_0.1.
    pretty = raw.replace("_", " ").strip()
    pretty = re.sub(r"\bbeta\b", "beta", pretty, flags=re.IGNORECASE)
    return pretty[:1].upper() + pretty[1:]


def ordered_unique(series: pd.Series, sort_alpha: bool = False) -> list[str]:
    values = [str(v) for v in series.dropna().tolist()]
    seen: set[str] = set()
    out: list[str] = []
    for v in values:
        if v not in seen:
            seen.add(v)
            out.append(v)
    if sort_alpha:
        out.sort(key=natural_key)
    return out


def load_results(input_path: str) -> pd.DataFrame:
    src = normalize_github_url(input_path)
    try:
        df = pd.read_csv(src)
    except Exception as exc:
        raise SystemExit(f"Could not read CSV from {src!r}: {exc}") from exc
    df.columns = [str(c).strip() for c in df.columns]
    return df


def ensure_columns(df: pd.DataFrame, required: Iterable[str]) -> None:
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise SystemExit(f"Missing required columns: {missing}. Available columns: {list(df.columns)}")


def filter_results(df: pd.DataFrame, args: argparse.Namespace) -> pd.DataFrame:
    ensure_columns(df, ["dataset", args.method_col])
    out = df.copy()

    if "is_primary" in out.columns and not args.keep_secondary:
        out = out[as_bool(out["is_primary"])]

    if "row_granularity" in out.columns and not args.keep_non_run_rows:
        out = out[out["row_granularity"].astype(str).str.strip().eq("run")]

    if "status" in out.columns and args.valid_status:
        ok = {s.strip() for s in args.valid_status}
        out = out[out["status"].astype(str).str.strip().isin(ok)]

    if args.experiment_group is not None:
        ensure_columns(out, ["experiment_group"])
        out = out[out["experiment_group"].astype(str).eq(args.experiment_group)]

    if args.experiment_group_contains is not None:
        ensure_columns(out, ["experiment_group"])
        out = out[out["experiment_group"].astype(str).str.contains(args.experiment_group_contains, regex=True, na=False)]

    if args.dataset_group is not None:
        ensure_columns(out, ["dataset_group"])
        out = out[out["dataset_group"].astype(str).eq(args.dataset_group)]

    if args.stage is not None:
        ensure_columns(out, ["stage"])
        out = out[out["stage"].astype(str).eq(args.stage)]

    if args.score_source is not None:
        ensure_columns(out, ["score_source"])
        out = out[out["score_source"].astype(str).eq(args.score_source)]

    if args.datasets:
        allowed = set(args.datasets)
        out = out[out["dataset"].astype(str).isin(allowed)]

    method_series = out[args.method_col].astype(str)
    if args.include_method_regex:
        out = out[method_series.str.contains(args.include_method_regex, regex=True, na=False)]
        method_series = out[args.method_col].astype(str)
    if args.exclude_method_regex:
        out = out[~method_series.str.contains(args.exclude_method_regex, regex=True, na=False)]

    for metric in args.metrics:
        if metric in out.columns:
            out[metric] = pd.to_numeric(out[metric], errors="coerce")

    metric_cols = [m for m in args.metrics if m in out.columns]
    if not metric_cols:
        raise SystemExit(f"None of the requested metrics are present: {args.metrics}")

    if args.primary_metric not in metric_cols:
        raise SystemExit(
            f"Primary metric {args.primary_metric!r} is not available after filtering. "
            f"Available requested metrics: {metric_cols}"
        )

    out = out[out[metric_cols].notna().any(axis=1)]
    out = out[out[args.primary_metric].notna()]

    if out.empty:
        raise SystemExit("No usable rows remain after filtering. Relax filters or check experiment_group/stage/status.")

    if "seed" in out.columns:
        out["seed"] = out["seed"].astype(str).str.strip()

    if not args.no_deduplicate and "seed" in out.columns:
        # Keep the last natural record_id/run if duplicate seeds exist for the same reported condition.
        dedup_cols = ["dataset", args.method_col, "seed"]
        for optional in ["experiment_group", "dataset_group", "stage", "score_source"]:
            if optional in out.columns:
                dedup_cols.append(optional)
        sort_cols = [c for c in ["record_id", "source_path", "run_dir"] if c in out.columns]
        if sort_cols:
            out = out.sort_values(sort_cols, key=lambda col: col.map(natural_key))
        out = out.drop_duplicates(subset=dedup_cols, keep="last")

    return out


def should_group_by_stage(df: pd.DataFrame, args: argparse.Namespace) -> bool:
    if "stage" not in df.columns:
        return False
    if args.split_stage:
        return True
    if args.no_auto_split_stage:
        return False
    if args.stage is not None:
        return False
    stages = [s for s in df["stage"].fillna("").astype(str).str.strip().unique().tolist() if s and s.lower() != "nan"]
    return len(stages) > 1


def aggregate(df: pd.DataFrame, args: argparse.Namespace) -> tuple[pd.DataFrame, list[str]]:
    metric_cols = [m for m in args.metrics if m in df.columns]
    group_cols = ["dataset"]
    if should_group_by_stage(df, args):
        group_cols.append("stage")
    group_cols.append(args.method_col)

    agg_spec: dict[str, list[str] | str] = {}
    for metric in metric_cols:
        agg_spec[metric] = ["mean", "std", "count"]
    if "runtime_seconds" in df.columns:
        df["runtime_seconds"] = pd.to_numeric(df["runtime_seconds"], errors="coerce")
        if df["runtime_seconds"].notna().any():
            agg_spec["runtime_seconds"] = ["mean", "std"]
    if "seed" in df.columns:
        agg_spec["seed"] = lambda x: ",".join(sorted(set(map(str, x)), key=natural_key))

    summary = df.groupby(group_cols, dropna=False, sort=False).agg(agg_spec).reset_index()
    summary.columns = [
        "_".join([str(x) for x in col if str(x) != ""]).rstrip("_")
        if isinstance(col, tuple)
        else str(col)
        for col in summary.columns
    ]

    # One clear n_seeds column based on the primary metric count.
    summary["n_seeds"] = summary[f"{args.primary_metric}_count"].astype(int)
    summary["method_raw"] = summary[args.method_col].astype(str)
    summary["method"] = summary["method_raw"].map(pretty_method)
    if "stage" in group_cols:
        summary["method"] = summary["stage"].astype(str).map(lambda s: f"[{s}] ") + summary["method"]

    return summary, metric_cols


def rank_flags(summary: pd.DataFrame, metric_cols: list[str]) -> dict[tuple[str, str, str], str]:
    """Return {(dataset, method_raw, metric): 'best'|'second'} for metric means."""
    flags: dict[tuple[str, str, str], str] = {}
    for dataset, sub in summary.groupby("dataset", sort=False):
        for metric in metric_cols:
            mean_col = f"{metric}_mean"
            vals = sub[["method_raw", mean_col]].dropna()
            if vals.empty:
                continue
            unique_vals = sorted(vals[mean_col].unique(), reverse=True)
            best = unique_vals[0]
            second = unique_vals[1] if len(unique_vals) > 1 else None
            for _, row in vals.iterrows():
                val = row[mean_col]
                key = (str(dataset), str(row["method_raw"]), metric)
                if math.isclose(val, best, rel_tol=1e-12, abs_tol=1e-12):
                    flags[key] = "best"
                elif second is not None and math.isclose(val, second, rel_tol=1e-12, abs_tol=1e-12):
                    flags[key] = "second"
    return flags


def fmt_value(mean: float, std: float, n: int, decimals: int, latex: bool, flag: str | None = None) -> str:
    if pd.isna(mean):
        return "--"
    base = f"{mean:.{decimals}f}"
    if n >= 2 and not pd.isna(std):
        pm = r"$\pm$" if latex else " ± "
        base = f"{base}{pm}{std:.{decimals}f}"
    if latex:
        if flag == "best":
            base = rf"\textbf{{{base}}}"
        elif flag == "second":
            base = rf"\underline{{{base}}}"
    else:
        if flag == "best":
            base = f"*{base}*"
        elif flag == "second":
            base = f"_{base}_"
    return base


def sort_summary(summary: pd.DataFrame, args: argparse.Namespace) -> pd.DataFrame:
    out = summary.copy()
    if args.sort_datasets:
        dataset_order = ordered_unique(out["dataset"], sort_alpha=True)
    else:
        dataset_order = ordered_unique(out["dataset"], sort_alpha=False)
    out["_dataset_order"] = out["dataset"].astype(str).map({d: i for i, d in enumerate(dataset_order)})

    if args.sort_methods_by_primary:
        out = out.sort_values(["_dataset_order", f"{args.primary_metric}_mean"], ascending=[True, False])
    else:
        # Preserve first-seen method order within the filtered file.
        method_order = ordered_unique(out["method_raw"], sort_alpha=False)
        out["_method_order"] = out["method_raw"].map({m: i for i, m in enumerate(method_order)})
        out = out.sort_values(["_dataset_order", "_method_order"])
        out = out.drop(columns=["_method_order"])
    return out.drop(columns=["_dataset_order"])


def build_formatted_table(summary: pd.DataFrame, metric_cols: list[str], args: argparse.Namespace, latex: bool) -> pd.DataFrame:
    flags = rank_flags(summary, metric_cols)
    rows = []
    for _, row in summary.iterrows():
        record = {
            "Dataset": str(row["dataset"]),
            "Method": latex_escape(row["method"]) if latex else row["method"],
            "#Seeds": int(row["n_seeds"]),
        }
        for metric in metric_cols:
            record[METRIC_LABELS.get(metric, metric)] = fmt_value(
                row[f"{metric}_mean"],
                row.get(f"{metric}_std", np.nan),
                int(row["n_seeds"]),
                args.decimals,
                latex=latex,
                flag=flags.get((str(row["dataset"]), str(row["method_raw"]), metric)),
            )
        rows.append(record)
    cols = ["Dataset", "Method", "#Seeds"] + [METRIC_LABELS.get(m, m) for m in metric_cols]
    return pd.DataFrame(rows, columns=cols)


def dataframe_to_latex_booktabs(table: pd.DataFrame, caption: str, label: str) -> str:
    col_spec = "llr" + "c" * max(0, table.shape[1] - 3)
    lines = [
        r"\begin{table*}[t]",
        r"\centering",
        r"\small",
        rf"\caption{{{latex_escape(caption)}}}",
        rf"\label{{{latex_escape(label)}}}",
        rf"\begin{{tabular}}{{{col_spec}}}",
        r"\toprule",
        " & ".join(latex_escape(c) for c in table.columns) + r" \\",
        r"\midrule",
    ]
    prev_dataset = None
    for _, row in table.iterrows():
        dataset = row["Dataset"]
        if prev_dataset is not None and dataset != prev_dataset:
            lines.append(r"\midrule")
        values = []
        for col in table.columns:
            val = row[col]
            if col == "Dataset":
                val = latex_escape(val)
            elif col == "Method":
                # Already escaped by build_formatted_table(latex=True).
                val = str(val)
            elif col == "#Seeds":
                val = str(val)
            else:
                val = str(val)
            values.append(val)
        lines.append(" & ".join(values) + r" \\")
        prev_dataset = dataset
    lines += [
        r"\bottomrule",
        r"\end{tabular}",
        r"\vspace{2pt}",
        r"\footnotesize Values are mean$\pm$std over random seeds. Best results per dataset and metric are bold; second-best results are underlined.",
        r"\end{table*}",
        "",
    ]
    return "\n".join(lines)


def build_compact_primary(summary: pd.DataFrame, args: argparse.Namespace, latex: bool) -> pd.DataFrame:
    metric = args.primary_metric
    flags = rank_flags(summary, [metric])
    datasets = ordered_unique(summary["dataset"], sort_alpha=args.sort_datasets)
    methods = ordered_unique(summary["method_raw"], sort_alpha=False)
    method_label = dict(zip(summary["method_raw"], summary["method"]))

    rows = []
    for method_raw in methods:
        sub = summary[summary["method_raw"].eq(method_raw)]
        if sub.empty:
            continue
        row_out = {"Method": latex_escape(method_label.get(method_raw, method_raw)) if latex else method_label.get(method_raw, method_raw)}
        values_for_avg = []
        for dataset in datasets:
            hit = sub[sub["dataset"].astype(str).eq(dataset)]
            if hit.empty:
                row_out[dataset] = "--"
                continue
            r = hit.iloc[0]
            row_out[dataset] = fmt_value(
                r[f"{metric}_mean"],
                r.get(f"{metric}_std", np.nan),
                int(r["n_seeds"]),
                args.decimals,
                latex=latex,
                flag=flags.get((str(dataset), str(method_raw), metric)),
            )
            values_for_avg.append(r[f"{metric}_mean"])
        row_out["Avg."] = f"{np.nanmean(values_for_avg):.{args.decimals}f}" if values_for_avg else "--"
        rows.append(row_out)
    return pd.DataFrame(rows, columns=["Method"] + datasets + ["Avg."])


def compact_to_latex_booktabs(table: pd.DataFrame, caption: str, label: str) -> str:
    col_spec = "l" + "c" * (table.shape[1] - 1)
    lines = [
        r"\begin{table*}[t]",
        r"\centering",
        r"\scriptsize",
        rf"\caption{{{latex_escape(caption)}}}",
        rf"\label{{{latex_escape(label)}}}",
        rf"\begin{{tabular}}{{{col_spec}}}",
        r"\toprule",
        " & ".join(latex_escape(c) for c in table.columns) + r" \\",
        r"\midrule",
    ]
    for _, row in table.iterrows():
        values = []
        for col in table.columns:
            if col == "Method":
                values.append(str(row[col]))
            else:
                values.append(str(row[col]))
        lines.append(" & ".join(values) + r" \\")
    lines += [
        r"\bottomrule",
        r"\end{tabular}",
        r"\vspace{2pt}",
        r"\footnotesize Compact view for the primary metric. Best per dataset is bold; second-best is underlined.",
        r"\end{table*}",
        "",
    ]
    return "\n".join(lines)


def write_excel(out_xlsx: Path, filtered: pd.DataFrame, summary: pd.DataFrame, formatted: pd.DataFrame, compact: pd.DataFrame) -> None:
    try:
        import openpyxl  # noqa: F401
    except Exception:
        print("[warn] openpyxl is not installed; skipping .xlsx output.", file=sys.stderr)
        return

    with pd.ExcelWriter(out_xlsx, engine="openpyxl") as writer:
        filtered.to_excel(writer, sheet_name="filtered_runs", index=False)
        summary.to_excel(writer, sheet_name="summary_numeric", index=False)
        formatted.to_excel(writer, sheet_name="paper_table", index=False)
        compact.to_excel(writer, sheet_name="compact_primary", index=False)

    # Light, non-intrusive formatting.
    wb = openpyxl.load_workbook(out_xlsx)
    for ws in wb.worksheets:
        ws.freeze_panes = "A2"
        ws.auto_filter.ref = ws.dimensions
        for cell in ws[1]:
            cell.font = openpyxl.styles.Font(bold=True)
        for col_cells in ws.columns:
            max_len = max(len(str(c.value)) if c.value is not None else 0 for c in col_cells)
            ws.column_dimensions[col_cells[0].column_letter].width = min(max(max_len + 2, 10), 45)
    wb.save(out_xlsx)


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    raw = load_results(args.input)
    filtered = filter_results(raw, args)
    summary, metric_cols = aggregate(filtered, args)
    summary = sort_summary(summary, args)

    formatted_csv = build_formatted_table(summary, metric_cols, args, latex=False)
    formatted_latex = build_formatted_table(summary, metric_cols, args, latex=True)
    compact_latex = build_compact_primary(summary, args, latex=True)
    compact_csv = build_compact_primary(summary, args, latex=False)

    numeric_path = outdir / "paper_benchmark_summary_numeric.csv"
    formatted_path = outdir / "paper_benchmark_table_formatted.csv"
    compact_path = outdir / "paper_benchmark_compact_primary_metric.csv"
    tex_path = outdir / "paper_benchmark_main.tex"
    compact_tex_path = outdir / "paper_benchmark_compact_primary_metric.tex"
    xlsx_path = outdir / "paper_benchmark.xlsx"

    summary.to_csv(numeric_path, index=False)
    formatted_csv.to_csv(formatted_path, index=False)
    compact_csv.to_csv(compact_path, index=False)

    caption = "Benchmark results across datasets."
    if args.experiment_group:
        caption += f" Experiment group: {args.experiment_group}."
    if args.stage:
        caption += f" Stage: {args.stage}."
    tex_path.write_text(
        dataframe_to_latex_booktabs(formatted_latex, caption, "tab:benchmark_main"),
        encoding="utf-8",
    )
    compact_caption = f"Compact benchmark comparison using {METRIC_LABELS.get(args.primary_metric, args.primary_metric)}."
    compact_tex_path.write_text(
        compact_to_latex_booktabs(compact_latex, compact_caption, "tab:benchmark_compact_primary"),
        encoding="utf-8",
    )

    write_excel(xlsx_path, filtered, summary, formatted_csv, compact_csv)

    print("Done. Wrote:")
    for p in [numeric_path, formatted_path, compact_path, tex_path, compact_tex_path, xlsx_path]:
        if p.exists():
            print(f"  {p}")

    # Useful run summary for sanity checking.
    print("\nFiltered runs:", len(filtered))
    print("Datasets:", ", ".join(ordered_unique(filtered["dataset"], sort_alpha=args.sort_datasets)))
    print("Methods:", ", ".join(ordered_unique(filtered[args.method_col], sort_alpha=False)))
    if "stage" in filtered.columns:
        stages = [s for s in ordered_unique(filtered["stage"], sort_alpha=False) if s and s.lower() != "nan"]
        if stages:
            print("Stages:", ", ".join(stages))


if __name__ == "__main__":
    main()
