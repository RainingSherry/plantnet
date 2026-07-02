#!/usr/bin/env python3
"""Build the latest non-ablation paper benchmark tables.

This script intentionally excludes mechanism/ablation experiments and keeps
only benchmark-like results from `results/` plus the latest compatible result
CSVs under `experiment_reports/`.
"""

from __future__ import annotations

import importlib.util
import json
import math
import subprocess
from collections import Counter
from pathlib import Path
from typing import Any

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = PROJECT_ROOT / "experiment_reports" / "paper_benchmark_latest_20260629"
MASTER_BUILDER = PROJECT_ROOT / "experiment_reports" / "results_master_20260619" / "build_results_master.py"
PAPER_SCRIPT = PROJECT_ROOT / "make_paper_benchmark_tables.py"

ALLOWED_RESULT_GROUPS = {
    "canonical/formal_benchmark_20260607_08",
    "canonical/scmae_11datasets_20260609_12",
    "experiments/scgpt_plantnet_20260615",
}

EXCLUDED_GROUP_HINTS = (
    "ablation",
    "beta_mechanism",
    "stochastic_regularization",
    "cutaware",
    "neighbormix_ra_rg",
    "rc_nm_checkpoint",
    "smoke",
)

SCORE_COLUMNS = [
    "acc",
    "nmi",
    "ari",
    "f1_macro",
    "macro_f1",
    "fmi",
    "v_measure",
    "homogeneity",
    "completeness",
    "silhouette",
    "n_pred_clusters",
]


def load_master_builder():
    spec = importlib.util.spec_from_file_location("results_master_builder", MASTER_BUILDER)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load {MASTER_BUILDER}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def clean(value: Any) -> Any:
    if value is None:
        return ""
    try:
        if isinstance(value, float) and math.isnan(value):
            return ""
    except Exception:
        pass
    return value


def dataset_name(value: Any) -> str:
    text = str(clean(value)).strip()
    for prefix in ("processed__", "other__"):
        if text.startswith(prefix):
            return text.removeprefix(prefix)
    return text


def row_template(columns: list[str]) -> dict[str, Any]:
    return {column: "" for column in columns}


def as_output_row(raw: dict[str, Any], columns: list[str]) -> dict[str, Any]:
    row = row_template(columns)
    for key, value in raw.items():
        if key in row:
            row[key] = clean(value)
    return row


def add_desc_scname_rows(columns: list[str]) -> list[dict[str, Any]]:
    path = PROJECT_ROOT / "experiment_reports" / "desc_scname_fix_benchmark_20260628" / "desc_scname_updated_values_all_runs.csv"
    if not path.exists():
        return []
    df = pd.read_csv(path)
    rows: list[dict[str, Any]] = []
    for idx, raw in df.iterrows():
        raw = raw.to_dict()
        row = row_template(columns)
        method = str(raw.get("method", "")).strip()
        row.update(
            {
                "record_kind": "experiment_report_csv",
                "row_granularity": "run",
                "experiment_group": "experiment_reports/desc_scname_fix_benchmark_20260628",
                "dataset_group": "latest_desc_scname_fix",
                "dataset": dataset_name(raw.get("dataset")),
                "method": method,
                "variant": method,
                "model_key": method,
                "seed": clean(raw.get("seed")),
                "status": clean(raw.get("status")),
                "return_code": clean(raw.get("return_code")),
                "error": clean(raw.get("error")),
                "source_path": str(path.relative_to(PROJECT_ROOT)),
                "run_dir": clean(raw.get("save_dir")),
                "gpu": clean(raw.get("gpu")),
                "elapsed_seconds": clean(raw.get("elapsed_seconds")),
                "runtime_seconds": clean(raw.get("elapsed_seconds")),
                "command": clean(raw.get("command")),
                "score_source": path.name,
                "extra_json": json.dumps(
                    {
                        "value_tag": clean(raw.get("value_tag")),
                        "data_path": clean(raw.get("data_path")),
                        "epochs": clean(raw.get("epochs")),
                        "pretrain_epochs": clean(raw.get("pretrain_epochs")),
                        "n_clusters": clean(raw.get("n_clusters")),
                        "latest_replacement": True,
                    },
                    sort_keys=True,
                    default=str,
                ),
            }
        )
        for metric in SCORE_COLUMNS:
            row[metric] = clean(raw.get(metric))
        if row["macro_f1"] == "" and row["f1_macro"] != "":
            row["macro_f1"] = row["f1_macro"]
        rows.append(row)
    return rows


def add_apa_rows(columns: list[str]) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    path = PROJECT_ROOT / "experiment_reports" / "apa_scmae_full_benchmark_partial_20260628" / "tables" / "apa_scmae_full_metrics_by_run.csv"
    if not path.exists():
        return [], []
    df = pd.read_csv(path)
    evidence_rows: list[dict[str, Any]] = []
    paper_rows: list[dict[str, Any]] = []
    for _, raw in df.iterrows():
        raw = raw.to_dict()
        cluster_method = str(raw.get("cluster_method", "")).strip()
        method = "apa_scmae" if cluster_method == "kmeans_known_k" else f"apa_scmae_{cluster_method}"
        row = row_template(columns)
        row.update(
            {
                "record_kind": "experiment_report_csv",
                "row_granularity": "run",
                "experiment_group": "experiment_reports/apa_scmae_full_benchmark_partial_20260628",
                "dataset_group": "latest_apa_scmae_full_partial",
                "dataset": clean(raw.get("dataset")),
                "method": method,
                "variant": method,
                "model_key": method,
                "seed": clean(raw.get("seed")),
                "status": "success",
                "source_path": str(path.relative_to(PROJECT_ROOT)),
                "score_source": path.name,
                "extra_json": json.dumps(
                    {
                        "cluster_method": cluster_method,
                        "primary_model_readout": cluster_method == "kmeans_known_k",
                        "uses_known_k": clean(raw.get("uses_known_k")),
                        "oracle-K": clean(raw.get("oracle-K")),
                        "resolution": clean(raw.get("resolution")),
                        "n_neighbors": clean(raw.get("n_neighbors")),
                    },
                    sort_keys=True,
                    default=str,
                ),
            }
        )
        for metric in SCORE_COLUMNS:
            row[metric] = clean(raw.get(metric))
        if row["macro_f1"] == "" and row["f1_macro"] != "":
            row["macro_f1"] = row["f1_macro"]
        evidence_rows.append(row)
        if cluster_method == "kmeans_known_k":
            paper_rows.append(row.copy())
    return evidence_rows, paper_rows


def collect_result_rows(builder, columns: list[str]) -> list[dict[str, Any]]:
    rows = [builder.extract_run_row(path) for path in builder.iter_run_dirs()]
    selected = []
    for raw in rows:
        group = str(raw.get("experiment_group", ""))
        if group not in ALLOWED_RESULT_GROUPS:
            continue
        if any(hint in group.lower() for hint in EXCLUDED_GROUP_HINTS):
            continue
        row = as_output_row(raw, columns)
        selected.append(row)
    return selected


def assign_record_ids(rows: list[dict[str, Any]]) -> None:
    for idx, row in enumerate(rows, start=1):
        row["record_id"] = f"PB{idx:06d}"
        row["is_primary"] = "true"
        row["duplicate_reason"] = ""


def write_status_tables(rows: list[dict[str, Any]]) -> None:
    df = pd.DataFrame(rows)
    status = (
        df.groupby(["dataset", "method", "status"], dropna=False)
        .size()
        .reset_index(name="n_runs")
        .sort_values(["dataset", "method", "status"])
    )
    status.to_csv(OUT_DIR / "latest_benchmark_status_by_dataset_method.csv", index=False)

    coverage = (
        df.groupby(["dataset", "method"], dropna=False)
        .agg(
            n_records=("record_id", "count"),
            n_success=("status", lambda x: int((x.astype(str) == "success").sum())),
            statuses=("status", lambda x: ";".join(f"{k}:{v}" for k, v in sorted(Counter(map(str, x)).items()))),
            seeds=("seed", lambda x: ",".join(sorted({str(v) for v in x if str(v) and str(v) != "nan"}))),
        )
        .reset_index()
        .sort_values(["dataset", "method"])
    )
    coverage.to_csv(OUT_DIR / "latest_benchmark_coverage_matrix.csv", index=False)


def write_readme(evidence_rows: list[dict[str, Any]], paper_rows: list[dict[str, Any]]) -> None:
    evidence = pd.DataFrame(evidence_rows)
    paper = pd.DataFrame(paper_rows)
    valid_status = {"success", "completed_no_status_json", "finished", "completed"}
    if paper.empty:
        usable = paper
    else:
        usable = paper[
            paper.get("status", pd.Series(dtype=str)).astype(str).isin(valid_status)
            & paper[["acc", "nmi", "ari", "f1_macro"]].notna().any(axis=1)
            & paper["acc"].notna()
        ]
    lines = [
        "# Latest paper benchmark 2026-06-29",
        "",
        "This package consolidates the latest non-ablation benchmark results from `results/` and `experiment_reports/`.",
        "",
        "## Included sources",
        "",
        "- `results/canonical/formal_benchmark_20260607_08`",
        "- `results/canonical/scmae_11datasets_20260609_12`",
        "- `results/experiments/scgpt_plantnet_20260615`",
        "- `experiment_reports/desc_scname_fix_benchmark_20260628/desc_scname_updated_values_all_runs.csv`",
        "- `experiment_reports/apa_scmae_full_benchmark_partial_20260628/tables/apa_scmae_full_metrics_by_run.csv`",
        "",
        "## Exclusion rule",
        "",
        "Ablation/mechanism/smoke result groups are excluded, including beta mechanism, stochastic regularization, cut-aware, RA/RG sensitivity, RC checkpoint, and APA v2 ablation outputs.",
        "",
        "APA_scMAE `kmeans_known_k` is used as the primary model readout in the paper table. The `leiden_fixed` readout is kept in the evidence table only.",
        "",
        "## Counts",
        "",
        f"- evidence rows: {len(evidence_rows)}",
        f"- paper candidate rows: {len(paper_rows)}",
        f"- paper formatter usable scored rows: {len(usable)}",
        f"- evidence datasets: {evidence['dataset'].nunique() if not evidence.empty else 0}",
        f"- paper datasets: {paper['dataset'].nunique() if not paper.empty else 0}",
        f"- paper scored datasets: {usable['dataset'].nunique() if not usable.empty else 0}",
        f"- paper methods: {paper['method'].nunique() if not paper.empty else 0}",
        f"- paper scored methods: {usable['method'].nunique() if not usable.empty else 0}",
        "",
        "## Files",
        "",
        "- `latest_benchmark_evidence_all_rows.csv`: all selected non-ablation evidence rows, including skipped runs and APA secondary readouts.",
        "- `latest_benchmark_paper_runs.csv`: rows used to build the paper benchmark tables.",
        "- `paper_tables/`: numeric, formatted CSV, LaTeX, and XLSX benchmark tables.",
        "- `latest_benchmark_status_by_dataset_method.csv`: status counts.",
        "- `latest_benchmark_coverage_matrix.csv`: dataset/method coverage.",
        "",
    ]
    (OUT_DIR / "README.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    builder = load_master_builder()
    columns = list(builder.OUTPUT_COLUMNS)
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    result_rows = collect_result_rows(builder, columns)
    desc_rows = add_desc_scname_rows(columns)
    apa_evidence_rows, apa_paper_rows = add_apa_rows(columns)

    evidence_rows = result_rows + desc_rows + apa_evidence_rows
    paper_rows = result_rows + desc_rows + apa_paper_rows
    assign_record_ids(evidence_rows)
    assign_record_ids(paper_rows)

    evidence_path = OUT_DIR / "latest_benchmark_evidence_all_rows.csv"
    paper_path = OUT_DIR / "latest_benchmark_paper_runs.csv"
    pd.DataFrame(evidence_rows, columns=columns).to_csv(evidence_path, index=False)
    pd.DataFrame(paper_rows, columns=columns).to_csv(paper_path, index=False)
    write_status_tables(evidence_rows)
    write_readme(evidence_rows, paper_rows)

    paper_outdir = OUT_DIR / "paper_tables"
    cmd = [
        "python",
        str(PAPER_SCRIPT),
        "--input",
        str(paper_path),
        "--outdir",
        str(paper_outdir),
        "--metrics",
        "acc",
        "nmi",
        "ari",
        "f1_macro",
        "--primary-metric",
        "acc",
        "--sort-datasets",
    ]
    subprocess.run(cmd, check=True, cwd=PROJECT_ROOT)

    print(f"wrote: {OUT_DIR}")
    print(f"evidence rows: {len(evidence_rows)}")
    print(f"paper candidate rows: {len(paper_rows)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
