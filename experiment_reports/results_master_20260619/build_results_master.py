#!/usr/bin/env python3
"""Build comprehensive result tables from the project results tree.

The goal is coverage, not statistical interpretation.  The script scans
run-like artifact directories and compatible summary CSV files, then writes a
primary run table plus an evidence table that keeps every extractable score row
with provenance.
"""

from __future__ import annotations

import csv
import json
import math
import os
import re
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]
RESULTS_ROOT = (PROJECT_ROOT / "results").resolve()
OUT_DIR = Path(__file__).resolve().parent

RUN_MARKERS = {
    "metrics.json",
    "summary.json",
    "eval_fixed.csv",
    "args.json",
    "status.json",
    "command.txt",
    "run.log",
    "run_config.json",
}

RUN_DIR_MARKERS = {
    "metrics.json",
    "summary.json",
    "eval_fixed.csv",
    "args.json",
    "status.json",
    "run_config.json",
}

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

MEAN_SCORE_COLUMNS = [
    "acc_mean",
    "acc_std",
    "nmi_mean",
    "nmi_std",
    "ari_mean",
    "ari_std",
    "f1_macro_mean",
    "f1_macro_std",
    "macro_f1_mean",
    "macro_f1_std",
    "fmi_mean",
    "fmi_std",
    "v_measure_mean",
    "v_measure_std",
]

BASE_COLUMNS = [
    "record_id",
    "is_primary",
    "duplicate_reason",
    "record_kind",
    "row_granularity",
    "experiment_group",
    "dataset_group",
    "stage",
    "sweep",
    "dataset",
    "method",
    "variant",
    "model_key",
    "seed",
    "status",
    "return_code",
    "error",
    "n_success",
    "n_total",
    "status_summary",
    "authenticity",
    "substitute_model_used",
    "source_path",
    "run_dir",
    "metrics_path",
    "summary_path",
    "eval_path",
    "args_path",
    "status_path",
    "commit_sha",
    "branch",
    "gpu",
    "elapsed_seconds",
    "runtime_seconds",
    "command",
    "score_source",
]

OUTPUT_COLUMNS = BASE_COLUMNS + SCORE_COLUMNS + MEAN_SCORE_COLUMNS + [
    "extra_json",
]


def rel(path: Path) -> str:
    try:
        return path.resolve().relative_to(PROJECT_ROOT.resolve()).as_posix()
    except Exception:
        return path.as_posix()


def rel_results(path: Path) -> str:
    try:
        return path.resolve().relative_to(RESULTS_ROOT).as_posix()
    except Exception:
        return rel(path)


def safe_json(path: Path) -> dict[str, Any]:
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return {}
    return data if isinstance(data, dict) else {}


def clean_value(value: Any) -> Any:
    if value is None:
        return ""
    if isinstance(value, float) and math.isnan(value):
        return ""
    return value


def as_text(value: Any) -> str:
    value = clean_value(value)
    if value == "":
        return ""
    return str(value)


def first_present(*values: Any) -> Any:
    for value in values:
        value = clean_value(value)
        if value != "":
            return value
    return ""


def parse_method_seed_timestamp(name: str) -> tuple[str, str, str]:
    pieces = name.split("__")
    if len(pieces) < 3:
        return name, "", ""
    return pieces[0], pieces[1].removeprefix("seed"), pieces[2]


def experiment_group_from_rel(parts: tuple[str, ...]) -> str:
    if len(parts) >= 2 and parts[0] in {"canonical", "experiments"}:
        return f"{parts[0]}/{parts[1]}"
    if len(parts) >= 2 and parts[0] == "scratch":
        return "scratch/" + parts[1]
    if len(parts) >= 2 and parts[0] == "analysis":
        return "analysis/" + parts[1]
    return parts[0] if parts else ""


def parse_identity_from_path(run_dir: Path) -> dict[str, Any]:
    parts = tuple(Path(rel_results(run_dir)).parts)
    data: dict[str, Any] = {
        "experiment_group": experiment_group_from_rel(parts),
        "dataset_group": "",
        "stage": "",
        "sweep": "",
        "dataset": "",
        "method": "",
        "variant": "",
        "seed": "",
    }

    if len(parts) >= 4 and parts[0] == "canonical":
        method, seed, _timestamp = parse_method_seed_timestamp(parts[-1])
        data.update({"dataset_group": parts[1], "dataset": parts[2], "method": method, "variant": method, "seed": seed})
        return data

    if len(parts) >= 6 and parts[:2] == ("experiments", "neighbormix_beta_mechanism_20260617"):
        data.update(
            {
                "stage": parts[2],
                "dataset": parts[3],
                "method": parts[4],
                "variant": parts[4],
                "seed": parts[5].removeprefix("seed"),
            }
        )
        return data

    if len(parts) >= 5 and parts[:2] == ("experiments", "neighbormix_stochastic_regularization_20260616"):
        data.update({"dataset": parts[2], "method": parts[3], "variant": parts[3], "seed": parts[4].removeprefix("seed")})
        return data

    if len(parts) >= 5 and parts[:2] == ("experiments", "cutaware_neighbormix_20260615"):
        data.update({"dataset": parts[2], "method": parts[3], "variant": parts[3], "seed": parts[4].removeprefix("seed")})
        return data

    if len(parts) >= 6 and parts[:2] == ("experiments", "rc_nm_checkpoint_v4_1"):
        data.update(
            {
                "stage": parts[2],
                "dataset": parts[3],
                "seed": parts[4].removeprefix("seed"),
                "method": "RC_NeighborMix",
                "variant": parts[5],
            }
        )
        return data

    if len(parts) >= 7 and parts[:3] == ("experiments", "neighbormix_ra_rg", "rg_phase2_sensitivity_e80"):
        data.update(
            {
                "stage": parts[2],
                "sweep": parts[3],
                "dataset": parts[4],
                "method": parts[5],
                "variant": parts[5],
                "seed": parts[6].removeprefix("seed"),
            }
        )
        return data

    if len(parts) >= 6 and parts[:2] == ("experiments", "neighbormix_ra_rg"):
        data.update(
            {
                "stage": parts[2],
                "dataset": parts[3],
                "method": parts[4],
                "variant": parts[4],
                "seed": parts[5].removeprefix("seed"),
            }
        )
        return data

    if len(parts) >= 5 and parts[0] == "experiments":
        data.update({"dataset": parts[2], "method": parts[3], "variant": parts[3], "seed": parts[-1].removeprefix("seed")})
        return data

    if len(parts) >= 3 and parts[:2] == ("scratch", "smoke"):
        data.update({"stage": "smoke", "method": parts[2], "variant": parts[2]})
        return data

    return data


def flatten_metric_dict(data: dict[str, Any]) -> tuple[dict[str, Any], str]:
    if not data:
        return {}, ""

    top = {key: data.get(key, "") for key in SCORE_COLUMNS if key in data}
    if any(clean_value(value) != "" for value in top.values()):
        return top, "metrics.json:top_level"

    for preferred in ("kmeans_known_k", "fixed", "metrics", "fixed_metrics"):
        value = data.get(preferred)
        if isinstance(value, dict):
            nested = {key: value.get(key, "") for key in SCORE_COLUMNS if key in value}
            if any(clean_value(v) != "" for v in nested.values()):
                return nested, f"metrics.json:{preferred}"

    for key, value in data.items():
        if isinstance(value, dict):
            nested = {metric: value.get(metric, "") for metric in SCORE_COLUMNS if metric in value}
            if any(clean_value(v) != "" for v in nested.values()):
                return nested, f"metrics.json:{key}"

    return {}, ""


def flatten_summary_metrics(data: dict[str, Any]) -> tuple[dict[str, Any], str]:
    fixed = data.get("fixed_metrics")
    if isinstance(fixed, dict):
        for preferred in ("kmeans_known_k", "fixed"):
            value = fixed.get(preferred)
            if isinstance(value, dict):
                nested = {key: value.get(key, "") for key in SCORE_COLUMNS if key in value}
                if any(clean_value(v) != "" for v in nested.values()):
                    return nested, f"summary.json:fixed_metrics.{preferred}"
        for key, value in fixed.items():
            if isinstance(value, dict):
                nested = {metric: value.get(metric, "") for metric in SCORE_COLUMNS if metric in value}
                if any(clean_value(v) != "" for v in nested.values()):
                    return nested, f"summary.json:fixed_metrics.{key}"

    return flatten_metric_dict(data)


def read_eval_fixed(path: Path) -> tuple[dict[str, Any], str]:
    try:
        df = pd.read_csv(path)
    except Exception:
        return {}, ""
    if df.empty:
        return {}, ""
    row = df.iloc[0].to_dict()
    metrics = {key: row.get(key, "") for key in SCORE_COLUMNS if key in row}
    return metrics, "eval_fixed.csv"


def extract_run_metrics(run_dir: Path) -> tuple[dict[str, Any], str]:
    eval_path = run_dir / "eval_fixed.csv"
    if eval_path.exists():
        metrics, source = read_eval_fixed(eval_path)
        if metrics:
            return metrics, source

    metrics_path = run_dir / "metrics.json"
    if metrics_path.exists():
        metrics, source = flatten_metric_dict(safe_json(metrics_path))
        if metrics:
            return metrics, source

    summary_path = run_dir / "summary.json"
    if summary_path.exists():
        metrics, source = flatten_summary_metrics(safe_json(summary_path))
        if metrics:
            return metrics, source

    return {}, ""


def iter_run_dirs() -> list[Path]:
    run_dirs: list[Path] = []
    for dirpath, dirnames, filenames in os.walk(RESULTS_ROOT, followlinks=True):
        current = Path(dirpath)
        if current.name == "_index":
            dirnames[:] = []
            continue
        if set(filenames) & RUN_DIR_MARKERS:
            run_dirs.append(current)
    return sorted(set(run_dirs), key=lambda path: rel_results(path))


def model_key(method: Any, variant: Any) -> str:
    method_s = as_text(method)
    variant_s = as_text(variant)
    if variant_s and variant_s.lower() not in {"nan", "none", "null"}:
        if method_s in {"", "RC_NeighborMix", "NeighborMix_scMAE", "CutAware_NeighborMix_scMAE"}:
            return variant_s
        if variant_s != method_s and method_s.lower() in {
            "nm_scmae_nomix",
            "neighbormix_scmae",
            "scmae",
        }:
            return method_s
    return method_s or variant_s


def row_base() -> dict[str, Any]:
    return {key: "" for key in OUTPUT_COLUMNS}


def extract_run_row(run_dir: Path) -> dict[str, Any]:
    args = safe_json(run_dir / "args.json") if (run_dir / "args.json").exists() else {}
    status = safe_json(run_dir / "status.json") if (run_dir / "status.json").exists() else {}
    summary = safe_json(run_dir / "summary.json") if (run_dir / "summary.json").exists() else {}
    identity = parse_identity_from_path(run_dir)

    inferred_dataset = ""
    data_path = first_present(args.get("data_path"))
    if data_path:
        inferred_dataset = Path(str(data_path)).stem
    dataset = first_present(status.get("dataset"), args.get("dataset_name"), summary.get("dataset"), identity["dataset"], inferred_dataset)
    method = first_present(status.get("method"), identity["method"], args.get("method_name"), args.get("ablation_method"), summary.get("method"))
    variant = first_present(identity["variant"], args.get("variant_name"), summary.get("variant"), method)
    seed = first_present(status.get("seed"), args.get("seed"), summary.get("seed"), identity["seed"])
    stage = first_present(identity["stage"], args.get("stage"))
    sweep = first_present(identity["sweep"], args.get("sweep"))

    metrics, score_source = extract_run_metrics(run_dir)
    markers = {path.name for path in run_dir.iterdir() if path.is_file() and path.name in RUN_MARKERS}
    status_value = first_present(status.get("status"), status.get("state"))
    if not status_value:
        status_value = "completed_no_status_json" if metrics else "no_score_marker"

    row = row_base()
    row.update(
        {
            "record_kind": "run_artifact",
            "row_granularity": "run",
            "experiment_group": identity["experiment_group"],
            "dataset_group": identity["dataset_group"],
            "stage": stage,
            "sweep": sweep,
            "dataset": dataset,
            "method": method,
            "variant": variant,
            "model_key": model_key(method, variant),
            "seed": seed,
            "status": status_value,
            "return_code": first_present(status.get("return_code"), status.get("returncode"), status.get("exit_code")),
            "error": first_present(status.get("error"), status.get("reason")),
            "source_path": rel(run_dir),
            "run_dir": rel(run_dir),
            "metrics_path": rel(run_dir / "metrics.json") if (run_dir / "metrics.json").exists() else "",
            "summary_path": rel(run_dir / "summary.json") if (run_dir / "summary.json").exists() else "",
            "eval_path": rel(run_dir / "eval_fixed.csv") if (run_dir / "eval_fixed.csv").exists() else "",
            "args_path": rel(run_dir / "args.json") if (run_dir / "args.json").exists() else "",
            "status_path": rel(run_dir / "status.json") if (run_dir / "status.json").exists() else "",
            "commit_sha": first_present(status.get("commit_sha"), status.get("git_commit"), args.get("git_commit")),
            "branch": first_present(status.get("branch"), args.get("branch")),
            "gpu": first_present(status.get("gpu"), args.get("gpu")),
            "elapsed_seconds": first_present(status.get("elapsed_seconds"), status.get("runtime_seconds")),
            "runtime_seconds": first_present(status.get("runtime_seconds"), status.get("elapsed_seconds")),
            "command": first_present(status.get("command")),
            "score_source": score_source,
            "extra_json": json.dumps({"markers": sorted(markers)}, sort_keys=True),
        }
    )
    for key in SCORE_COLUMNS:
        row[key] = clean_value(metrics.get(key, ""))
    if row["macro_f1"] == "" and row["f1_macro"] != "":
        row["macro_f1"] = row["f1_macro"]
    return row


def csv_has_result_shape(path: Path, columns: set[str]) -> bool:
    lower_name = path.name.lower()
    if lower_name in {"per_cell_type_metrics.csv", "confusion_matrix_raw.csv", "confusion_matrix_mapped.csv"}:
        return False
    has_dataset = "dataset" in columns
    has_score = bool(columns & set(SCORE_COLUMNS + MEAN_SCORE_COLUMNS + ["ari", "acc", "nmi"]))
    has_method = bool(columns & {"method", "method_key", "variant", "mix_mode", "run", "label"})
    return has_dataset and has_score and has_method


def iter_summary_csv_rows() -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for path in sorted(RESULTS_ROOT.rglob("*.csv"), key=lambda p: rel_results(p)):
        parts = tuple(Path(rel_results(path)).parts)
        if parts and parts[0] == "_index":
            continue
        if path.name == "eval_fixed.csv":
            continue
        try:
            header = pd.read_csv(path, nrows=0)
        except Exception:
            continue
        columns = set(header.columns)
        if not csv_has_result_shape(path, columns):
            continue
        try:
            df = pd.read_csv(path)
        except Exception:
            continue
        if df.empty:
            continue
        experiment_group = experiment_group_from_rel(parts)
        for _, series in df.iterrows():
            raw = series.to_dict()
            dataset = first_present(raw.get("dataset"))
            method = first_present(raw.get("method"), raw.get("method_key"), raw.get("mix_mode"), raw.get("run"), raw.get("label"))
            variant = first_present(raw.get("variant"), raw.get("method_key"), raw.get("mix_mode"), raw.get("run"), method)
            seed = first_present(raw.get("seed"))

            if seed != "":
                granularity = "run_summary_csv"
            elif first_present(raw.get("n_success"), raw.get("n_total")) != "":
                granularity = "dataset_model_summary_csv"
            else:
                granularity = "summary_csv"

            row = row_base()
            row.update(
                {
                    "record_kind": "summary_csv",
                    "row_granularity": granularity,
                    "experiment_group": experiment_group,
                    "dataset_group": first_present(raw.get("dataset_group")),
                    "stage": first_present(raw.get("stage")),
                    "sweep": first_present(raw.get("sweep")),
                    "dataset": dataset,
                    "method": method,
                    "variant": variant,
                    "model_key": model_key(method, variant),
                    "seed": seed,
                    "status": first_present(raw.get("status"), "summary_only"),
                    "n_success": first_present(raw.get("n_success")),
                    "n_total": first_present(raw.get("n_total")),
                    "status_summary": first_present(raw.get("status_summary")),
                    "authenticity": first_present(raw.get("authenticity")),
                    "substitute_model_used": first_present(raw.get("substitute_model_used")),
                    "source_path": rel(path),
                    "run_dir": first_present(raw.get("run_dir"), raw.get("run_path"), raw.get("path"), raw.get("save_dir")),
                    "runtime_seconds": first_present(raw.get("runtime_seconds")),
                    "gpu": first_present(raw.get("gpu")),
                    "commit_sha": first_present(raw.get("commit_sha")),
                    "branch": first_present(raw.get("branch")),
                    "score_source": path.name,
                }
            )
            for key in SCORE_COLUMNS + MEAN_SCORE_COLUMNS:
                row[key] = clean_value(raw.get(key, ""))
            if row["macro_f1"] == "" and row["f1_macro"] != "":
                row["macro_f1"] = row["f1_macro"]
            row["extra_json"] = json.dumps(
                {
                    key: clean_value(value)
                    for key, value in raw.items()
                    if key not in set(OUTPUT_COLUMNS) and clean_value(value) != ""
                },
                ensure_ascii=True,
                sort_keys=True,
                default=str,
            )
            rows.append(row)
    return rows


def norm_key(row: dict[str, Any], include_score: bool = False) -> tuple[Any, ...]:
    score = ""
    if include_score:
        for key in ("ari", "ari_mean", "nmi", "acc"):
            value = row.get(key, "")
            if value != "":
                try:
                    score = round(float(value), 8)
                except Exception:
                    score = str(value)
                break
    return (
        row.get("experiment_group", ""),
        row.get("dataset", ""),
        row.get("model_key", ""),
        str(row.get("seed", "")),
        score,
    )


def global_norm_key(row: dict[str, Any], include_score: bool = False) -> tuple[Any, ...]:
    score = ""
    if include_score:
        for key in ("ari", "ari_mean", "nmi", "acc"):
            value = row.get(key, "")
            if value != "":
                try:
                    score = round(float(value), 8)
                except Exception:
                    score = str(value)
                break
    return (
        row.get("dataset", ""),
        row.get("model_key", ""),
        str(row.get("seed", "")),
        score,
    )


def mark_primary(rows: list[dict[str, Any]]) -> None:
    run_keys = {norm_key(row, include_score=False) for row in rows if row["record_kind"] == "run_artifact"}
    run_score_keys = {norm_key(row, include_score=True) for row in rows if row["record_kind"] == "run_artifact"}
    global_run_keys = {global_norm_key(row, include_score=False) for row in rows if row["record_kind"] == "run_artifact"}
    global_run_score_keys = {global_norm_key(row, include_score=True) for row in rows if row["record_kind"] == "run_artifact"}

    for row in rows:
        row["is_primary"] = "true"
        row["duplicate_reason"] = ""
        if row["record_kind"] == "summary_csv":
            key = norm_key(row, include_score=False)
            score_key = norm_key(row, include_score=True)
            global_key = global_norm_key(row, include_score=False)
            global_score_key = global_norm_key(row, include_score=True)
            if row["row_granularity"] == "dataset_model_summary_csv":
                row["is_primary"] = "false"
                row["duplicate_reason"] = "aggregate_summary_kept_in_evidence_table"
            elif key in run_keys or score_key in run_score_keys:
                row["is_primary"] = "false"
                row["duplicate_reason"] = "covered_by_run_artifact"
            elif global_key in global_run_keys or global_score_key in global_run_score_keys:
                row["is_primary"] = "false"
                row["duplicate_reason"] = "covered_by_run_artifact_global_key"


def assign_record_ids(rows: list[dict[str, Any]]) -> None:
    for idx, row in enumerate(rows, start=1):
        row["record_id"] = f"R{idx:06d}"


def write_csv(path: Path, rows: list[dict[str, Any]], columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: clean_value(row.get(key, "")) for key in columns})


def numeric(value: Any) -> float | None:
    try:
        if value == "":
            return None
        out = float(value)
    except Exception:
        return None
    if math.isnan(out):
        return None
    return out


def build_dataset_model_summary(primary_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    groups: dict[tuple[str, str, str], list[dict[str, Any]]] = defaultdict(list)
    for row in primary_rows:
        if row["row_granularity"] != "run":
            continue
        key = (row["experiment_group"], row["dataset"], row["model_key"])
        groups[key].append(row)

    out: list[dict[str, Any]] = []
    for (experiment_group, dataset, model), group_rows in sorted(groups.items()):
        scores: dict[str, list[float]] = {key: [] for key in SCORE_COLUMNS}
        status_counts = Counter(as_text(row.get("status")) for row in group_rows)
        seeds = sorted({as_text(row.get("seed")) for row in group_rows if as_text(row.get("seed"))})
        for row in group_rows:
            for key in SCORE_COLUMNS:
                value = numeric(row.get(key))
                if value is not None:
                    scores[key].append(value)
        item: dict[str, Any] = {
            "experiment_group": experiment_group,
            "dataset": dataset,
            "model_key": model,
            "n_runs": len(group_rows),
            "n_scored": len(scores["ari"]) if scores["ari"] else max(len(v) for v in scores.values()),
            "seeds": ";".join(seeds),
            "status_summary": ";".join(f"{k}:{v}" for k, v in sorted(status_counts.items())),
        }
        for key in SCORE_COLUMNS:
            values = scores[key]
            if values:
                item[f"{key}_mean"] = sum(values) / len(values)
                item[f"{key}_min"] = min(values)
                item[f"{key}_max"] = max(values)
            else:
                item[f"{key}_mean"] = ""
                item[f"{key}_min"] = ""
                item[f"{key}_max"] = ""
        out.append(item)
    return out


def build_coverage(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    groups: dict[tuple[str, str], list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        if row["is_primary"] != "true":
            continue
        groups[(row["experiment_group"], row["dataset"])].append(row)
    out = []
    for (experiment_group, dataset), group_rows in sorted(groups.items()):
        models = sorted({as_text(row.get("model_key")) for row in group_rows if as_text(row.get("model_key"))})
        statuses = Counter(as_text(row.get("status")) for row in group_rows)
        out.append(
            {
                "experiment_group": experiment_group,
                "dataset": dataset,
                "n_primary_records": len(group_rows),
                "n_models": len(models),
                "models": ";".join(models),
                "status_summary": ";".join(f"{k}:{v}" for k, v in sorted(statuses.items())),
            }
        )
    return out


def write_readme(
    all_rows: list[dict[str, Any]],
    primary_rows: list[dict[str, Any]],
    summary_rows: list[dict[str, Any]],
    dataset_summary_rows: list[dict[str, Any]],
) -> None:
    kinds = Counter(row["record_kind"] for row in all_rows)
    statuses = Counter(as_text(row.get("status")) for row in primary_rows)
    groups = Counter(as_text(row.get("experiment_group")) for row in primary_rows)
    lines = [
        "# Results master table 2026-06-19",
        "",
        "Generated from `results/` with symlinks followed.",
        "",
        "## Files",
        "",
        "- `all_results_master_table.csv`: primary table. It contains one row per run artifact plus summary-only rows that are not covered by a run artifact.",
        "- `all_results_evidence_table.csv`: exhaustive evidence table. It keeps every compatible run artifact and every compatible summary CSV row, including duplicate aggregate rows.",
        "- `dataset_model_status_summary.csv`: primary run-level aggregation by experiment group, dataset, and model key.",
        "- `dataset_coverage_matrix.csv`: dataset-level coverage and status counts.",
        "- `build_results_master.py`: reproducible generator.",
        "",
        "## Counts",
        "",
        f"- total evidence rows: {len(all_rows)}",
        f"- primary master rows: {len(primary_rows)}",
        f"- summary CSV evidence rows: {len(summary_rows)}",
        f"- dataset/model summary rows: {len(dataset_summary_rows)}",
        "",
        "## Primary Status Counts",
        "",
        "| status | count |",
        "| --- | ---: |",
    ]
    for key, count in statuses.most_common():
        lines.append(f"| `{key}` | {count} |")
    lines.extend(["", "## Primary Experiment Groups", "", "| group | count |", "| --- | ---: |"])
    for key, count in groups.most_common():
        lines.append(f"| `{key}` | {count} |")
    lines.extend(["", "## Evidence Row Kinds", "", "| kind | count |", "| --- | ---: |"])
    for key, count in kinds.most_common():
        lines.append(f"| `{key}` | {count} |")
    lines.append("")
    (OUT_DIR / "README.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    run_rows = [extract_run_row(path) for path in iter_run_dirs()]
    summary_rows = iter_summary_csv_rows()
    all_rows = run_rows + summary_rows
    mark_primary(all_rows)
    assign_record_ids(all_rows)

    primary_rows = [row for row in all_rows if row["is_primary"] == "true"]
    dataset_summary_rows = build_dataset_model_summary(primary_rows)
    coverage_rows = build_coverage(primary_rows)

    write_csv(OUT_DIR / "all_results_evidence_table.csv", all_rows, OUTPUT_COLUMNS)
    write_csv(OUT_DIR / "all_results_master_table.csv", primary_rows, OUTPUT_COLUMNS)

    summary_columns = [
        "experiment_group",
        "dataset",
        "model_key",
        "n_runs",
        "n_scored",
        "seeds",
        "status_summary",
    ]
    for key in SCORE_COLUMNS:
        summary_columns.extend([f"{key}_mean", f"{key}_min", f"{key}_max"])
    write_csv(OUT_DIR / "dataset_model_status_summary.csv", dataset_summary_rows, summary_columns)
    write_csv(
        OUT_DIR / "dataset_coverage_matrix.csv",
        coverage_rows,
        ["experiment_group", "dataset", "n_primary_records", "n_models", "models", "status_summary"],
    )
    write_readme(all_rows, primary_rows, summary_rows, dataset_summary_rows)

    print(f"run artifact rows: {len(run_rows)}")
    print(f"summary csv evidence rows: {len(summary_rows)}")
    print(f"primary master rows: {len(primary_rows)}")
    print(f"wrote: {OUT_DIR}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
