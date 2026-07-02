#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import re
from pathlib import Path
from statistics import mean
from typing import Any


ROLE_BY_METHOD = {
    "caam_scmae_control": "control",
    "caam_scmae_axial": "axial",
    "caam_scmae_advmask": "advmask",
    "caam_scmae_full": "full",
    "caam_scmae_mlp_parammatched": "mlp_parammatched",
}

ROLE_BY_PREFIX = {
    "control": "control",
    "axial": "axial",
    "advmask": "advmask",
    "full": "full",
    "mlp_parammatched": "mlp_parammatched",
}

REQUIRED_ROLES = ("control", "axial", "advmask", "full", "mlp_parammatched")
ROLE_TO_FACTORIAL = {
    "control": "Y00",
    "axial": "Y10",
    "advmask": "Y01",
    "full": "Y11",
}

METRIC_KEYS = (
    "kmeans_known_k.acc",
    "kmeans_known_k.nmi",
    "kmeans_known_k.ari",
    "kmeans_known_k.f1_macro",
    "leiden_fixed.acc",
    "leiden_fixed.nmi",
    "leiden_fixed.ari",
    "leiden_fixed.f1_macro",
)
PRIMARY_METRIC = "kmeans_known_k.ari"


def load_json(path: Path) -> dict[str, Any]:
    with open(path, "r", encoding="utf-8") as handle:
        return json.load(handle)


def infer_role(run_dir: Path, artifact: dict[str, Any]) -> str | None:
    method = str(artifact.get("method", ""))
    if method in ROLE_BY_METHOD:
        return ROLE_BY_METHOD[method]
    prefix = run_dir.name.split("__", 1)[0]
    return ROLE_BY_PREFIX.get(prefix)


def infer_seed(run_dir: Path, artifact: dict[str, Any]) -> int | None:
    if artifact.get("seed") is not None:
        return int(artifact["seed"])
    match = re.search(r"__seed(\d+)", run_dir.name)
    return int(match.group(1)) if match else None


def flatten_metrics(metrics: dict[str, Any]) -> dict[str, float]:
    out: dict[str, float] = {}
    for metric_key in METRIC_KEYS:
        group, name = metric_key.split(".", 1)
        value = metrics.get(group, {}).get(name)
        out[metric_key] = float(value) if value is not None else float("nan")
    return out


def discover_runs(root: Path) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for metrics_path in sorted(root.rglob("metrics.json")):
        run_dir = metrics_path.parent
        artifact_path = run_dir / "artifact_manifest.json"
        artifact = load_json(artifact_path) if artifact_path.exists() else {}
        role = infer_role(run_dir, artifact)
        if role is None:
            continue
        metrics = load_json(metrics_path)
        row: dict[str, Any] = {
            "role": role,
            "factorial_symbol": ROLE_TO_FACTORIAL.get(role, ""),
            "seed": infer_seed(run_dir, artifact),
            "method": artifact.get("method", ""),
            "variant": artifact.get("variant", ""),
            "run_dir": str(run_dir),
        }
        row.update(flatten_metrics(metrics))
        rows.append(row)
    return rows


def mean_by_role(rows: list[dict[str, Any]]) -> dict[str, dict[str, float]]:
    grouped: dict[str, list[dict[str, Any]]] = {role: [] for role in REQUIRED_ROLES}
    for row in rows:
        if row["role"] in grouped:
            grouped[row["role"]].append(row)
    missing = [role for role, role_rows in grouped.items() if not role_rows]
    if missing:
        raise ValueError(f"Missing required ablation role(s): {', '.join(missing)}")
    return {
        role: {metric: mean(float(row[metric]) for row in role_rows) for metric in METRIC_KEYS}
        for role, role_rows in grouped.items()
    }


def build_interaction_report(role_means: dict[str, dict[str, float]], rows: list[dict[str, Any]]) -> dict[str, Any]:
    y00 = role_means["control"]
    y10 = role_means["axial"]
    y01 = role_means["advmask"]
    y11 = role_means["full"]
    param = role_means["mlp_parammatched"]
    delta_ab = {metric: y11[metric] - y10[metric] - y01[metric] + y00[metric] for metric in METRIC_KEYS}
    full_minus_axial = {metric: y11[metric] - y10[metric] for metric in METRIC_KEYS}
    full_minus_advmask = {metric: y11[metric] - y01[metric] for metric in METRIC_KEYS}
    full_minus_parammatched = {metric: y11[metric] - param[metric] for metric in METRIC_KEYS}

    claim_status_by_metric = {}
    for metric in METRIC_KEYS:
        positive = full_minus_axial[metric] > 0.0 and full_minus_advmask[metric] > 0.0 and delta_ab[metric] > 0.0
        claim_status_by_metric[metric] = "candidate_positive_interaction" if positive else "no_positive_interaction"

    return {
        "primary_metric": PRIMARY_METRIC,
        "claim_status": claim_status_by_metric[PRIMARY_METRIC],
        "claim_status_by_metric": claim_status_by_metric,
        "n_runs_by_role": {role: sum(1 for row in rows if row["role"] == role) for role in REQUIRED_ROLES},
        "Y00": y00,
        "Y10": y10,
        "Y01": y01,
        "Y11": y11,
        "parameter_matched_mlp": param,
        "delta_AB": delta_ab,
        "full_minus_axial": full_minus_axial,
        "full_minus_advmask": full_minus_advmask,
        "full_minus_parammatched_mlp": full_minus_parammatched,
        "note": "candidate_positive_interaction is not a synergy confirmation; paired CI is still required.",
    }


def write_summary_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    fieldnames = ["role", "factorial_symbol", "seed", "method", "variant", "run_dir", *METRIC_KEYS]
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in fieldnames})


def summarize(root: Path, output_dir: Path | None = None) -> tuple[Path, Path]:
    rows = discover_runs(root)
    role_means = mean_by_role(rows)
    report = build_interaction_report(role_means, rows)
    out_dir = output_dir or root
    out_dir.mkdir(parents=True, exist_ok=True)
    summary_path = out_dir / "ablation_summary.csv"
    report_path = out_dir / "interaction_report.json"
    write_summary_csv(summary_path, rows)
    with open(report_path, "w", encoding="utf-8") as handle:
        json.dump(report, handle, indent=2, sort_keys=True)
    return summary_path, report_path


def main() -> int:
    parser = argparse.ArgumentParser(description="Summarize CAAM-scMAE internal ablation runs.")
    parser.add_argument("ablation_root", type=Path, help="Directory produced by run_ablation.py.")
    parser.add_argument("--output_dir", type=Path, default=None, help="Output directory; defaults to ablation_root.")
    args = parser.parse_args()
    try:
        summary_path, report_path = summarize(args.ablation_root, args.output_dir)
    except Exception as exc:
        print(f"CAAM ablation summary failed: {exc}")
        return 1
    print(f"Wrote {summary_path}")
    print(f"Wrote {report_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
