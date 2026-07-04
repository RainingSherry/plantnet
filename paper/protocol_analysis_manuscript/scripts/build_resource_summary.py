#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import math
from collections import defaultdict
from pathlib import Path
from statistics import mean, stdev


SCRIPT_DIR = Path(__file__).resolve().parent
MANUSCRIPT_ROOT = SCRIPT_DIR.parent
PROJECT_ROOT = MANUSCRIPT_ROOT.parents[1]

PHASE13_ROOT = PROJECT_ROOT / "results/CAAM_scMAE_correction/corruption_triad/formal"
PHASE14_ROOT = PROJECT_ROOT / "results/CAAM_scMAE_correction/advmask_triage/formal"
ATTENTION_ROOT = Path("/tmp/caam_attention_context_smoke/dev_20260626")

WALL_CLOCK_KEYS = (
    "wall_clock_seconds",
    "elapsed_seconds",
    "duration_seconds",
    "runtime_seconds",
    "train_seconds",
    "train_time_seconds",
)

STUDY_ORDER = {
    "phase13_corruption_triad": 0,
    "phase14_advmask_triage": 1,
    "attention_context_smoke": 2,
}
CONDITION_ORDER = {
    "scmae_shuffle": 0,
    "matched_donor": 1,
    "nonzero_aware_donor": 2,
    "control": 3,
    "advmask": 4,
    "axial": 5,
    "mlp_parammatched": 6,
}
DISPLAY_NAMES = {
    "phase13_corruption_triad": "Phase 13",
    "phase14_advmask_triage": "Phase 14",
    "attention_context_smoke": "Attention smoke",
    "scmae_shuffle": "scMAE shuffle",
    "matched_donor": "matched donor",
    "nonzero_aware_donor": "nonzero-aware donor",
    "control": "control",
    "advmask": "AdvMask",
    "axial": "Axial",
    "mlp_parammatched": "parameter-matched MLP",
}


def load_json(path: Path) -> dict:
    with open(path, "r", encoding="utf-8") as handle:
        return json.load(handle)


def parse_seed(token: str) -> int:
    if not token.startswith("seed"):
        raise ValueError(f"Expected seed token, got {token!r}")
    return int(token.removeprefix("seed"))


def parse_epochs(token: str) -> int:
    if not token.startswith("epochs"):
        raise ValueError(f"Expected epochs token, got {token!r}")
    return int(token.removeprefix("epochs"))


def metric_value(metrics_path: Path, group: str, metric: str) -> float:
    metrics = load_json(metrics_path)
    return float(metrics[group][metric])


def wall_clock_seconds(runtime: dict) -> float | None:
    present_keys = [key for key in WALL_CLOCK_KEYS if key in runtime]
    if not present_keys:
        return None
    if len(present_keys) > 1:
        raise ValueError(f"Multiple wall-clock keys found in runtime.json: {present_keys}")
    value = runtime[present_keys[0]]
    if value is None:
        return None
    return float(value)


def resource_row(
    *,
    study: str,
    dataset: str,
    corruption: str,
    condition: str,
    seed: int,
    epochs: int,
    run_dir: Path,
) -> dict:
    run_manifest = load_json(run_dir / "run_manifest.json")
    runtime = load_json(run_dir / "runtime.json")
    metrics_path = run_dir / "metrics.json"
    student_params = int(run_manifest["student_trainable_params"])
    generator_params = int(run_manifest["generator_trainable_params"])
    wall_seconds = wall_clock_seconds(runtime)
    return {
        "study": study,
        "study_label": DISPLAY_NAMES[study],
        "dataset": dataset,
        "corruption": corruption,
        "condition": condition,
        "condition_label": DISPLAY_NAMES[condition],
        "seed": seed,
        "epochs": epochs,
        "student_trainable_params": student_params,
        "generator_trainable_params": generator_params,
        "total_trainable_params": student_params + generator_params,
        "kmeans_known_k.ari": metric_value(metrics_path, "kmeans_known_k", "ari"),
        "leiden_fixed.ari": metric_value(metrics_path, "leiden_fixed", "ari"),
        "logical_device": str(runtime["logical_device"]),
        "physical_gpu": "" if runtime["physical_gpu"] is None else str(runtime["physical_gpu"]),
        "cuda_visible_devices": str(runtime["cuda_visible_devices"]),
        "amp": bool(runtime["amp"]),
        "num_workers": int(runtime["num_workers"]),
        "wall_clock_seconds": "" if wall_seconds is None else wall_seconds,
        "wall_clock_recorded": wall_seconds is not None,
        "run_dir": str(run_dir),
    }


def discover_phase13(root: Path) -> list[dict]:
    rows = []
    for manifest_path in sorted(root.glob("*__*__seed*__epochs*/run_manifest.json")):
        run_dir = manifest_path.parent
        parts = run_dir.name.split("__")
        if len(parts) != 4:
            raise ValueError(f"Unexpected Phase 13 run directory name: {run_dir.name}")
        dataset, corruption, seed_token, epoch_token = parts
        rows.append(
            resource_row(
                study="phase13_corruption_triad",
                dataset=dataset,
                corruption=corruption,
                condition=corruption,
                seed=parse_seed(seed_token),
                epochs=parse_epochs(epoch_token),
                run_dir=run_dir,
            )
        )
    return rows


def discover_phase14(root: Path) -> list[dict]:
    rows = []
    for manifest_path in sorted(root.glob("*__*__*__seed*__epochs*/run_manifest.json")):
        run_dir = manifest_path.parent
        parts = run_dir.name.split("__")
        if len(parts) != 5:
            raise ValueError(f"Unexpected Phase 14 run directory name: {run_dir.name}")
        dataset, corruption, variant, seed_token, epoch_token = parts
        rows.append(
            resource_row(
                study="phase14_advmask_triage",
                dataset=dataset,
                corruption=corruption,
                condition=variant,
                seed=parse_seed(seed_token),
                epochs=parse_epochs(epoch_token),
                run_dir=run_dir,
            )
        )
    return rows


def discover_attention(root: Path) -> list[dict]:
    rows = []
    for manifest_path in sorted(root.glob("*__*__*__seed*__epochs*/run_manifest.json")):
        run_dir = manifest_path.parent
        parts = run_dir.name.split("__")
        if len(parts) != 5:
            raise ValueError(f"Unexpected attention smoke run directory name: {run_dir.name}")
        dataset, corruption, role, seed_token, epoch_token = parts
        rows.append(
            resource_row(
                study="attention_context_smoke",
                dataset=dataset,
                corruption=corruption,
                condition=role,
                seed=parse_seed(seed_token),
                epochs=parse_epochs(epoch_token),
                run_dir=run_dir,
            )
        )
    return rows


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        raise ValueError(f"No rows to write for {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(rows[0].keys())
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def numeric_summary(values: list[float]) -> tuple[float, float, float, float]:
    return (
        min(values),
        max(values),
        mean(values),
        stdev(values) if len(values) > 1 else 0.0,
    )


def unique_join(values: list[object]) -> str:
    unique_values = sorted({str(value) for value in values if str(value) != ""})
    return ";".join(unique_values)


def aggregate(rows: list[dict], keys: tuple[str, ...]) -> list[dict]:
    grouped = defaultdict(list)
    for row in rows:
        grouped[tuple(row[key] for key in keys)].append(row)

    out = []
    for key_values, group_rows in grouped.items():
        row = {key: value for key, value in zip(keys, key_values)}
        row["study_label"] = DISPLAY_NAMES[str(row["study"])]
        row["condition_label"] = DISPLAY_NAMES[str(row["condition"])]
        row["n_runs"] = len(group_rows)
        row["n_datasets"] = len({item["dataset"] for item in group_rows})
        row["epochs"] = unique_join([item["epochs"] for item in group_rows])
        for field in (
            "student_trainable_params",
            "generator_trainable_params",
            "total_trainable_params",
            "kmeans_known_k.ari",
            "leiden_fixed.ari",
        ):
            values = [float(item[field]) for item in group_rows]
            min_value, max_value, mean_value, std_value = numeric_summary(values)
            row[f"{field}.min"] = min_value
            row[f"{field}.max"] = max_value
            row[f"{field}.mean"] = mean_value
            row[f"{field}.std"] = std_value
        recorded_wall_times = [float(item["wall_clock_seconds"]) for item in group_rows if item["wall_clock_recorded"]]
        row["wall_clock_recorded_runs"] = len(recorded_wall_times)
        row["wall_clock_status"] = f"{len(recorded_wall_times)}/{len(group_rows)} recorded"
        if recorded_wall_times:
            min_value, max_value, mean_value, std_value = numeric_summary(recorded_wall_times)
            row["wall_clock_seconds.min"] = min_value
            row["wall_clock_seconds.max"] = max_value
            row["wall_clock_seconds.mean"] = mean_value
            row["wall_clock_seconds.std"] = std_value
        else:
            row["wall_clock_seconds.min"] = ""
            row["wall_clock_seconds.max"] = ""
            row["wall_clock_seconds.mean"] = ""
            row["wall_clock_seconds.std"] = ""
        row["logical_devices"] = unique_join([item["logical_device"] for item in group_rows])
        row["physical_gpus"] = unique_join([item["physical_gpu"] for item in group_rows])
        row["cuda_visible_devices"] = unique_join([item["cuda_visible_devices"] for item in group_rows])
        row["amp_values"] = unique_join([item["amp"] for item in group_rows])
        row["num_workers_values"] = unique_join([item["num_workers"] for item in group_rows])
        out.append(row)

    out.sort(key=lambda row: (STUDY_ORDER[str(row["study"])], CONDITION_ORDER[str(row["condition"])]))
    return out


def latex_escape(value: object) -> str:
    return str(value).replace("_", r"\_")


def fmt_float(value: object, digits: int = 6) -> str:
    if value == "":
        return "--"
    number = float(value)
    if not math.isfinite(number):
        return "--"
    return f"{number:.{digits}f}"


def fmt_param_range(row: dict, field: str) -> str:
    min_value = int(row[f"{field}.min"])
    max_value = int(row[f"{field}.max"])
    if min_value == max_value:
        return f"{min_value:,}"
    return f"{min_value:,}--{max_value:,}"


def fmt_wall_time(row: dict) -> str:
    if row["wall_clock_recorded_runs"] == 0:
        return "not recorded"
    return fmt_float(row["wall_clock_seconds.mean"], digits=2)


def write_latex_table(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        r"\begin{table}[ht]",
        r"\centering",
        r"\caption{Trainable-parameter and recorded-runtime summary from existing development artifacts. Wall-clock and memory values were not recorded by the current source runtime manifests.}",
        r"\label{tab:generated-resource-summary}",
        r"\scriptsize",
        r"\begin{tabular}{llrrrrrl}",
        r"\toprule",
        r"Study & Condition & Runs & Student params & Gen. params & Known-K ARI & Leiden ARI & Wall time \\",
        r"\midrule",
    ]
    for row in rows:
        lines.append(
            " & ".join(
                [
                    latex_escape(row["study_label"]),
                    latex_escape(row["condition_label"]),
                    str(row["n_runs"]),
                    fmt_param_range(row, "student_trainable_params"),
                    fmt_param_range(row, "generator_trainable_params"),
                    fmt_float(row["kmeans_known_k.ari.mean"]),
                    fmt_float(row["leiden_fixed.ari.mean"]),
                    latex_escape(fmt_wall_time(row)),
                ]
            )
            + r" \\"
        )
    lines.extend([r"\bottomrule", r"\end{tabular}", r"\end{table}", ""])
    path.write_text("\n".join(lines), encoding="utf-8")


def build_outputs(args: argparse.Namespace) -> None:
    phase13_rows = discover_phase13(args.phase13_root)
    phase14_rows = discover_phase14(args.phase14_root)
    attention_rows = discover_attention(args.attention_root)

    if len(phase13_rows) != 27:
        raise ValueError(f"Expected 27 Phase 13 runs, found {len(phase13_rows)}")
    if len(phase14_rows) != 18:
        raise ValueError(f"Expected 18 Phase 14 runs, found {len(phase14_rows)}")
    if len(attention_rows) != 9:
        raise ValueError(f"Expected 9 attention smoke runs, found {len(attention_rows)}")

    rows = phase13_rows + phase14_rows + attention_rows
    rows.sort(
        key=lambda row: (
            STUDY_ORDER[row["study"]],
            row["dataset"],
            CONDITION_ORDER[row["condition"]],
            row["seed"],
        )
    )

    output_dir = args.output_dir
    data_dir = output_dir / "data"
    tables_dir = output_dir / "tables"
    data_dir.mkdir(parents=True, exist_ok=True)
    tables_dir.mkdir(parents=True, exist_ok=True)

    global_summary = aggregate(rows, ("study", "condition"))
    dataset_summary = aggregate(rows, ("study", "dataset", "condition"))

    write_csv(data_dir / "resource_runs.csv", rows)
    write_csv(data_dir / "resource_summary.csv", global_summary)
    write_csv(data_dir / "resource_summary_by_dataset.csv", dataset_summary)
    write_latex_table(tables_dir / "resource_summary.tex", global_summary)

    manifest = {
        "phase13_root": str(args.phase13_root),
        "phase14_root": str(args.phase14_root),
        "attention_root": str(args.attention_root),
        "phase13_runs": len(phase13_rows),
        "phase14_runs": len(phase14_rows),
        "attention_runs": len(attention_rows),
        "claim_scope": "development evidence only; no validation or sealed test data",
        "resource_scope": "trainable parameters and device metadata from existing manifests",
        "wall_clock_seconds_status": "not recorded by current source runtime.json files",
        "memory_status": "not recorded by current source manifests",
        "known_k_metric": "kmeans_known_k.ari",
        "non_oracle_metric": "leiden_fixed.ari",
    }
    (output_dir / "artifact_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description="Build resource summary tables from existing CAAM development artifacts.")
    parser.add_argument("--phase13-root", type=Path, default=PHASE13_ROOT)
    parser.add_argument("--phase14-root", type=Path, default=PHASE14_ROOT)
    parser.add_argument("--attention-root", type=Path, default=ATTENTION_ROOT)
    parser.add_argument("--output-dir", type=Path, default=MANUSCRIPT_ROOT / "generated/resource_summary")
    args = parser.parse_args()
    build_outputs(args)
    print(f"Wrote resource summary assets to {args.output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
