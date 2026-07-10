from __future__ import annotations

import argparse
import json
import subprocess
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from .config import DATASETS, PAPER_ROOT, PROTOCOL_VERSION, REMOTE_RESULT_ROOT, SEEDS, SPLIT_SEEDS, VARIANTS
from .downstream_orchestrate import LABEL_FRACTIONS, _parse_sha256sums
from .io_utils import sha256_file
from .remote_store import RemoteStore
from .statistics import CONTRASTS, bootstrap_ci, holm_adjust, sign_flip_pvalue


def sync_metadata(local_root: Path) -> None:
    store = RemoteStore()
    local_root.mkdir(parents=True, exist_ok=True)
    command = [
        "rsync", "-a", "--prune-empty-dirs", "-e", store.rsync_shell,
        "--include=*/", "--include=downstream_summary.json",
        "--include=downstream_metadata.json", "--include=COMPLETED",
        "--include=SHA256SUMS", "--exclude=*",
        f"{store.target}:{REMOTE_RESULT_ROOT}/downstream/{PROTOCOL_VERSION}/",
        f"{local_root}/",
    ]
    subprocess.run(command, check=True)


def verify_selected_files(run_dir: Path) -> None:
    checksums = _parse_sha256sums(run_dir / "SHA256SUMS")
    for name in ("downstream_summary.json", "downstream_metadata.json", "COMPLETED"):
        if name not in checksums or not (run_dir / name).is_file():
            raise FileNotFoundError(f"Missing checksummed downstream artifact: {run_dir / name}")
        if sha256_file(run_dir / name) != checksums[name]:
            raise ValueError(f"Downstream checksum mismatch: {run_dir / name}")


def scalar_metrics(payload: dict[str, Any]) -> dict[str, float]:
    return {
        key: float(value) for key, value in payload.items()
        if isinstance(value, (int, float)) and not isinstance(value, bool)
    }


def collect(index_root: Path, run_master: Path) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    formal = pd.read_csv(run_master)
    formal = formal[formal["execution_mode"] == "formal"].copy()
    expected = len(DATASETS) * len(VARIANTS) * len(SEEDS)
    if len(formal) != expected or len(formal[["dataset", "variant", "seed"]].drop_duplicates()) != expected:
        raise RuntimeError("Downstream aggregation requires the exact complete 108-run frozen matrix")
    lookup = formal.set_index("run_id")[["dataset", "variant", "seed", "config_hash"]].to_dict("index")
    marker_rows: list[dict] = []
    probe_rows: list[dict] = []
    recall_rows: list[dict] = []
    seen: set[str] = set()
    for summary_path in sorted(index_root.rglob("downstream_summary.json")):
        run_dir = summary_path.parent
        verify_selected_files(run_dir)
        metadata = json.loads((run_dir / "downstream_metadata.json").read_text(encoding="utf-8"))
        run_id = str(metadata["run_id"])
        if run_id not in lookup or run_id in seen:
            raise ValueError(f"Unknown or duplicate downstream run_id: {run_id}")
        seen.add(run_id)
        key = lookup[run_id]
        if (metadata.get("dataset"), metadata.get("variant"), int(metadata.get("model_seed"))) != (
            key["dataset"], key["variant"], int(key["seed"])
        ):
            raise ValueError(f"Downstream metadata disagrees with frozen run master: {run_id}")
        summary = json.loads(summary_path.read_text(encoding="utf-8"))
        if [int(row["split_seed"]) for row in summary["marker_splits"]] != list(SPLIT_SEEDS):
            raise ValueError(f"Marker split grid changed: {run_id}")
        probe_keys = [(float(row["label_fraction"]), int(row["split_seed"])) for row in summary["linear_probe"]]
        if probe_keys != [(fraction, seed) for fraction in LABEL_FRACTIONS for seed in SPLIT_SEEDS]:
            raise ValueError(f"Probe split grid changed: {run_id}")
        base = {"run_id": run_id, "dataset": key["dataset"], "variant": key["variant"], "model_seed": int(key["seed"])}
        for row in summary["marker_splits"]:
            split = {**base, "split_seed": int(row["split_seed"])}
            marker_rows.append({
                **split,
                **{f"recovery_{k}": v for k, v in scalar_metrics(row["marker_recovery"]).items()},
                **{f"annotation_{k}": v for k, v in scalar_metrics(row["marker_annotation"]).items()},
                **{f"oracle_{k}": v for k, v in scalar_metrics(row["oracle_upper_bound"]).items()},
            })
            for scope in ("marker_annotation", "oracle_upper_bound"):
                for cell_type, recall in row[scope].get("per_class_recall", {}).items():
                    recall_rows.append({**split, "scope": scope, "cell_type": cell_type, "recall": float(recall)})
        for row in summary["linear_probe"]:
            probe = {
                **base, "split_seed": int(row["split_seed"]),
                "label_fraction": float(row["label_fraction"]),
                **{f"probe_{k}": v for k, v in scalar_metrics(row).items() if k not in {"split_seed", "label_fraction"}},
            }
            probe_rows.append(probe)
            for cell_type, recall in row.get("per_class_recall", {}).items():
                recall_rows.append({
                    **base, "split_seed": int(row["split_seed"]), "scope": "linear_probe",
                    "label_fraction": float(row["label_fraction"]), "cell_type": cell_type, "recall": float(recall),
                })
    if seen != set(lookup):
        raise RuntimeError(f"Downstream results incomplete: found {len(seen)}/{expected} formal run IDs")
    return pd.DataFrame(marker_rows), pd.DataFrame(probe_rows), pd.DataFrame(recall_rows)


def dataset_contrasts(dataset_level: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict] = []
    rng = np.random.default_rng(20260710)
    id_columns = {"dataset", "variant", "task", "label_fraction"}
    for group_key, group in dataset_level.groupby(["task", "label_fraction"], dropna=False):
        task, fraction = group_key
        for metric in [
            column for column in group.columns
            if column not in id_columns and group[column].notna().any()
        ]:
            pivot = group.pivot(index="dataset", columns="variant", values=metric)
            metric_rows = []
            for left, right, name in CONTRASTS:
                paired = pivot[[left, right]].dropna() if {left, right}.issubset(pivot.columns) else pd.DataFrame()
                delta = paired[left].to_numpy() - paired[right].to_numpy() if not paired.empty else np.array([])
                lo, hi = bootstrap_ci(delta, 10000, rng)
                metric_rows.append({
                    "task": task, "label_fraction": fraction, "metric": metric, "contrast": name,
                    "n_datasets": len(delta), "mean_delta": float(delta.mean()) if len(delta) else np.nan,
                    "ci95_low": lo, "ci95_high": hi, "wins": int((delta > 0).sum()),
                    "ties": int((delta == 0).sum()), "losses": int((delta < 0).sum()),
                    "permutation_p": sign_flip_pvalue(delta, 10000, rng),
                })
            adjusted = holm_adjust([row["permutation_p"] for row in metric_rows])
            for row, value in zip(metric_rows, adjusted):
                row["holm_p"] = value
                rows.append(row)
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description="Aggregate complete downstream results without pseudoreplication")
    parser.add_argument("--index-root", type=Path, default=PAPER_ROOT / ".staging/remote_downstream_index")
    parser.add_argument("--run-master", type=Path, default=PAPER_ROOT / "experiments/protocol_v1/run_master.csv")
    parser.add_argument("--output-dir", type=Path, default=PAPER_ROOT / "experiments/downstream_v1")
    parser.add_argument("--no-sync", action="store_true")
    args = parser.parse_args()
    if not args.no_sync:
        sync_metadata(args.index_root)
    marker, probe, recall = collect(args.index_root, args.run_master)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    marker.to_csv(args.output_dir / "marker_split_metrics.csv", index=False)
    probe.to_csv(args.output_dir / "probe_split_metrics.csv", index=False)
    recall.to_csv(args.output_dir / "per_class_recall.csv", index=False)
    marker_run = marker.drop(columns=["split_seed"]).groupby(
        ["run_id", "dataset", "variant", "model_seed"], as_index=False
    ).mean(numeric_only=True)
    probe_run = probe.drop(columns=["split_seed"]).groupby(
        ["run_id", "dataset", "variant", "model_seed", "label_fraction"], as_index=False
    ).mean(numeric_only=True)
    marker_dataset = marker_run.drop(columns=["run_id", "model_seed"]).groupby(
        ["dataset", "variant"], as_index=False
    ).mean(numeric_only=True).assign(task="marker", label_fraction=np.nan)
    probe_dataset = probe_run.drop(columns=["run_id", "model_seed"]).groupby(
        ["dataset", "variant", "label_fraction"], as_index=False
    ).mean(numeric_only=True).assign(task="linear_probe")
    dataset_level = pd.concat([marker_dataset, probe_dataset], ignore_index=True, sort=False)
    marker_run.to_csv(args.output_dir / "marker_run_metrics.csv", index=False)
    probe_run.to_csv(args.output_dir / "probe_run_metrics.csv", index=False)
    dataset_level.to_csv(args.output_dir / "dataset_variant_metrics.csv", index=False)
    dataset_contrasts(dataset_level).to_csv(args.output_dir / "downstream_contrasts.csv", index=False)
    print(f"marker_splits={len(marker)} probe_splits={len(probe)} run_ids={marker['run_id'].nunique()}")


if __name__ == "__main__":
    main()
