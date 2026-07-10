from __future__ import annotations

import argparse
import json
import subprocess
from itertools import product
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from methods.DeepLearning import scMAE_family as family

from .config import PAPER_ROOT, REMOTE_RESULT_ROOT, SEEDS
from .downstream_orchestrate import _parse_sha256sums
from .io_utils import sha256_file, utc_now, write_checksum_manifest, write_json
from .prepare_full_label_sensitivity import VERSION
from .remote_store import RemoteStore
from .run_full_label_sensitivity import SENSITIVITY_VARIANTS
from .secondary_aggregate import (
    CONFIRMATORY_METRICS, ITERATIONS, RNG_SEED, TIE_TOLERANCE,
    hierarchical_ci, write_frame,
)
from .statistics import CONTRASTS, holm_adjust, sign_flip_pvalue


ALL_METRICS = ("ari", "nmi", "acc", "f1_macro", "fmi", "v_measure", "homogeneity", "completeness")
SELECTED_FILES = ("run_metadata.json", "metrics.json", "clusters.npz", "run.log", "COMPLETED")


def sync_index(root: Path) -> None:
    store = RemoteStore(); root.mkdir(parents=True, exist_ok=True)
    command = ["rsync", "-a", "--prune-empty-dirs", "-e", store.rsync_shell, "--include=*/"]
    command.extend(f"--include={name}" for name in (*SELECTED_FILES, "SHA256SUMS"))
    command.extend([
        "--exclude=*", f"{store.target}:{REMOTE_RESULT_ROOT}/runs/{VERSION}/", f"{root}/",
    ])
    subprocess.run(command, check=True)


def verify_selected(run_dir: Path) -> None:
    checksums = _parse_sha256sums(run_dir / "SHA256SUMS")
    for name in SELECTED_FILES:
        path = run_dir / name
        if name not in checksums or not path.is_file():
            raise FileNotFoundError(f"Missing checksummed sensitivity artifact: {path}")
        if sha256_file(path) != checksums[name]:
            raise ValueError(f"Sensitivity checksum mismatch: {path}")


def recompute(run_dir: Path, n_cells: int, n_labels: int) -> tuple[dict[str, float], int]:
    with np.load(run_dir / "clusters.npz", allow_pickle=False) as bundle:
        if set(bundle.files) != {"labels", "predicted", "mapped"}:
            raise ValueError("Unexpected sensitivity clusters.npz schema")
        labels, pred, mapped = bundle["labels"], bundle["predicted"], bundle["mapped"]
    for array in (labels, pred, mapped):
        if array.ndim != 1 or len(array) != n_cells or array.dtype.kind not in "iu":
            raise ValueError("Invalid sensitivity cluster array")
    if len(np.unique(labels)) != n_labels:
        raise ValueError("Sensitivity gold-label K differs from frozen full-label K")
    if not 1 <= len(np.unique(pred)) <= n_labels:
        raise ValueError("Sensitivity predicted cluster count is outside the known-K contract")
    recalculated, expected_mapped = family.compute_kmeans_metrics(labels, pred)
    if not np.array_equal(mapped.astype(np.int64), expected_mapped):
        raise ValueError("Sensitivity Hungarian mapping mismatch")
    return {metric: float(recalculated[metric]) for metric in ALL_METRICS}, int(recalculated["n_pred_clusters"])


def collect(index_root: Path, freeze_path: Path, manifest_path: Path) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    freeze = json.loads(freeze_path.read_text(encoding="utf-8"))
    manifest_rows = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest = {row["name"]: row for row in manifest_rows}
    if freeze["dataset_manifest_sha256"] != sha256_file(manifest_path):
        raise ValueError("Sensitivity manifest changed after freeze")
    expected_keys = set(product(manifest, SENSITIVITY_VARIANTS, SEEDS))
    rows, rejected, seen = [], [], set()
    for meta_path in sorted(index_root.rglob("run_metadata.json")):
        run_dir = meta_path.parent
        try:
            verify_selected(run_dir)
            metadata = json.loads(meta_path.read_text(encoding="utf-8"))
            key = (metadata["dataset"], metadata["variant"], int(metadata["seed"]))
            if key not in expected_keys or key in seen:
                raise ValueError(f"Unknown or duplicate sensitivity key: {key}")
            dataset = manifest[key[0]]
            if run_dir.name != metadata.get("config_hash") or not str(metadata.get("run_id", "")).endswith(
                f"--{metadata.get('config_hash')}"
            ):
                raise ValueError("Sensitivity run path/run_id/config mismatch")
            if metadata.get("sensitivity_freeze_hash") != freeze["freeze_hash"]:
                raise ValueError("Sensitivity freeze mismatch")
            if metadata.get("parent_primary_freeze_hash") != freeze["parent_primary_freeze_hash"]:
                raise ValueError("Sensitivity primary parent mismatch")
            if metadata.get("data", {}).get("sha256") != dataset["sha256"]:
                raise ValueError("Sensitivity data hash mismatch")
            if metadata.get("code_sha256") != freeze["code_sha256"]:
                raise ValueError("Sensitivity code identity mismatch")
            if metadata["protocol"].get("n_clusters") != int(dataset["n_labels"]):
                raise ValueError("Sensitivity K mismatch")
            if metadata["protocol"].get("protocol_version") != VERSION:
                raise ValueError("Sensitivity protocol version mismatch")
            if not str(metadata.get("remote_dir", "")).endswith(f"/{metadata['config_hash']}"):
                raise ValueError("Sensitivity remote/config hash mismatch")
            if (run_dir / "run.log").read_text(encoding="utf-8", errors="replace").splitlines()[0] != "Using device: cuda:0":
                raise ValueError("Sensitivity run did not declare cuda:0")
            values, n_pred = recompute(run_dir, int(dataset["n_cells"]), int(dataset["n_labels"]))
            recorded = metadata["validation"]["metrics"]
            metrics_json = json.loads((run_dir / "metrics.json").read_text(encoding="utf-8"))["kmeans_known_k"]
            for metric in ALL_METRICS:
                if not np.isclose(values[metric], float(recorded[metric]), rtol=1e-9, atol=1e-12):
                    raise ValueError(f"Sensitivity metadata metric mismatch: {metric}")
                if not np.isclose(values[metric], float(metrics_json[metric]), rtol=1e-9, atol=1e-12):
                    raise ValueError(f"Sensitivity metrics.json mismatch: {metric}")
            seen.add(key)
            rows.append({
                "run_id": metadata["run_id"], "dataset": key[0], "variant": key[1],
                "model_seed": key[2], "config_hash": metadata["config_hash"],
                "n_cells": int(dataset["n_cells"]), "n_labels": int(dataset["n_labels"]),
                "data_sha256": dataset["sha256"], "sensitivity_freeze_hash": freeze["freeze_hash"],
                "runtime_seconds": float(metadata["runtime_seconds"]), "n_pred_clusters": n_pred,
                "remote_dir": metadata["remote_dir"], **values,
            })
        except Exception as exc:
            rejected.append({"run_id": run_dir.name, "reason": str(exc)})
    missing = expected_keys - seen
    if rejected or missing:
        raise RuntimeError(
            f"Sensitivity evidence incomplete/invalid: accepted={len(seen)}/54 "
            f"rejected={len(rejected)} missing={len(missing)}"
        )
    return pd.DataFrame(rows).sort_values(["dataset", "variant", "model_seed"]), pd.DataFrame(rejected), freeze


def summarize(run_level: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    dataset_rows = []
    for (dataset, variant), group in run_level.groupby(["dataset", "variant"], sort=True):
        if len(group) != 3 or set(group["model_seed"]) != set(SEEDS):
            raise RuntimeError(f"Sensitivity seed grid incomplete: {dataset}/{variant}")
        row: dict[str, Any] = {"dataset": dataset, "variant": variant, "n_seeds": 3}
        for metric in ALL_METRICS:
            values = group[metric].to_numpy(float)
            row.update({f"{metric}_mean": values.mean(), f"{metric}_sd": values.std(ddof=1)})
        dataset_rows.append(row)
    dataset_level = pd.DataFrame(dataset_rows)
    overall_rows = []
    for variant, group in dataset_level.groupby("variant", sort=True):
        row = {"variant": variant, "n_datasets": len(group)}
        for metric in ALL_METRICS:
            values = group[f"{metric}_mean"].to_numpy(float)
            row.update({f"{metric}_mean": values.mean(), f"{metric}_sd": values.std(ddof=1)})
        overall_rows.append(row)
    return dataset_level, pd.DataFrame(overall_rows)


def contrast_table(run_level: pd.DataFrame, dataset_level: pd.DataFrame) -> pd.DataFrame:
    rng = np.random.default_rng(RNG_SEED); rows = []
    for metric in CONFIRMATORY_METRICS:
        pivot = dataset_level.pivot(index="dataset", columns="variant", values=f"{metric}_mean")
        family_rows = []
        for left, right, name in CONTRASTS:
            delta = (pivot[left] - pivot[right]).to_numpy(float)
            lo, hi = hierarchical_ci(run_level, left, right, metric, rng)
            family_rows.append({
                "metric": metric, "contrast": name, "n_datasets": len(delta),
                "mean_delta": delta.mean(), "median_delta": np.median(delta),
                "ci95_low": lo, "ci95_high": hi,
                "wins": int((delta > TIE_TOLERANCE).sum()),
                "ties": int((np.abs(delta) <= TIE_TOLERANCE).sum()),
                "losses": int((delta < -TIE_TOLERANCE).sum()),
                "permutation_p": sign_flip_pvalue(delta, ITERATIONS, rng),
            })
        for row, adjusted in zip(family_rows, holm_adjust([row["permutation_p"] for row in family_rows])):
            row["holm_p"] = adjusted; rows.append(row)
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description="Aggregate independently verified full-label sensitivity runs")
    parser.add_argument("--index-root", type=Path, default=PAPER_ROOT / f".staging/remote_{VERSION}_index")
    parser.add_argument("--freeze", type=Path, default=PAPER_ROOT / f"experiments/{VERSION}/source_freeze.json")
    parser.add_argument("--manifest", type=Path, default=PAPER_ROOT / f"manifests/{VERSION}/dataset_manifest.json")
    parser.add_argument("--output-dir", type=Path, default=PAPER_ROOT / f"experiments/{VERSION}/aggregate")
    parser.add_argument("--no-sync", action="store_true")
    args = parser.parse_args()
    if not args.no_sync:
        sync_index(args.index_root)
    run_level, rejected, freeze = collect(args.index_root, args.freeze, args.manifest)
    dataset_level, overall = summarize(run_level)
    contrasts = contrast_table(run_level, dataset_level)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for frame, name in ((run_level, "run_metrics"), (dataset_level, "dataset_variant_metrics"),
                        (overall, "variant_overall_metrics"), (contrasts, "contrasts"),
                        (rejected, "rejected_runs")):
        write_frame(frame, args.output_dir / name)
    write_json(args.output_dir / "aggregation_metadata.json", {
        "generated_utc": utc_now(), "version": VERSION, "accepted_runs": len(run_level),
        "expected_runs": 54, "freeze_hash": freeze["freeze_hash"],
        "dataset_is_independent_unit": True, "no_best_seed_selection": True,
        "bootstrap_iterations": ITERATIONS, "rng_seed": RNG_SEED,
        "aggregator_sha256": sha256_file(Path(__file__).resolve()),
    })
    write_checksum_manifest(args.output_dir)
    print(f"sensitivity_runs={len(run_level)} datasets={run_level.dataset.nunique()} contrasts={len(contrasts)}")


if __name__ == "__main__":
    main()
