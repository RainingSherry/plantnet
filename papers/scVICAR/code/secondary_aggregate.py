from __future__ import annotations

import argparse
import json
import subprocess
from datetime import datetime
from itertools import product
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from methods.DeepLearning import scMAE_family as family

from .config import DATASETS, PAPER_ROOT, PROTOCOL_VERSION, REMOTE_RESULT_ROOT, SEEDS, VARIANTS
from .downstream_orchestrate import _parse_sha256sums
from .io_utils import sha256_file, utc_now, write_checksum_manifest, write_json
from .remote_store import RemoteStore
from .statistics import CONTRASTS, holm_adjust, sign_flip_pvalue


SECONDARY_VERSION = "leiden_fixed_v1"
METRICS = ("ari", "nmi", "acc", "f1_macro", "fmi", "v_measure", "homogeneity", "completeness")
CONFIRMATORY_METRICS = ("ari", "nmi", "acc", "f1_macro")
PAYLOAD_FILES = {"COMPLETED", "leiden_clusters.npz", "metrics.json", "secondary_metadata.json"}
SOURCE_ARTIFACTS = {"clusters.npz", "embedding_float32.npz", "run_metadata.json"}
TIE_TOLERANCE = 1e-12
RNG_SEED = 20260710
ITERATIONS = 10000


def sync_index(local_root: Path) -> None:
    store = RemoteStore()
    local_root.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        [
            "rsync", "-a", "--prune-empty-dirs", "-e", store.rsync_shell,
            "--include=*/", "--include=metrics.json", "--include=secondary_metadata.json",
            "--include=leiden_clusters.npz", "--include=COMPLETED", "--include=SHA256SUMS",
            "--exclude=*",
            f"{store.target}:{REMOTE_RESULT_ROOT}/downstream/{SECONDARY_VERSION}/{PROTOCOL_VERSION}/",
            f"{local_root}/",
        ],
        check=True,
    )


def verify_manifest(run_dir: Path) -> dict[str, str]:
    manifest = run_dir / "SHA256SUMS"
    if not manifest.is_file():
        raise FileNotFoundError(f"Missing Leiden checksum manifest: {manifest}")
    checksums = _parse_sha256sums(manifest)
    if set(checksums) != PAYLOAD_FILES:
        raise ValueError(f"Unexpected Leiden payload manifest in {run_dir}: {sorted(checksums)}")
    actual = {path.name for path in run_dir.iterdir() if path.is_file()}
    if actual != PAYLOAD_FILES | {"SHA256SUMS"}:
        raise ValueError(f"Unexpected or missing Leiden files in {run_dir}: {sorted(actual)}")
    for name, expected in checksums.items():
        if sha256_file(run_dir / name) != expected:
            raise ValueError(f"Leiden checksum mismatch: {run_dir / name}")
    return checksums


def _assert_close(name: str, observed: Any, expected: Any, run_id: str) -> None:
    if isinstance(expected, float) and np.isnan(expected):
        if not (isinstance(observed, (int, float)) and np.isnan(float(observed))):
            raise ValueError(f"{name} must be unavailable (NaN): {run_id}")
    elif not np.isclose(float(observed), float(expected), rtol=1e-9, atol=1e-12):
        raise ValueError(f"Recomputed {name} disagrees with metrics.json: {run_id}")


def verify_arrays_and_metrics(
    run_dir: Path, metrics: dict[str, Any], expected_cells: int, expected_classes: int, run_id: str
) -> tuple[dict[str, float], int]:
    with np.load(run_dir / "leiden_clusters.npz", allow_pickle=False) as bundle:
        if set(bundle.files) != {"labels", "predicted", "mapped"}:
            raise ValueError(f"Unexpected Leiden NPZ arrays: {run_id}")
        labels = bundle["labels"]
        predicted = bundle["predicted"]
        mapped = bundle["mapped"]
    for name, array in (("labels", labels), ("predicted", predicted), ("mapped", mapped)):
        if array.ndim != 1 or len(array) != expected_cells or array.dtype.kind not in "iu":
            raise ValueError(f"Invalid Leiden {name} array: {run_id}")
    if len(np.unique(labels)) != expected_classes:
        raise ValueError(f"Gold label class count changed: {run_id}")
    pred_values = np.unique(predicted)
    if not np.array_equal(pred_values, np.arange(len(pred_values))):
        raise ValueError(f"Leiden predictions are not contiguous nonnegative codes: {run_id}")
    if not set(np.unique(mapped)).issubset(set(np.unique(labels))):
        raise ValueError(f"Mapped predictions leave the gold-label domain: {run_id}")
    recomputed, expected_mapped = family.compute_kmeans_metrics(labels, predicted)
    if not np.array_equal(mapped.astype(np.int64), expected_mapped):
        raise ValueError(f"Hungarian mapping changed: {run_id}")
    required = set(METRICS) | {
        "n_pred_clusters", "silhouette", "protocol", "cluster_method", "uses_known_k",
        "resolution", "n_neighbors",
    }
    missing = required - set(metrics)
    if missing:
        raise ValueError(f"Missing Leiden metrics {sorted(missing)}: {run_id}")
    for metric in (*METRICS, "silhouette"):
        _assert_close(metric, metrics[metric], recomputed[metric], run_id)
    if int(metrics["n_pred_clusters"]) != len(pred_values):
        raise ValueError(f"Leiden cluster count changed: {run_id}")
    if not (-1.0 <= float(metrics["ari"]) <= 1.0):
        raise ValueError(f"ARI outside [-1,1]: {run_id}")
    for metric in set(METRICS) - {"ari"}:
        if not (0.0 <= float(metrics[metric]) <= 1.0):
            raise ValueError(f"{metric} outside [0,1]: {run_id}")
    if (
        metrics["protocol"] != "fixed"
        or metrics["cluster_method"] != "leiden_fixed_resolution_1p0"
        or bool(metrics["uses_known_k"])
        or float(metrics["resolution"]) != 1.0
        or int(metrics["n_neighbors"]) != 15
    ):
        raise ValueError(f"Leiden metric protocol changed: {run_id}")
    return {metric: float(metrics[metric]) for metric in METRICS}, len(pred_values)


def _load_frozen_inputs(run_master: Path, freeze_path: Path, dataset_manifest: Path):
    formal = pd.read_csv(run_master)
    formal = formal[formal["execution_mode"] == "formal"].copy()
    expected_keys = set(product(DATASETS, VARIANTS, SEEDS))
    observed_keys = set(map(tuple, formal[["dataset", "variant", "seed"]].itertuples(index=False, name=None)))
    if len(formal) != len(expected_keys) or observed_keys != expected_keys or formal["run_id"].duplicated().any():
        raise RuntimeError("Leiden aggregation requires the exact frozen 6 x 6 x 3 Cartesian matrix")
    primary_freeze = json.loads(
        (PAPER_ROOT / f"experiments/{PROTOCOL_VERSION}/source_freeze.json").read_text(encoding="utf-8")
    )
    secondary_freeze = json.loads(freeze_path.read_text(encoding="utf-8"))
    if secondary_freeze["parent_primary_freeze_hash"] != primary_freeze["freeze_hash"]:
        raise ValueError("Secondary freeze does not bind the active primary freeze")
    if secondary_freeze["formal_run_master_sha256"] != sha256_file(run_master):
        raise ValueError("Formal run master changed after secondary freeze")
    if formal["code_hash"].nunique() != 1:
        raise ValueError("Primary formal runs contain multiple code hashes")
    manifest_rows = json.loads(dataset_manifest.read_text(encoding="utf-8"))
    datasets = {row["name"]: row for row in manifest_rows}
    if set(datasets) != set(DATASETS):
        raise ValueError("Canonical dataset manifest changed")
    for row in formal.itertuples(index=False):
        if row.data_sha256 != datasets[row.dataset]["sha256"]:
            raise ValueError(f"Primary data hash disagrees with canonical manifest: {row.run_id}")
        if not str(row.remote_dir).endswith(f"/{row.config_hash}") or not str(row.run_id).endswith(f"--{row.config_hash}"):
            raise ValueError(f"Primary run identity/config suffix changed: {row.run_id}")
    return formal, primary_freeze, secondary_freeze, datasets


def collect(
    index_root: Path, run_master: Path, freeze_path: Path, dataset_manifest: Path
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    formal, primary_freeze, secondary_freeze, datasets = _load_frozen_inputs(
        run_master, freeze_path, dataset_manifest
    )
    lookup = formal.set_index("run_id")[[
        "dataset", "variant", "seed", "config_hash", "data_sha256", "remote_dir", "code_hash"
    ]].to_dict("index")
    rows: list[dict[str, Any]] = []
    rejected: list[dict[str, str]] = []
    seen: set[str] = set()
    for metadata_path in sorted(index_root.rglob("secondary_metadata.json")):
        run_dir = metadata_path.parent
        run_id = run_dir.name
        try:
            checksums = verify_manifest(run_dir)
            metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
            metrics = json.loads((run_dir / "metrics.json").read_text(encoding="utf-8"))
            if metadata.get("run_id") != run_id or run_id not in lookup or run_id in seen:
                raise ValueError(f"Unknown, duplicate, or path-mismatched run_id: {run_id}")
            expected = lookup[run_id]
            if (metadata.get("dataset"), metadata.get("variant"), int(metadata.get("model_seed", -1))) != (
                expected["dataset"], expected["variant"], int(expected["seed"])
            ):
                raise ValueError(f"Leiden metadata disagrees with frozen run master: {run_id}")
            if metadata.get("version") != SECONDARY_VERSION:
                raise ValueError(f"Unexpected Leiden version: {run_id}")
            if metadata.get("source_remote_dir") != expected["remote_dir"] or metadata.get("source_data_sha256") != expected["data_sha256"]:
                raise ValueError(f"Leiden source provenance changed: {run_id}")
            if set(metadata.get("source_artifacts", {})) != SOURCE_ARTIFACTS:
                raise ValueError(f"Leiden source artifact contract changed: {run_id}")
            if metadata.get("runner_sha256") != secondary_freeze["runner_sha256"]:
                raise ValueError(f"Unapproved Leiden runner: {run_id}")
            datetime.fromisoformat(str(metadata["completed_utc"]).replace("Z", "+00:00"))
            protocol = metadata.get("protocol", {})
            if protocol != {
                "n_neighbors": 15, "resolution": 1.0, "random_state": int(expected["seed"]),
                "oracle_sweep": False, "known_k": False,
            }:
                raise ValueError(f"Leiden protocol changed: {run_id}")
            values, n_pred = verify_arrays_and_metrics(
                run_dir, metrics, int(datasets[expected["dataset"]]["retained_cells"]),
                int(datasets[expected["dataset"]]["expected_clusters"]), run_id,
            )
            seen.add(run_id)
            rows.append({
                "run_id": run_id, "dataset": expected["dataset"], "variant": expected["variant"],
                "model_seed": int(expected["seed"]), "config_hash": expected["config_hash"],
                "data_sha256": expected["data_sha256"], "primary_code_hash": expected["code_hash"],
                "primary_freeze_hash": primary_freeze["freeze_hash"],
                "secondary_freeze_hash": secondary_freeze["freeze_hash"],
                "runner_sha256": metadata["runner_sha256"], "n_neighbors": 15, "resolution": 1.0,
                "uses_known_k": False, "oracle_sweep": False, "silhouette_status": "not_computed",
                "n_pred_clusters": n_pred, "completed_utc": metadata["completed_utc"],
                "remote_dir": f"{REMOTE_RESULT_ROOT}/downstream/{SECONDARY_VERSION}/{PROTOCOL_VERSION}/{run_id}",
                "payload_manifest_sha256": sha256_file(run_dir / "SHA256SUMS"), **values,
            })
        except Exception as exc:
            rejected.append({"run_id": run_id, "local_index_dir": str(run_dir), "reason": str(exc)})
    missing = sorted(set(lookup) - seen)
    if rejected or missing:
        detail = f"accepted={len(seen)}/108 rejected={len(rejected)} missing={len(missing)}"
        raise RuntimeError(f"Leiden evidence incomplete or invalid: {detail}")
    run_level = pd.DataFrame(rows).sort_values(["dataset", "variant", "model_seed"]).reset_index(drop=True)
    return run_level, pd.DataFrame(rejected, columns=["run_id", "local_index_dir", "reason"]), {
        "primary_freeze": primary_freeze, "secondary_freeze": secondary_freeze,
        "dataset_manifest_sha256": sha256_file(dataset_manifest),
    }


def summarize(run_level: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    records = []
    for (dataset, variant), group in run_level.groupby(["dataset", "variant"], sort=True):
        if set(group["model_seed"]) != set(SEEDS) or len(group) != len(SEEDS):
            raise RuntimeError(f"Seed grid incomplete for {dataset}/{variant}")
        row: dict[str, Any] = {"dataset": dataset, "variant": variant, "n_seeds": len(group)}
        for metric in METRICS:
            values = group[metric].to_numpy(float)
            row.update({
                f"{metric}_mean": float(values.mean()), f"{metric}_sd": float(values.std(ddof=1)),
                f"{metric}_median": float(np.median(values)), f"{metric}_min": float(values.min()),
                f"{metric}_max": float(values.max()),
            })
        records.append(row)
    dataset_level = pd.DataFrame(records)
    overall_records = []
    for variant, group in dataset_level.groupby("variant", sort=True):
        if len(group) != len(DATASETS):
            raise RuntimeError(f"Dataset grid incomplete for {variant}")
        row = {"variant": variant, "n_datasets": len(group)}
        for metric in METRICS:
            values = group[f"{metric}_mean"].to_numpy(float)
            row.update({
                f"{metric}_mean": float(values.mean()), f"{metric}_sd": float(values.std(ddof=1)),
                f"{metric}_median": float(np.median(values)),
            })
        overall_records.append(row)
    return dataset_level, pd.DataFrame(overall_records)


def hierarchical_ci(
    run_level: pd.DataFrame, left: str, right: str, metric: str, rng: np.random.Generator
) -> tuple[float, float]:
    paired = run_level.pivot(index=["dataset", "model_seed"], columns="variant", values=metric)[[left, right]].dropna()
    by_dataset = {dataset: group[left].to_numpy() - group[right].to_numpy() for dataset, group in paired.groupby(level=0)}
    datasets = np.asarray(sorted(by_dataset), dtype=object)
    draws = np.empty(ITERATIONS, dtype=float)
    for index in range(ITERATIONS):
        sampled_datasets = rng.choice(datasets, size=len(datasets), replace=True)
        dataset_means = [rng.choice(by_dataset[name], size=len(by_dataset[name]), replace=True).mean() for name in sampled_datasets]
        draws[index] = np.mean(dataset_means)
    return tuple(np.percentile(draws, [2.5, 97.5]).astype(float))


def contrast_table(run_level: pd.DataFrame, dataset_level: pd.DataFrame) -> pd.DataFrame:
    rng = np.random.default_rng(RNG_SEED)
    rows: list[dict[str, Any]] = []
    for metric in CONFIRMATORY_METRICS:
        pivot = dataset_level.pivot(index="dataset", columns="variant", values=f"{metric}_mean")
        metric_rows = []
        for left, right, name in CONTRASTS:
            delta = (pivot[left] - pivot[right]).to_numpy(float)
            if len(delta) != len(DATASETS):
                raise RuntimeError(f"Contrast {name}/{metric} does not contain six datasets")
            lo, hi = hierarchical_ci(run_level, left, right, metric, rng)
            metric_rows.append({
                "metric": metric, "contrast": name, "n_datasets": len(delta),
                "mean_delta": float(delta.mean()), "median_delta": float(np.median(delta)),
                "ci95_low": lo, "ci95_high": hi,
                "wins": int((delta > TIE_TOLERANCE).sum()),
                "ties": int((np.abs(delta) <= TIE_TOLERANCE).sum()),
                "losses": int((delta < -TIE_TOLERANCE).sum()),
                "permutation_p": sign_flip_pvalue(delta, ITERATIONS, rng),
                "bootstrap": "paired hierarchical: datasets then paired model seeds",
            })
        for row, adjusted in zip(metric_rows, holm_adjust([row["permutation_p"] for row in metric_rows])):
            row["holm_p"] = adjusted
            rows.append(row)
    return pd.DataFrame(rows)


def write_frame(frame: pd.DataFrame, stem: Path) -> None:
    frame.to_csv(stem.with_suffix(".csv"), index=False)
    stem.with_suffix(".json").write_text(frame.to_json(orient="records", indent=2), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description="Aggregate and independently recompute fixed Leiden evidence")
    parser.add_argument("--index-root", type=Path, default=PAPER_ROOT / f".staging/remote_secondary_index/{SECONDARY_VERSION}")
    parser.add_argument("--run-master", type=Path, default=PAPER_ROOT / f"experiments/{PROTOCOL_VERSION}/run_master.csv")
    parser.add_argument("--freeze", type=Path, default=PAPER_ROOT / f"experiments/{SECONDARY_VERSION}/source_freeze.json")
    parser.add_argument("--dataset-manifest", type=Path, default=PAPER_ROOT / "manifests/dataset_upload/dataset_manifest.json")
    parser.add_argument(
        "--output-dir", type=Path,
        default=PAPER_ROOT / f"experiments/{SECONDARY_VERSION}/aggregate",
    )
    parser.add_argument("--no-sync", action="store_true")
    args = parser.parse_args()
    if not args.no_sync:
        sync_index(args.index_root)
    run_level, rejected, provenance = collect(args.index_root, args.run_master, args.freeze, args.dataset_manifest)
    dataset_level, overall = summarize(run_level)
    contrasts = contrast_table(run_level, dataset_level)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    write_frame(run_level, args.output_dir / "run_metrics")
    write_frame(dataset_level, args.output_dir / "dataset_variant_metrics")
    write_frame(overall, args.output_dir / "variant_overall_metrics")
    write_frame(contrasts, args.output_dir / "contrasts")
    write_frame(rejected, args.output_dir / "rejected_or_superseded_runs")
    write_json(args.output_dir / "aggregation_metadata.json", {
        "generated_utc": utc_now(), "secondary_version": SECONDARY_VERSION,
        "accepted_runs": len(run_level), "rejected_runs": len(rejected), "expected_runs": 108,
        "dataset_is_independent_unit": True, "no_best_seed_selection": True,
        "confirmatory_metrics": list(CONFIRMATORY_METRICS),
        "descriptive_metrics": sorted(set(METRICS) - set(CONFIRMATORY_METRICS)),
        "bootstrap_iterations": ITERATIONS, "permutation_iterations": ITERATIONS,
        "rng_seed": RNG_SEED, "tie_tolerance": TIE_TOLERANCE,
        "silhouette": "not computed by the frozen secondary runner", **provenance,
    })
    write_checksum_manifest(args.output_dir)
    print(f"leiden_runs={len(run_level)} datasets={run_level['dataset'].nunique()} contrasts={len(contrasts)}")


if __name__ == "__main__":
    main()
