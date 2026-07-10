from __future__ import annotations

import argparse
import json
import shutil
from pathlib import Path, PurePosixPath

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.metrics import recall_score, roc_auc_score

from .config import PAPER_ROOT, REMOTE_RESULT_ROOT
from .io_utils import require_disk_space, sha256_file, write_json
from .remote_store import RemoteStore
from .run_stress import STRESS_VERSION


def safe_spearman(x: np.ndarray, y: np.ndarray) -> float:
    x = np.asarray(x, dtype=float); y = np.asarray(y, dtype=float)
    valid = np.isfinite(x) & np.isfinite(y)
    if valid.sum() < 3 or np.unique(x[valid]).size < 2 or np.unique(y[valid]).size < 2:
        return float("nan")
    return float(spearmanr(x[valid], y[valid]).statistic)


def safe_auroc(y: np.ndarray, score: np.ndarray) -> float:
    y = np.asarray(y, dtype=bool).ravel(); score = np.asarray(score, dtype=float).ravel()
    valid = np.isfinite(score)
    if valid.sum() < 2 or np.unique(y[valid]).size < 2:
        return float("nan")
    return float(roc_auc_score(y[valid], score[valid]))


def parse_manifest(path: Path) -> dict[str, str]:
    result = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.strip():
            digest, name = line.split("  ", 1); result[name.lstrip("*")] = digest
    return result


def download_verified(store: RemoteStore, remote_dir: str, name: str, local_dir: Path, checksums: dict[str, str]) -> Path:
    if name not in checksums:
        raise FileNotFoundError(f"{name} is absent from {remote_dir}/SHA256SUMS")
    target = local_dir / name
    store.download_file(f"{remote_dir}/{name}", target)
    if sha256_file(target) != checksums[name]:
        raise ValueError(f"Checksum mismatch: {remote_dir}/{name}")
    return target


def analyze_one(metadata: dict, bundle_path: Path) -> tuple[dict, list[dict]]:
    with np.load(bundle_path) as z:
        labels = z["labels"].astype(int)
        mapped = z["mapped"].astype(int)
        indices = z["indices"].astype(int)
        base = z["base_probs"].astype(float)
        rel = z["edge_reliability"].astype(float)
        gate = z["node_gate"].astype(float)
    affinity = base * rel
    affinity /= np.clip(affinity.sum(axis=1, keepdims=True), 1e-12, None)
    same = labels[indices] == labels[:, None]
    purity = np.sum(affinity * same, axis=1)
    class_values, counts = np.unique(labels, return_counts=True)
    frequency = dict(zip(class_values.tolist(), counts.tolist()))
    cell_frequency = np.asarray([frequency[value] for value in labels], dtype=float)
    recall = recall_score(labels, mapped, labels=class_values, average=None, zero_division=0)
    class_rows = []
    for value, count, rec in zip(class_values, counts, recall):
        mask = labels == value
        class_rows.append(
            {
                "run_id": metadata["run_id"], "class_id": int(value), "n_cells": int(count),
                "recall": float(rec), "mean_gate": float(gate[mask].mean()),
                "mean_neighbor_purity": float(purity[mask].mean()),
            }
        )
    metrics = metadata["validation"]["metrics"]
    protocol = metadata["protocol"]
    row = {
        "run_id": metadata["run_id"], "dataset": metadata["dataset"],
        "variant": metadata["variant"], "seed": metadata["seed"],
        "contamination": float(protocol["stress_bad_edge_ratio"]),
        "estimator": protocol["neighbor_estimator"], "config_hash": metadata["config_hash"],
        "ari": metrics.get("ari"), "nmi": metrics.get("nmi"), "acc": metrics.get("acc"),
        "f1_macro": metrics.get("f1_macro"), "runtime_seconds": metadata.get("runtime_seconds"),
        "same_edge_fraction": float(same.mean()), "weighted_same_edge_purity": float(purity.mean()),
        "affinity_same_edge_auroc": safe_auroc(same, affinity),
        "affinity_same_edge_spearman": safe_spearman(affinity.ravel(), same.astype(float).ravel()),
        "gate_purity_spearman": safe_spearman(gate, purity),
        "gate_class_frequency_spearman": safe_spearman(gate, cell_frequency),
        "class_frequency_recall_spearman": safe_spearman(counts.astype(float), recall),
        "source_freeze_hash": metadata.get("source_freeze_hash"),
    }
    return row, class_rows


def main() -> None:
    parser = argparse.ArgumentParser(description="Stream and aggregate frozen scVICAR stress diagnostics")
    parser.add_argument("--staging", type=Path, default=PAPER_ROOT / ".staging/stress_analysis")
    parser.add_argument("--output-dir", type=Path, default=PAPER_ROOT / "experiments/stress_v1")
    parser.add_argument("--require-complete", action="store_true")
    args = parser.parse_args()
    require_disk_space(args.staging, 5.0)
    store = RemoteStore()
    root = f"{REMOTE_RESULT_ROOT}/runs/{STRESS_VERSION}"
    found = store.run(f"find {root} -type f -name COMPLETED -print", check=False)
    remote_dirs = sorted({str(PurePosixPath(line.strip()).parent) for line in found.stdout.splitlines() if line.strip()})
    expected = 126
    if args.require_complete and len(remote_dirs) != expected:
        raise RuntimeError(f"Stress matrix incomplete: {len(remote_dirs)}/{expected}")
    rows: list[dict] = []; class_rows: list[dict] = []
    for index, remote_dir in enumerate(remote_dirs):
        require_disk_space(args.staging, 5.0)
        local = args.staging / f"run_{index:04d}"; local.mkdir(parents=True, exist_ok=True)
        manifest_path = local / "SHA256SUMS"; store.download_file(f"{remote_dir}/SHA256SUMS", manifest_path)
        checksums = parse_manifest(manifest_path)
        metadata_path = download_verified(store, remote_dir, "run_metadata.json", local, checksums)
        bundle_path = download_verified(store, remote_dir, "topology_diagnostics.npz", local, checksums)
        metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
        row, classes = analyze_one(metadata, bundle_path); rows.append(row); class_rows.extend(classes)
        shutil.rmtree(local)
    frame = pd.DataFrame(rows); classes = pd.DataFrame(class_rows)
    if not frame.empty:
        baseline = frame[(frame["contamination"] == 0) & (frame["estimator"] == "current")][
            ["dataset", "variant", "seed", "ari"]
        ].rename(columns={"ari": "ari_clean"})
        frame = frame.merge(baseline, on=["dataset", "variant", "seed"], how="left")
        frame["ari_degradation"] = frame["ari"] - frame["ari_clean"]
    args.output_dir.mkdir(parents=True, exist_ok=True)
    frame.to_csv(args.output_dir / "stress_runs.csv", index=False)
    classes.to_csv(args.output_dir / "stress_class_diagnostics.csv", index=False)
    if not frame.empty:
        summary = frame.groupby(["dataset", "variant", "contamination", "estimator"], as_index=False).agg(
            ari_mean=("ari", "mean"), ari_std=("ari", "std"),
            degradation_mean=("ari_degradation", "mean"),
            purity_mean=("weighted_same_edge_purity", "mean"),
            affinity_auroc_mean=("affinity_same_edge_auroc", "mean"),
            gate_purity_spearman_mean=("gate_purity_spearman", "mean"),
        )
        summary.to_csv(args.output_dir / "stress_summary.csv", index=False)
    write_json(args.output_dir / "stress_traceability.json", {"remote_runs": len(remote_dirs), "expected": expected, "remote_dirs": remote_dirs})


if __name__ == "__main__":
    main()
