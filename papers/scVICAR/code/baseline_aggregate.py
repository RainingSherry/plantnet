from __future__ import annotations

import argparse
import json
import subprocess
from pathlib import Path

import pandas as pd

from .config import DATASETS, PAPER_ROOT, REMOTE_RESULT_ROOT, SEEDS, sha256_payload
from .downstream_orchestrate import _parse_sha256sums
from .io_utils import sha256_file
from .remote_store import RemoteStore
from .run_baseline import BASELINES, BASELINE_VERSION


SYNC_NAMES = ("run_metadata.json", "metrics.json", "COMPLETED", "SHA256SUMS")


def verify_selected_files(run_dir: Path) -> None:
    manifest = run_dir / "SHA256SUMS"
    if not manifest.is_file():
        raise FileNotFoundError(manifest)
    checksums = _parse_sha256sums(manifest)
    for name in ("run_metadata.json", "metrics.json", "COMPLETED"):
        path = run_dir / name
        if name not in checksums or not path.is_file():
            raise FileNotFoundError(f"Missing checksummed baseline artifact: {path}")
        if sha256_file(path) != checksums[name]:
            raise ValueError(f"Baseline checksum mismatch: {path}")


def sync_metadata(local_root: Path) -> None:
    store = RemoteStore()
    local_root.mkdir(parents=True, exist_ok=True)
    command = ["rsync", "-a", "--prune-empty-dirs", "-e", store.rsync_shell, "--include=*/"]
    command.extend(f"--include={name}" for name in SYNC_NAMES)
    command.extend([
        "--exclude=*",
        f"{store.target}:{REMOTE_RESULT_ROOT}/runs/{BASELINE_VERSION}/",
        f"{local_root}/",
    ])
    subprocess.run(command, check=True)


def collect(
    local_root: Path,
    freeze: dict,
    primary_freeze: dict,
    scdeepcluster_repair: dict,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    accepted: list[dict] = []
    rejected: list[dict] = []
    for meta_path in sorted(local_root.rglob("run_metadata.json")):
        run_dir = meta_path.parent
        try:
            verify_selected_files(run_dir)
            metadata = json.loads(meta_path.read_text(encoding="utf-8"))
        except Exception as exc:
            rejected.append({"run_id": run_dir.name, "reason": f"artifact_validation: {exc}"})
            continue
        record = {
            "run_id": metadata.get("run_id"), "dataset": metadata.get("dataset"),
            "method": metadata.get("method"), "seed": metadata.get("seed"),
            "config_hash": metadata.get("config_hash"), "runtime_seconds": metadata.get("runtime_seconds"),
            "remote_dir": metadata.get("remote_dir"),
            "known_k": metadata.get("protocol", {}).get("known_k"),
            "known_k_training": metadata.get("protocol", {}).get("known_k_training"),
            "per_cell_labels_in_optimization": metadata.get("protocol", {}).get("per_cell_labels_in_optimization"),
        }
        reason = None
        method = metadata.get("method")
        protocol = metadata.get("protocol", {})
        repaired_scdeepcluster = method == "scdeepcluster"
        expected_baseline_freeze = (
            scdeepcluster_repair["freeze_hash"] if repaired_scdeepcluster else freeze["freeze_hash"]
        )
        identity_keys = (
            "python", "python_sha256", "distribution_lock", "distribution_lock_hash",
            "script", "gpu", "runner_args", "code_sha256",
        )
        protocol_identity = {key: protocol.get(key) for key in identity_keys}
        if not (run_dir / "COMPLETED").is_file():
            reason = "missing_completed"
        elif metadata.get("baseline_freeze_hash") != expected_baseline_freeze:
            reason = "baseline_freeze_mismatch"
        elif metadata.get("source_freeze_hash") != primary_freeze["freeze_hash"]:
            reason = "primary_freeze_mismatch"
        elif metadata.get("data", {}).get("sha256") != primary_freeze["dataset_sha256"].get(metadata.get("dataset")):
            reason = "dataset_sha256_mismatch"
        elif metadata.get("protocol", {}).get("baseline_version") != BASELINE_VERSION:
            reason = "baseline_version_mismatch"
        elif metadata.get("protocol", {}).get("method") != metadata.get("method"):
            reason = "method_protocol_mismatch"
        elif metadata.get("protocol", {}).get("seed") != metadata.get("seed"):
            reason = "seed_protocol_mismatch"
        elif repaired_scdeepcluster and sha256_payload(protocol_identity) != scdeepcluster_repair["identity_hash"]:
            reason = "scdeepcluster_repair_identity_mismatch"
        elif not repaired_scdeepcluster and protocol.get("code_sha256") != freeze["methods"].get(method, {}).get("code_sha256"):
            reason = "method_code_mismatch"
        elif not metadata.get("validation", {}).get("cell_set_complete"):
            reason = "cell_set_incomplete"
        elif not metadata.get("validation", {}).get("cell_id_order_exact"):
            reason = "cell_id_order_mismatch"
        elif not str(metadata.get("run_id", "")).endswith(f"--{metadata.get('config_hash')}"):
            reason = "run_id_config_hash_mismatch"
        elif not str(metadata.get("remote_dir", "")).endswith(f"/{metadata.get('config_hash')}"):
            reason = "remote_config_hash_mismatch"
        if reason:
            rejected.append({**record, "reason": reason})
            continue
        metrics = metadata.get("validation", {}).get("metrics", {})
        accepted.append({
            **record,
            "data_sha256": metadata.get("data", {}).get("sha256"),
            "baseline_freeze_hash": metadata.get("baseline_freeze_hash"),
            "source_freeze_hash": metadata.get("source_freeze_hash"),
            "cell_id_order_exact": metadata.get("validation", {}).get("cell_id_order_exact"),
            **{key: metrics.get(key) for key in (
                "ari", "nmi", "acc", "f1_macro", "silhouette", "n_pred_clusters"
            )},
        })
    return pd.DataFrame(accepted), pd.DataFrame(rejected)


def main() -> None:
    parser = argparse.ArgumentParser(description="Aggregate only frozen, complete external baseline runs")
    parser.add_argument("--index-root", type=Path, default=PAPER_ROOT / ".staging/remote_baseline_index")
    parser.add_argument("--output-dir", type=Path, default=PAPER_ROOT / "experiments/baselines_v1")
    parser.add_argument("--no-sync", action="store_true")
    parser.add_argument("--require-complete", action="store_true")
    args = parser.parse_args()
    if not args.no_sync:
        sync_metadata(args.index_root)
    baseline_freeze = json.loads((args.output_dir / "source_freeze.json").read_text(encoding="utf-8"))
    scdeepcluster_repair = json.loads(
        (args.output_dir / "source_freeze_scdeepcluster_repair_v2.json").read_text(encoding="utf-8")
    )
    primary_freeze = json.loads(
        (PAPER_ROOT / "experiments/protocol_v1/source_freeze.json").read_text(encoding="utf-8")
    )
    if baseline_freeze["primary_freeze_hash"] != primary_freeze["freeze_hash"]:
        raise ValueError("Baseline freeze does not bind the active primary freeze")
    repair_payload = {key: value for key, value in scdeepcluster_repair.items() if key != "freeze_hash"}
    if sha256_payload(repair_payload) != scdeepcluster_repair["freeze_hash"]:
        raise ValueError("scDeepCluster repair freeze hash is invalid")
    frame, rejected = collect(
        args.index_root, baseline_freeze, primary_freeze, scdeepcluster_repair
    )
    args.output_dir.mkdir(parents=True, exist_ok=True)
    frame.to_csv(args.output_dir / "run_master.csv", index=False)
    rejected.to_csv(args.output_dir / "superseded_or_invalid_runs.csv", index=False)
    expected = len(DATASETS) * len(BASELINES) * len(SEEDS)
    unique = len(frame[["dataset", "method", "seed"]].drop_duplicates()) if not frame.empty else 0
    duplicate_run_ids = int(frame["run_id"].duplicated().sum()) if not frame.empty else 0
    expected_keys = {
        (dataset, method, seed)
        for dataset in DATASETS for method in BASELINES for seed in SEEDS
    }
    observed_keys = set(map(tuple, frame[["dataset", "method", "seed"]].itertuples(index=False, name=None))) if not frame.empty else set()
    if args.require_complete and (
        len(frame) != expected or unique != expected or duplicate_run_ids or observed_keys != expected_keys
    ):
        raise RuntimeError(f"Baseline matrix incomplete: rows={len(frame)}, unique={unique}, expected={expected}")
    if not frame.empty:
        frame.groupby("method")[["ari", "nmi", "acc", "f1_macro", "runtime_seconds"]].agg(
            ["mean", "std", "count"]
        ).to_csv(args.output_dir / "summary_by_method.csv")
        frame.groupby(["dataset", "method"])[["ari", "nmi", "acc", "f1_macro"]].agg(
            ["mean", "std", "count"]
        ).to_csv(args.output_dir / "summary_by_dataset_method.csv")
    print(f"accepted={len(frame)} rejected={len(rejected)}")


if __name__ == "__main__":
    main()
