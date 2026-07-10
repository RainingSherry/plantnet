from __future__ import annotations

import argparse
import json
import subprocess
from pathlib import Path

import pandas as pd

from .config import DATASETS, PAPER_ROOT, REMOTE_RESULT_ROOT, SEEDS, VARIANTS, sha256_payload
from .io_utils import write_json
from .remote_store import RemoteStore


SYNC_NAMES = [
    "run_metadata.json", "summary.json", "metrics.json", "eval_fixed.csv", "gate_summary.json",
    "edge_weight_summary.json", "embedding_geometry_summary.csv",
    "rare_cell_effect_summary.csv", "COMPLETED",
]


def sync_metadata(local_root: Path) -> None:
    store = RemoteStore()
    local_root.mkdir(parents=True, exist_ok=True)
    command = ["rsync", "-a", "--prune-empty-dirs", "-e", store.rsync_shell]
    command.extend(["--include=*/"])
    command.extend(f"--include={name}" for name in SYNC_NAMES)
    command.extend(["--exclude=*", f"{store.target}:{REMOTE_RESULT_ROOT}/runs/protocol_v1/", f"{local_root}/"])
    subprocess.run(command, check=True)


def collect(local_root: Path) -> pd.DataFrame:
    rows: list[dict] = []
    for meta_path in sorted(local_root.rglob("run_metadata.json")):
        run_dir = meta_path.parent
        if not (run_dir / "COMPLETED").is_file():
            continue
        metadata = json.loads(meta_path.read_text(encoding="utf-8"))
        validation = metadata.get("validation", {})
        metrics = validation.get("metrics", {})
        row = {
            "run_id": metadata.get("run_id"),
            "dataset": metadata.get("dataset"),
            "variant": metadata.get("variant"),
            "seed": metadata.get("seed"),
            "config_hash": metadata.get("config_hash"),
            "runtime_seconds": metadata.get("runtime_seconds"),
            "gpu_memory_peak_mib": metadata.get("resources", {}).get("device_memory_peak_mib"),
            "gpu_memory_peak_delta_mib": metadata.get("resources", {}).get("device_memory_peak_delta_mib"),
            "remote_dir": metadata.get("remote_dir"),
            "data_sha256": metadata.get("data", {}).get("sha256"),
            "execution_mode": metadata.get("protocol", {}).get("execution_mode"),
            "code_hash": sha256_payload(metadata.get("code_sha256", {})),
            **{key: metrics.get(key) for key in ["ari", "nmi", "acc", "f1_macro", "silhouette", "n_pred_clusters"]},
        }
        gate_path = run_dir / "gate_summary.json"
        if gate_path.is_file():
            gate = json.loads(gate_path.read_text(encoding="utf-8"))
            row.update({f"gate_{key}": gate.get(key) for key in ["mean_node_gate", "min_node_gate", "max_node_gate", "mean_perturb_proxy"]})
        edge_path = run_dir / "edge_weight_summary.json"
        if edge_path.is_file():
            edge = json.loads(edge_path.read_text(encoding="utf-8"))
            row.update({f"edge_{key}": edge.get(key) for key in ["edge_weight_entropy", "effective_neighbor_count", "max_edge_weight_mean"]})
        geometry_path = run_dir / "embedding_geometry_summary.csv"
        if geometry_path.is_file():
            geometry = pd.read_csv(geometry_path).iloc[0].to_dict()
            row.update({f"geometry_{key}": geometry.get(key) for key in [
                "within_class_distance", "between_class_distance", "between_within_ratio"
            ]})
        rows.append(row)
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--index-root", type=Path, default=PAPER_ROOT / ".staging" / "remote_index")
    parser.add_argument("--output-dir", type=Path, default=PAPER_ROOT / "experiments" / "protocol_v1")
    parser.add_argument("--no-sync", action="store_true")
    parser.add_argument(
        "--source-freeze",
        type=Path,
        default=PAPER_ROOT / "experiments" / "protocol_v1" / "source_freeze.json",
    )
    parser.add_argument("--require-complete", action="store_true")
    args = parser.parse_args()
    if not args.no_sync:
        sync_metadata(args.index_root)
    frame = collect(args.index_root)
    if not args.source_freeze.is_file():
        raise FileNotFoundError(f"Missing frozen source identity: {args.source_freeze}")
    freeze = json.loads(args.source_freeze.read_text(encoding="utf-8"))
    frozen_code_hash = sha256_payload(freeze["code_sha256"])
    superseded = frame[frame["code_hash"] != frozen_code_hash].copy() if not frame.empty else frame.copy()
    frame = frame[frame["code_hash"] == frozen_code_hash].copy() if not frame.empty else frame
    args.output_dir.mkdir(parents=True, exist_ok=True)
    superseded.to_csv(args.output_dir / "superseded_runs.csv", index=False)
    frame.to_csv(args.output_dir / "run_master.csv", index=False)
    write_json(args.output_dir / "run_master.json", frame.where(pd.notna(frame), None).to_dict(orient="records"))
    if not frame.empty:
        formal = frame[frame["execution_mode"] == "formal"].copy()
        if not formal.empty:
            key_count = len(formal[["dataset", "variant", "seed"]].drop_duplicates())
            expected = len(DATASETS) * len(VARIANTS) * len(SEEDS)
            if args.require_complete and (key_count != expected or len(formal) != expected):
                raise RuntimeError(
                    f"Frozen formal matrix is not exactly complete: rows={len(formal)}, "
                    f"unique_keys={key_count}, expected={expected}"
                )
            summary = formal.groupby("variant")[[
                "ari", "nmi", "acc", "f1_macro", "runtime_seconds",
                "gpu_memory_peak_mib", "gpu_memory_peak_delta_mib",
            ]].agg(["mean", "std", "count"])
            summary.to_csv(args.output_dir / "summary_by_variant.csv")
            per_dataset = formal.groupby(["dataset", "variant"])[["ari", "nmi", "acc", "f1_macro"]].agg(["mean", "std", "count"])
            per_dataset.to_csv(args.output_dir / "summary_by_dataset_variant.csv")
    print(f"collected={len(frame)}")


if __name__ == "__main__":
    main()
