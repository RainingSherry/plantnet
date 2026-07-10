from __future__ import annotations

import argparse
import json
import shutil
import traceback
from pathlib import Path

import anndata as ad
import numpy as np
import scanpy as sc

from methods.DeepLearning import scMAE_family as family

from .config import DATASETS, PAPER_ROOT, PROTOCOL_VERSION, REMOTE_RESULT_ROOT
from .downstream_orchestrate import (
    FormalRun,
    download_run_inputs,
    load_dataset_manifest,
    load_formal_runs,
    require_minimum_disk,
    verify_remote_directory,
)
from .io_utils import sha256_file, utc_now, write_checksum_manifest, write_json
from .remote_store import RemoteStore


SECONDARY_VERSION = "leiden_fixed_v1"


def process(store: RemoteStore, run: FormalRun, data_sha256: str, staging: Path) -> dict:
    remote = f"{REMOTE_RESULT_ROOT}/downstream/{SECONDARY_VERSION}/{PROTOCOL_VERSION}/{run.run_id}"
    if store.exists(f"{remote}/COMPLETED"):
        verify_remote_directory(store, remote)
        return {"run_id": run.run_id, "status": "remote_complete", "remote_dir": remote}
    require_minimum_disk(staging)
    root = staging / run.run_id
    inputs = root / "inputs"
    output = root / "result"
    try:
        source = download_run_inputs(store, run, inputs, data_sha256)
        embedding = np.load(inputs / "embedding_float32.npz")["embedding"].astype(np.float32)
        labels = np.load(inputs / "clusters.npz")["labels"].astype(np.int64)
        work = ad.AnnData(X=embedding)
        sc.pp.neighbors(
            work,
            n_neighbors=min(15, max(1, work.n_obs - 1)),
            use_rep="X",
            random_state=run.seed,
        )
        sc.tl.leiden(
            work,
            resolution=1.0,
            random_state=run.seed,
            key_added="leiden_fixed_1p0",
        )
        pred = work.obs["leiden_fixed_1p0"].cat.codes.to_numpy().astype(np.int32)
        metrics, mapped = family.compute_kmeans_metrics(labels, pred)
        metrics.update(
            {
                "cluster_method": "leiden_fixed_resolution_1p0",
                "uses_known_k": False,
                "resolution": 1.0,
                "n_neighbors": 15,
            }
        )
        output.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(
            output / "leiden_clusters.npz",
            labels=labels,
            predicted=pred,
            mapped=mapped.astype(np.int32),
        )
        write_json(output / "metrics.json", metrics)
        write_json(
            output / "secondary_metadata.json",
            {
                "version": SECONDARY_VERSION,
                "run_id": run.run_id,
                "dataset": run.dataset,
                "variant": run.variant,
                "model_seed": run.seed,
                "source_remote_dir": run.remote_dir,
                "source_data_sha256": data_sha256,
                "source_artifacts": source["source_artifacts"],
                "runner_sha256": sha256_file(Path(__file__).resolve()),
                "protocol": {
                    "n_neighbors": 15,
                    "resolution": 1.0,
                    "random_state": run.seed,
                    "oracle_sweep": False,
                    "known_k": False,
                },
                "completed_utc": utc_now(),
            },
        )
        (output / "COMPLETED").write_text(utc_now() + "\n", encoding="utf-8")
        write_checksum_manifest(output)
        store.upload_directory_atomic(output, remote)
        verify_remote_directory(store, remote)
        shutil.rmtree(root)
        return {"run_id": run.run_id, "status": "complete", "remote_dir": remote}
    except Exception as exc:
        root.mkdir(parents=True, exist_ok=True)
        write_json(
            root / "FAILED.json",
            {
                "run_id": run.run_id,
                "failed_utc": utc_now(),
                "error": str(exc),
                "traceback": traceback.format_exc(),
                "staging_preserved": True,
            },
        )
        raise


def main() -> int:
    parser = argparse.ArgumentParser(description="Fixed-resolution non-oracle Leiden evaluation")
    parser.add_argument(
        "--run-master",
        type=Path,
        default=PAPER_ROOT / f"experiments/{PROTOCOL_VERSION}/run_master.json",
    )
    parser.add_argument(
        "--dataset-manifest",
        type=Path,
        default=PAPER_ROOT / "manifests/dataset_upload/dataset_manifest.json",
    )
    parser.add_argument(
        "--staging",
        type=Path,
        default=PAPER_ROOT / f".staging/secondary/{SECONDARY_VERSION}",
    )
    parser.add_argument(
        "--status",
        type=Path,
        default=PAPER_ROOT / f"experiments/{SECONDARY_VERSION}/status.json",
    )
    args = parser.parse_args()
    runs = load_formal_runs(args.run_master)
    expected = len(DATASETS) * 6 * 3
    if len(runs) != expected:
        raise RuntimeError(f"Secondary evaluation requires exactly {expected} primary run IDs")
    manifest = load_dataset_manifest(args.dataset_manifest)
    store = RemoteStore()
    rows = []
    for run in runs:
        row = process(store, run, str(manifest[run.dataset]["sha256"]), args.staging)
        rows.append(row)
        write_json(args.status, {"updated_utc": utc_now(), "results": rows})
    print(json.dumps({"completed": len(rows)}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
