from __future__ import annotations

import argparse
import csv
import json
import os
import shlex
import shutil
import subprocess
import sys
import traceback
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Iterable, Mapping, Sequence

import numpy as np
import pandas as pd

from .config import (
    DATASETS,
    PAPER_ROOT,
    PROTOCOL_VERSION,
    REMOTE_DATA_ROOT,
    REMOTE_RESULT_ROOT,
    SEEDS,
    SPLIT_SEEDS,
    VARIANTS,
)
from .io_utils import disk_free_gib, sha256_file, utc_now, verify_checksum_manifest, write_checksum_manifest, write_json
from .remote_store import RemoteStore


DOWNSTREAM_TASKS = ("marker_recovery", "marker_annotation", "linear_probe")
LABEL_FRACTIONS = (0.1, 0.3)
MINIMUM_FREE_GIB = 5.0

# The first item in every tuple is the formal run_scVICAR.py contract.  The
# aliases are download-only compatibility for results created during early
# protocol development.
RUN_ARTIFACT_CANDIDATES: dict[str, tuple[str, ...]] = {
    "embedding_float32.npz": ("embedding_float32.npz", "embedding_final.npz"),
    "clusters.npz": ("clusters.npz", "clusters.csv"),
    "run_metadata.json": ("run_metadata.json", "metadata.json"),
}


class LowDiskSpaceError(RuntimeError):
    """Raised before another run is scheduled when the 5-GiB floor is crossed."""


@dataclass(frozen=True)
class FormalRun:
    run_id: str
    dataset: str
    variant: str
    seed: int
    remote_dir: str
    data_sha256: str | None = None

    @classmethod
    def from_row(cls, row: Mapping[str, Any]) -> "FormalRun":
        missing = [key for key in ("run_id", "dataset", "variant", "seed", "remote_dir") if row.get(key) in (None, "")]
        if missing:
            raise ValueError(f"Formal run is missing fields {missing}: {dict(row)}")
        run_id = str(row["run_id"])
        if run_id in {".", ".."} or "/" in run_id or "\\" in run_id:
            raise ValueError(f"Unsafe run_id: {run_id!r}")
        dataset = str(row["dataset"])
        variant = str(row["variant"])
        seed = int(row["seed"])
        if dataset not in DATASETS:
            raise ValueError(f"Unknown formal dataset {dataset!r}")
        if variant not in VARIANTS:
            raise ValueError(f"Unknown formal variant {variant!r}")
        if seed not in SEEDS:
            raise ValueError(f"Seed {seed} is outside the preregistered set {SEEDS}")
        remote_dir = str(row["remote_dir"]).rstrip("/")
        public_root = "<SCVICAR_RESULT_ROOT>"
        if remote_dir == public_root or remote_dir.startswith(public_root + "/"):
            remote_dir = REMOTE_RESULT_ROOT + remote_dir[len(public_root):]
        return cls(
            run_id=run_id,
            dataset=dataset,
            variant=variant,
            seed=seed,
            remote_dir=remote_dir,
            data_sha256=str(row["data_sha256"]) if row.get("data_sha256") else None,
        )


def _read_rows(path: Path) -> list[dict[str, Any]]:
    if path.suffix.lower() == ".csv":
        with path.open("r", encoding="utf-8", newline="") as handle:
            return [dict(row) for row in csv.DictReader(handle)]
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, list) or not all(isinstance(row, dict) for row in payload):
        raise ValueError(f"Run master must contain a JSON list of objects: {path}")
    return payload


def load_formal_runs(path: Path) -> list[FormalRun]:
    """Load every formal run; metrics are intentionally ignored."""

    rows = _read_rows(path)
    formal_rows = [row for row in rows if str(row.get("execution_mode", "")).lower() == "formal"]
    if not formal_rows:
        raise ValueError(f"No execution_mode=formal rows found in {path}")
    runs = [FormalRun.from_row(row) for row in formal_rows]
    run_ids = [run.run_id for run in runs]
    if len(run_ids) != len(set(run_ids)):
        duplicates = sorted({item for item in run_ids if run_ids.count(item) > 1})
        raise ValueError(f"Duplicate formal run_id values: {duplicates}")
    return sorted(runs, key=lambda run: (run.dataset, run.variant, run.seed, run.run_id))


def load_dataset_manifest(path: Path) -> dict[str, dict[str, Any]]:
    rows = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(rows, list):
        raise ValueError(f"Dataset manifest must be a JSON list: {path}")
    manifest = {str(row["name"]): dict(row) for row in rows}
    missing = sorted(set(DATASETS) - set(manifest))
    if missing:
        raise ValueError(f"Dataset manifest is missing preregistered datasets: {missing}")
    return manifest


def require_minimum_disk(path: Path, minimum_gib: float = MINIMUM_FREE_GIB) -> None:
    path.mkdir(parents=True, exist_ok=True)
    free = disk_free_gib(path)
    if free < minimum_gib:
        raise LowDiskSpaceError(f"Only {free:.2f} GiB free at {path}; stopping before scheduling below {minimum_gib:.2f} GiB")


def _parse_sha256sums(path: Path) -> dict[str, str]:
    checksums: dict[str, str] = {}
    for line_number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), start=1):
        if not line.strip():
            continue
        try:
            digest, relative = line.split("  ", 1)
        except ValueError as exc:
            raise ValueError(f"Malformed SHA256SUMS line {line_number}: {line!r}") from exc
        relative = relative.lstrip("*")
        if len(digest) != 64 or any(char not in "0123456789abcdefABCDEF" for char in digest):
            raise ValueError(f"Malformed SHA-256 digest on line {line_number}")
        if Path(relative).is_absolute() or ".." in Path(relative).parts:
            raise ValueError(f"Unsafe checksum path on line {line_number}: {relative!r}")
        checksums[relative] = digest.lower()
    return checksums


def _download_verified(store: RemoteStore, remote_path: str, local_path: Path, expected_sha256: str) -> None:
    if local_path.is_file() and sha256_file(local_path) == expected_sha256:
        return
    partial = local_path.with_name(local_path.name + ".partial")
    partial.unlink(missing_ok=True)
    store.download_file(remote_path, partial)
    observed = sha256_file(partial)
    if observed != expected_sha256:
        raise ValueError(f"Downloaded checksum failed for {remote_path}: expected {expected_sha256}, got {observed}")
    local_path.parent.mkdir(parents=True, exist_ok=True)
    partial.replace(local_path)


def ensure_cached_dataset(
    store: RemoteStore,
    dataset: str,
    manifest: Mapping[str, Mapping[str, Any]],
    cache_dir: Path,
) -> Path:
    require_minimum_disk(cache_dir)
    if dataset not in manifest:
        raise KeyError(f"Dataset {dataset!r} is absent from the canonical manifest")
    record = manifest[dataset]
    expected = str(record["sha256"])
    filename = str(record.get("canonical_filename") or f"{dataset}.h5ad")
    remote_path = str(record["remote_path"])
    public_root = "<SCVICAR_DATA_ROOT>"
    if remote_path == public_root or remote_path.startswith(public_root + "/"):
        remote_path = REMOTE_DATA_ROOT + remote_path[len(public_root):]
    target = cache_dir / filename
    _download_verified(store, remote_path, target, expected)
    if sha256_file(target) != expected:
        raise ValueError(f"Canonical dataset cache failed checksum verification: {dataset}")
    return target


def _select_source_artifacts(checksums: Mapping[str, str]) -> dict[str, str]:
    selected: dict[str, str] = {}
    for canonical, candidates in RUN_ARTIFACT_CANDIDATES.items():
        source = next((name for name in candidates if name in checksums), None)
        if source is None:
            raise FileNotFoundError(f"Remote SHA256SUMS has none of {candidates}")
        selected[canonical] = source
    return selected


def _normalise_compatibility_artifact(source: Path, canonical: Path) -> None:
    if source.name == canonical.name:
        return
    if canonical.name == "embedding_float32.npz":
        with np.load(source) as archive:
            key = "embedding" if "embedding" in archive.files else archive.files[0]
            embedding = np.asarray(archive[key], dtype=np.float32)
        np.savez_compressed(canonical, embedding=embedding)
        return
    if canonical.name == "clusters.npz":
        frame = pd.read_csv(source)
        candidates = ("predicted", "cluster", "predicted_cluster", "cluster_id")
        column = next((name for name in candidates if name in frame.columns), None)
        if column is None:
            usable = [name for name in frame.columns if not str(name).lower().startswith("unnamed")]
            if len(usable) != 1:
                raise ValueError(f"Cannot identify the predicted-cluster column in {source}")
            column = usable[0]
        np.savez_compressed(canonical, predicted=frame[column].to_numpy())
        return
    if canonical.name == "run_metadata.json":
        shutil.copy2(source, canonical)
        return
    raise ValueError(f"Unsupported compatibility artifact: {source.name}")


def validate_run_inputs(run: FormalRun, run_dir: Path, expected_data_sha256: str) -> dict[str, Any]:
    with np.load(run_dir / "embedding_float32.npz") as archive:
        if "embedding" not in archive.files:
            raise KeyError("embedding_float32.npz must contain 'embedding'")
        embedding = np.asarray(archive["embedding"])
    with np.load(run_dir / "clusters.npz") as archive:
        if "predicted" not in archive.files:
            raise KeyError("clusters.npz must contain 'predicted'")
        predicted = np.asarray(archive["predicted"])
    if embedding.ndim != 2 or predicted.shape != (embedding.shape[0],):
        raise ValueError("Downloaded embedding and predicted clusters have inconsistent shapes")
    if not np.isfinite(embedding).all():
        raise ValueError("Downloaded embedding contains NaN or Inf")
    metadata = json.loads((run_dir / "run_metadata.json").read_text(encoding="utf-8"))
    expected_fields = {
        "run_id": run.run_id,
        "dataset": run.dataset,
        "variant": run.variant,
        "seed": run.seed,
    }
    for field, expected in expected_fields.items():
        if metadata.get(field) != expected:
            raise ValueError(f"Run metadata mismatch for {field}: expected {expected!r}, got {metadata.get(field)!r}")
    if metadata.get("protocol", {}).get("execution_mode") != "formal":
        raise ValueError(f"Refusing non-formal remote run {run.run_id}")
    metadata_digest = metadata.get("data", {}).get("sha256")
    if metadata_digest != expected_data_sha256:
        raise ValueError(f"Dataset SHA mismatch in run metadata for {run.run_id}")
    return {"n_cells": int(embedding.shape[0]), "embedding_dim": int(embedding.shape[1])}


def download_run_inputs(
    store: RemoteStore,
    run: FormalRun,
    target_dir: Path,
    expected_data_sha256: str,
) -> dict[str, Any]:
    require_minimum_disk(target_dir.parent)
    target_dir.mkdir(parents=True, exist_ok=True)
    source_manifest = target_dir / "SOURCE_SHA256SUMS"
    store.download_file(f"{run.remote_dir}/SHA256SUMS", source_manifest)
    checksums = _parse_sha256sums(source_manifest)
    selected = _select_source_artifacts(checksums)
    source_digests: dict[str, str] = {}
    for canonical_name, source_name in selected.items():
        canonical_path = target_dir / canonical_name
        source_path = canonical_path if source_name == canonical_name else target_dir / f"source__{source_name}"
        _download_verified(store, f"{run.remote_dir}/{source_name}", source_path, checksums[source_name])
        source_digests[source_name] = checksums[source_name]
        _normalise_compatibility_artifact(source_path, canonical_path)
    validation = validate_run_inputs(run, target_dir, expected_data_sha256)
    return {"source_artifacts": selected, "source_sha256": source_digests, "validation": validation}


def downstream_command(data_path: Path, run_dir: Path, output_dir: Path) -> list[str]:
    return [
        sys.executable,
        "-m",
        "papers.scVICAR.code.downstream",
        "--data-path",
        str(data_path),
        "--run-dir",
        str(run_dir),
        "--output-dir",
        str(output_dir),
        "--task",
        *DOWNSTREAM_TASKS,
        "--reference-fraction",
        "0.5",
        "--label-fractions",
        ",".join(map(str, LABEL_FRACTIONS)),
        "--split-seeds",
        ",".join(map(str, SPLIT_SEEDS)),
    ]


def invoke_downstream(command: Sequence[str], log_path: Path) -> None:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    env.setdefault("MPLCONFIGDIR", "/tmp/scvicar-matplotlib")
    env.setdefault("NUMBA_CACHE_DIR", "/tmp/scvicar-numba")
    env.setdefault("OMP_NUM_THREADS", "4")
    env.setdefault("MKL_NUM_THREADS", "4")
    env.setdefault("OPENBLAS_NUM_THREADS", "4")
    with log_path.open("w", encoding="utf-8") as log:
        result = subprocess.run(command, env=env, stdout=log, stderr=subprocess.STDOUT, text=True, check=False)
    if result.returncode != 0:
        raise RuntimeError(f"downstream.py exited with code {result.returncode}; see {log_path}")


def validate_downstream_summary(output_dir: Path) -> dict[str, Any]:
    summary_path = output_dir / "downstream_summary.json"
    if not summary_path.is_file():
        raise FileNotFoundError(summary_path)
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    if tuple(summary.get("tasks", ())) != DOWNSTREAM_TASKS:
        raise ValueError(f"Downstream task contract changed: {summary.get('tasks')}")
    marker_seeds = [int(row["split_seed"]) for row in summary.get("marker_splits", [])]
    if marker_seeds != list(SPLIT_SEEDS):
        raise ValueError(f"Marker split seeds must be exactly {SPLIT_SEEDS}, got {marker_seeds}")
    probes = summary.get("linear_probe", [])
    observed_probe_keys = [(float(row["label_fraction"]), int(row["split_seed"])) for row in probes]
    expected_probe_keys = [(fraction, seed) for fraction in LABEL_FRACTIONS for seed in SPLIT_SEEDS]
    if observed_probe_keys != expected_probe_keys:
        raise ValueError("Linear-probe outputs do not contain the fixed fraction/seed Cartesian product")
    return summary


def downstream_remote_dir(run: FormalRun) -> str:
    return f"{REMOTE_RESULT_ROOT}/downstream/{PROTOCOL_VERSION}/{run.run_id}"


def verify_remote_directory(store: RemoteStore, remote_dir: str) -> None:
    result = store.run(
        f"cd {shlex.quote(remote_dir)} && test -f COMPLETED && sha256sum -c SHA256SUMS",
        check=False,
    )
    if result.returncode != 0:
        raise ValueError(f"Remote downstream checksum verification failed for {remote_dir}: {result.stderr}")


def _complete_local_result(output_dir: Path) -> bool:
    if not (output_dir / "COMPLETED").is_file() or not (output_dir / "SHA256SUMS").is_file():
        return False
    try:
        verify_checksum_manifest(output_dir)
        validate_downstream_summary(output_dir)
    except (FileNotFoundError, KeyError, ValueError):
        return False
    return True


def process_run(
    store: RemoteStore,
    run: FormalRun,
    dataset_path: Path,
    expected_data_sha256: str,
    staging_root: Path,
    runner: Callable[[Sequence[str], Path], None] = invoke_downstream,
) -> dict[str, Any]:
    """Process one run, preserving all local state unless remote verification succeeds."""

    require_minimum_disk(staging_root)
    stage_dir = staging_root / run.run_id
    input_dir = stage_dir / "run_input"
    output_dir = stage_dir / "result"
    remote_dir = downstream_remote_dir(run)
    if store.exists(f"{remote_dir}/COMPLETED"):
        verify_remote_directory(store, remote_dir)
        if stage_dir.exists():
            shutil.rmtree(stage_dir)
        return {"run_id": run.run_id, "status": "remote_complete", "remote_dir": remote_dir}

    try:
        source = download_run_inputs(store, run, input_dir, expected_data_sha256)
        if not _complete_local_result(output_dir):
            output_dir.mkdir(parents=True, exist_ok=True)
            command = downstream_command(dataset_path, input_dir, output_dir)
            runner(command, output_dir / "downstream.log")
            validate_downstream_summary(output_dir)
            write_json(
                output_dir / "downstream_metadata.json",
                {
                    "protocol_version": PROTOCOL_VERSION,
                    "run_id": run.run_id,
                    "dataset": run.dataset,
                    "variant": run.variant,
                    "model_seed": run.seed,
                    "source_remote_dir": run.remote_dir,
                    "target_remote_dir": remote_dir,
                    "source_data_sha256": expected_data_sha256,
                    "source_artifacts": source["source_artifacts"],
                    "source_sha256": source["source_sha256"],
                    "tasks": list(DOWNSTREAM_TASKS),
                    "split_seeds": list(SPLIT_SEEDS),
                    "label_fractions": list(LABEL_FRACTIONS),
                    "seed_selection": "none; every formal model seed is processed independently",
                    "completed_utc": utc_now(),
                    "input_validation": source["validation"],
                },
            )
            (output_dir / "COMPLETED").write_text(utc_now() + "\n", encoding="utf-8")
            write_checksum_manifest(output_dir)
            verify_checksum_manifest(output_dir)
        store.upload_directory_atomic(output_dir, remote_dir)
        verify_remote_directory(store, remote_dir)
        shutil.rmtree(stage_dir)
        return {"run_id": run.run_id, "status": "complete", "remote_dir": remote_dir}
    except Exception as exc:
        stage_dir.mkdir(parents=True, exist_ok=True)
        write_json(
            stage_dir / "FAILED.json",
            {
                "run_id": run.run_id,
                "status": "failed",
                "failed_utc": utc_now(),
                "error": str(exc),
                "traceback": traceback.format_exc(),
                "staging_preserved": True,
            },
        )
        raise


def run_scheduler(
    store: RemoteStore,
    runs: Iterable[FormalRun],
    dataset_manifest: Mapping[str, Mapping[str, Any]],
    cache_dir: Path,
    staging_root: Path,
    fail_fast: bool = True,
) -> list[dict[str, Any]]:
    results: list[dict[str, Any]] = []
    current_dataset: str | None = None
    current_dataset_path: Path | None = None
    status_path = staging_root / "scheduler_status.json"
    for run in runs:
        # Disk pressure is a scheduler-level stop, never a per-run failure to
        # skip over.  Check both filesystems before any new download or work.
        require_minimum_disk(cache_dir)
        require_minimum_disk(staging_root)
        record = dataset_manifest[run.dataset]
        expected_data_sha256 = str(record["sha256"])
        if run.data_sha256 and run.data_sha256 != expected_data_sha256:
            raise ValueError(f"Run-master dataset SHA does not match canonical manifest for {run.run_id}")
        if run.dataset != current_dataset:
            if current_dataset_path is not None:
                current_dataset_path.unlink(missing_ok=True)
            current_dataset = run.dataset
            current_dataset_path = ensure_cached_dataset(store, run.dataset, dataset_manifest, cache_dir)
        assert current_dataset_path is not None
        try:
            row = process_run(store, run, current_dataset_path, expected_data_sha256, staging_root)
        except LowDiskSpaceError:
            raise
        except Exception as exc:
            row = {"run_id": run.run_id, "status": "failed", "error": str(exc), "staging_preserved": True}
            results.append(row)
            write_json(status_path, {"updated_utc": utc_now(), "results": results})
            # A failed computation or upload retains its staging directory and
            # pauses the scheduler before any later run is started.
            raise
        results.append(row)
        write_json(status_path, {"updated_utc": utc_now(), "results": results})
    if current_dataset_path is not None:
        current_dataset_path.unlink(missing_ok=True)
    return results


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Recoverable protocol_v1 downstream scheduler for every formal run_id")
    parser.add_argument(
        "--run-master",
        type=Path,
        default=PAPER_ROOT / "experiments" / PROTOCOL_VERSION / "run_master.json",
    )
    parser.add_argument(
        "--dataset-manifest",
        type=Path,
        default=PAPER_ROOT / "manifests" / "dataset_upload" / "dataset_manifest.json",
    )
    parser.add_argument("--cache-dir", type=Path, default=PAPER_ROOT / ".staging" / "data")
    parser.add_argument(
        "--staging-root",
        type=Path,
        default=PAPER_ROOT / ".staging" / "downstream" / PROTOCOL_VERSION,
    )
    parser.add_argument("--fail-fast", action="store_true")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    runs = load_formal_runs(args.run_master)
    manifest = load_dataset_manifest(args.dataset_manifest)
    store = RemoteStore()
    store.ensure_layout()
    results = run_scheduler(store, runs, manifest, args.cache_dir, args.staging_root, args.fail_fast)
    completed = sum(row["status"] in {"complete", "remote_complete"} for row in results)
    failures = sum(row["status"] == "failed" for row in results)
    print(json.dumps({"formal_run_ids": len(runs), "completed": completed, "failed": failures}, indent=2))
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
