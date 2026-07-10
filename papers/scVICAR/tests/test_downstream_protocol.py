from __future__ import annotations

import json
import shutil
import subprocess
from pathlib import Path
from typing import Sequence

import numpy as np
import pytest

from papers.scVICAR.code import downstream_orchestrate as scheduler
from papers.scVICAR.code import downstream
from papers.scVICAR.code.config import SEEDS, SPLIT_SEEDS
from papers.scVICAR.code.io_utils import verify_checksum_manifest, write_checksum_manifest, write_json


class FakeStore:
    def __init__(self, files: dict[str, Path], upload_root: Path) -> None:
        self.files = files
        self.upload_root = upload_root
        self.downloads: list[str] = []
        self.uploads: list[str] = []
        self.completed: set[str] = set()

    def exists(self, remote_path: str) -> bool:
        return remote_path in self.completed

    def download_file(self, remote_path: str, local_path: Path) -> None:
        self.downloads.append(remote_path)
        local_path.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(self.files[remote_path], local_path)

    def upload_directory_atomic(self, local_dir: Path, remote_dir: str) -> None:
        verify_checksum_manifest(local_dir)
        self.uploads.append(remote_dir)
        target = self.upload_root / Path(remote_dir).name
        shutil.copytree(local_dir, target)
        self.completed.add(f"{remote_dir}/COMPLETED")

    def run(self, command: str, check: bool = True) -> subprocess.CompletedProcess[str]:
        return subprocess.CompletedProcess(["fake-ssh", command], 0, stdout="verified\n", stderr="")


def formal_run(seed: int = 42) -> scheduler.FormalRun:
    return scheduler.FormalRun(
        run_id=f"Mouse_Pancreas_1--nomix--seed{seed}--abcdef0123456789",
        dataset="Mouse_Pancreas_1",
        variant="nomix",
        seed=seed,
        remote_dir=f"/remote/runs/seed_{seed}/abcdef0123456789",
        data_sha256="dataset-sha",
    )


def make_remote_run(tmp_path: Path, run: scheduler.FormalRun) -> dict[str, Path]:
    source = tmp_path / "remote_source"
    source.mkdir()
    np.savez_compressed(source / "embedding_float32.npz", embedding=np.arange(12, dtype=np.float32).reshape(4, 3))
    np.savez_compressed(source / "clusters.npz", predicted=np.asarray([0, 0, 1, 1]))
    write_json(
        source / "run_metadata.json",
        {
            "run_id": run.run_id,
            "dataset": run.dataset,
            "variant": run.variant,
            "seed": run.seed,
            "protocol": {"execution_mode": "formal"},
            "data": {"sha256": "dataset-sha"},
        },
    )
    write_checksum_manifest(source)
    return {f"{run.remote_dir}/{path.name}": path for path in source.iterdir()}


def fake_downstream(command: Sequence[str], log_path: Path) -> None:
    command = list(command)
    output_dir = Path(command[command.index("--output-dir") + 1])
    output_dir.mkdir(parents=True, exist_ok=True)
    write_json(
        output_dir / "downstream_summary.json",
        {
            "tasks": list(scheduler.DOWNSTREAM_TASKS),
            "marker_splits": [{"split_seed": seed} for seed in SPLIT_SEEDS],
            "linear_probe": [
                {"label_fraction": fraction, "split_seed": seed}
                for fraction in scheduler.LABEL_FRACTIONS
                for seed in SPLIT_SEEDS
            ],
        },
    )
    log_path.write_text("lightweight fake; no differential expression run\n", encoding="utf-8")


def test_load_formal_runs_keeps_all_seeds_and_ignores_metrics(tmp_path: Path) -> None:
    rows = [
        {
            "run_id": formal_run(seed).run_id,
            "dataset": "Mouse_Pancreas_1",
            "variant": "nomix",
            "seed": seed,
            "remote_dir": formal_run(seed).remote_dir,
            "data_sha256": "dataset-sha",
            "execution_mode": "formal",
            "ari": -100.0 if seed == 42 else 100.0,
        }
        for seed in SEEDS
    ]
    rows.append(
        {
            "run_id": "smoke-with-best-score",
            "dataset": "Mouse_Pancreas_1",
            "variant": "nomix",
            "seed": 42,
            "remote_dir": "/remote/smoke",
            "execution_mode": "smoke",
            "ari": 999.0,
        }
    )
    path = tmp_path / "run_master.json"
    path.write_text(json.dumps(rows), encoding="utf-8")

    runs = scheduler.load_formal_runs(path)

    assert {run.seed for run in runs} == set(SEEDS)
    assert len(runs) == len(SEEDS)
    assert all("smoke" not in run.run_id for run in runs)


def test_canonical_artifacts_are_preferred_over_compatibility_aliases() -> None:
    checksums = {
        "embedding_float32.npz": "a" * 64,
        "embedding_final.npz": "b" * 64,
        "clusters.npz": "c" * 64,
        "clusters.csv": "d" * 64,
        "run_metadata.json": "e" * 64,
        "metadata.json": "f" * 64,
    }
    assert scheduler._select_source_artifacts(checksums) == {
        "embedding_float32.npz": "embedding_float32.npz",
        "clusters.npz": "clusters.npz",
        "run_metadata.json": "run_metadata.json",
    }


def test_command_freezes_all_tasks_and_split_seeds(tmp_path: Path) -> None:
    command = scheduler.downstream_command(tmp_path / "data.h5ad", tmp_path / "run", tmp_path / "out")
    task_start = command.index("--task") + 1
    task_end = command.index("--reference-fraction")
    assert command[task_start:task_end] == list(scheduler.DOWNSTREAM_TASKS)
    assert command[command.index("--split-seeds") + 1] == ",".join(map(str, SPLIT_SEEDS))


def test_verified_run_is_atomically_uploaded_then_cleaned(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(scheduler, "disk_free_gib", lambda _: 100.0)
    run = formal_run()
    store = FakeStore(make_remote_run(tmp_path, run), tmp_path / "uploaded")
    stage = tmp_path / "stage"
    dataset = tmp_path / "data.h5ad"
    dataset.write_bytes(b"not opened by fake downstream")

    result = scheduler.process_run(store, run, dataset, "dataset-sha", stage, runner=fake_downstream)

    expected_remote = scheduler.downstream_remote_dir(run)
    assert result["status"] == "complete"
    assert store.uploads == [expected_remote]
    assert not (stage / run.run_id).exists()
    uploaded = store.upload_root / run.run_id
    verify_checksum_manifest(uploaded)
    metadata = json.loads((uploaded / "downstream_metadata.json").read_text(encoding="utf-8"))
    assert metadata["split_seeds"] == list(SPLIT_SEEDS)
    assert metadata["seed_selection"].startswith("none")


def test_failure_preserves_downloaded_staging(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(scheduler, "disk_free_gib", lambda _: 100.0)
    run = formal_run()
    store = FakeStore(make_remote_run(tmp_path, run), tmp_path / "uploaded")
    dataset = tmp_path / "data.h5ad"
    dataset.write_bytes(b"fixture")

    def fail_runner(command: Sequence[str], log_path: Path) -> None:
        raise RuntimeError("intentional unit-test failure")

    with pytest.raises(RuntimeError, match="intentional"):
        scheduler.process_run(store, run, dataset, "dataset-sha", tmp_path / "stage", runner=fail_runner)

    preserved = tmp_path / "stage" / run.run_id
    assert (preserved / "run_input" / "embedding_float32.npz").is_file()
    assert (preserved / "FAILED.json").is_file()
    assert store.uploads == []


def test_dataset_download_is_cached_and_sha_verified(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(scheduler, "disk_free_gib", lambda _: 100.0)
    remote_data = tmp_path / "canonical.h5ad"
    remote_data.write_bytes(b"canonical dataset bytes")
    digest = scheduler.sha256_file(remote_data)
    remote_path = "/remote/data/Mouse_Pancreas_1.h5ad"
    store = FakeStore({remote_path: remote_data}, tmp_path / "uploaded")
    manifest = {
        "Mouse_Pancreas_1": {
            "sha256": digest,
            "canonical_filename": "Mouse_Pancreas_1.h5ad",
            "remote_path": remote_path,
        }
    }

    first = scheduler.ensure_cached_dataset(store, "Mouse_Pancreas_1", manifest, tmp_path / "cache")
    second = scheduler.ensure_cached_dataset(store, "Mouse_Pancreas_1", manifest, tmp_path / "cache")

    assert first == second
    assert scheduler.sha256_file(first) == digest
    assert store.downloads == [remote_path]


def test_scheduler_stops_below_five_gib(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(scheduler, "disk_free_gib", lambda _: 4.99)
    with pytest.raises(scheduler.LowDiskSpaceError, match="stopping before scheduling"):
        scheduler.require_minimum_disk(tmp_path)


def test_missing_reference_type_counts_as_zero_marker_recovery() -> None:
    import pandas as pd

    reference = {
        "A": pd.DataFrame({"names": ["g1", "g2"]}),
        "B": pd.DataFrame({"names": ["g3", "g4"]}),
    }
    predicted = {"0": pd.DataFrame({"names": ["g1", "g2"]})}
    result = downstream.matched_marker_metrics(reference, predicted)
    assert result["recovery_at_20"] == pytest.approx(0.05)


def test_linear_probe_runs_with_frozen_protocol_on_current_sklearn() -> None:
    rng = np.random.default_rng(7)
    labels = np.repeat(np.asarray(["A", "B", "C"]), 20)
    embedding = rng.normal(size=(60, 5)).astype(np.float32)
    result = downstream.run_linear_probe(embedding, labels, fraction=0.3, split_seed=11)
    assert result["train_cells"] == 18
    assert result["test_cells"] == 42
    assert set(result["per_class_recall"]) == {"A", "B", "C"}
