from __future__ import annotations

import argparse
import json
import shlex
from pathlib import Path, PurePosixPath

from .config import PAPER_ROOT, REMOTE_RESULT_ROOT
from .io_utils import utc_now, write_json
from .remote_store import RemoteStore
from .run_stress import STRESS_DATASETS, STRESS_VERSION
from .stress_orchestrate import tasks_for_dataset


def main() -> int:
    parser = argparse.ArgumentParser(description="Audit stress runs for checksums, frozen identity, and actual CUDA use")
    parser.add_argument("--require-complete", action="store_true")
    parser.add_argument(
        "--output", type=Path,
        default=PAPER_ROOT / "experiments/stress_v1/runtime_audit.json",
    )
    args = parser.parse_args()
    primary = json.loads((PAPER_ROOT / "experiments/protocol_v1/source_freeze.json").read_text(encoding="utf-8"))
    stress = json.loads((PAPER_ROOT / "experiments/stress_v1/source_freeze.json").read_text(encoding="utf-8"))
    expected_keys = {
        (dataset, variant, seed, float(contamination), estimator)
        for dataset in STRESS_DATASETS
        for variant, seed, contamination, estimator in tasks_for_dataset(dataset)
    }
    store = RemoteStore()
    root = f"{REMOTE_RESULT_ROOT}/runs/{STRESS_VERSION}"
    found = store.run(f"find {shlex.quote(root)} -type f -name COMPLETED -print", check=False)
    remote_dirs = sorted({str(PurePosixPath(line).parent) for line in found.stdout.splitlines() if line.strip()})
    accepted = []
    rejected = []
    observed_keys = set()
    for remote_dir in remote_dirs:
        try:
            verify = store.run(
                f"cd {shlex.quote(remote_dir)} && sha256sum -c SHA256SUMS",
                check=False,
            )
            if verify.returncode != 0:
                raise ValueError("remote SHA256SUMS validation failed")
            first_line = store.run(f"head -n 1 {shlex.quote(remote_dir + '/run.log')}", check=False)
            if first_line.returncode != 0 or first_line.stdout.strip() != "Using device: cuda:0":
                raise ValueError(f"actual device is not cuda:0: {first_line.stdout.strip()!r}")
            metadata_result = store.run(f"cat {shlex.quote(remote_dir + '/run_metadata.json')}", check=False)
            if metadata_result.returncode != 0:
                raise FileNotFoundError("run_metadata.json")
            metadata = json.loads(metadata_result.stdout)
            protocol = metadata["protocol"]
            key = (
                metadata["dataset"], metadata["variant"], int(metadata["seed"]),
                float(protocol["stress_bad_edge_ratio"]), protocol["neighbor_estimator"],
            )
            if key not in expected_keys or key in observed_keys:
                raise ValueError(f"unknown or duplicate stress task key: {key}")
            if metadata.get("source_freeze_hash") != primary["freeze_hash"]:
                raise ValueError("primary freeze mismatch")
            if metadata.get("stress_freeze_hash") != stress["freeze_hash"]:
                raise ValueError("stress freeze mismatch")
            if metadata.get("data", {}).get("sha256") != primary["dataset_sha256"][metadata["dataset"]]:
                raise ValueError("dataset hash mismatch")
            resources = metadata.get("resources", {})
            if resources.get("measurement") != "nvidia-smi device memory sampled every 0.5 s on the isolated physical GPU":
                raise ValueError("GPU measurement contract mismatch")
            if resources.get("device_memory_peak_delta_mib") is None or float(resources["device_memory_peak_delta_mib"]) < 0:
                raise ValueError("invalid GPU peak-memory delta")
            if metadata.get("remote_dir") != remote_dir:
                raise ValueError("metadata remote path mismatch")
            observed_keys.add(key)
            accepted.append({"run_id": metadata["run_id"], "remote_dir": remote_dir, "task_key": list(key)})
        except Exception as exc:
            rejected.append({"remote_dir": remote_dir, "reason": str(exc)})
    missing = sorted(expected_keys - observed_keys)
    passed = not rejected and (not args.require_complete or not missing) and observed_keys.issubset(expected_keys)
    payload = {
        "created_utc": utc_now(), "passed": passed, "expected_runs": len(expected_keys),
        "accepted_runs": len(accepted), "rejected_runs": len(rejected), "missing_runs": len(missing),
        "device_requirement": "checksummed run.log first line equals 'Using device: cuda:0'",
        "accepted": accepted, "rejected": rejected, "missing_task_keys": [list(key) for key in missing],
    }
    write_json(args.output, payload)
    print(json.dumps({key: payload[key] for key in ("passed", "accepted_runs", "rejected_runs", "missing_runs")}, indent=2))
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
