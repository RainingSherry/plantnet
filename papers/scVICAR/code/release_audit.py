from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from .config import (
    DATASETS, PAPER_ROOT, PROTOCOL_VERSION, REMOTE_DATA_ROOT, REMOTE_RESULT_ROOT, SEEDS, VARIANTS,
)
from .io_utils import sha256_file, utc_now, write_json
from .run_baseline import BASELINES
from .run_scVICAR import code_digests
from .remote_store import RemoteStore


def main() -> int:
    parser = argparse.ArgumentParser(description="Fail-closed scVICAR manuscript release audit")
    parser.add_argument("--remote", action="store_true")
    parser.add_argument("--output", type=Path, default=PAPER_ROOT / "manifests/release_audit.json")
    args = parser.parse_args()
    checks: list[dict] = []

    def record(name: str, passed: bool, detail: str) -> None:
        checks.append({"name": name, "passed": bool(passed), "detail": detail})

    primary_freeze_path = PAPER_ROOT / "experiments/protocol_v1/source_freeze.json"
    primary_freeze = json.loads(primary_freeze_path.read_text(encoding="utf-8"))
    record(
        "primary_source_freeze",
        primary_freeze["code_sha256"] == code_digests(),
        primary_freeze["freeze_hash"],
    )
    baseline_freeze = json.loads(
        (PAPER_ROOT / "experiments/baselines_v1/source_freeze.json").read_text(encoding="utf-8")
    )
    record(
        "baseline_freeze_parent",
        baseline_freeze["primary_freeze_hash"] == primary_freeze["freeze_hash"],
        baseline_freeze["freeze_hash"],
    )
    secondary_freeze_path = PAPER_ROOT / "experiments/leiden_fixed_v1/source_freeze.json"
    secondary_freeze = json.loads(secondary_freeze_path.read_text(encoding="utf-8"))
    record(
        "secondary_freeze_parent",
        secondary_freeze["parent_primary_freeze_hash"] == primary_freeze["freeze_hash"],
        secondary_freeze["freeze_hash"],
    )
    sensitivity_freeze_path = PAPER_ROOT / "experiments/sensitivity_full_labels_v1/source_freeze.json"
    if sensitivity_freeze_path.is_file():
        sensitivity_freeze = json.loads(sensitivity_freeze_path.read_text(encoding="utf-8"))
        record(
            "full_label_sensitivity_freeze_parent",
            sensitivity_freeze["parent_primary_freeze_hash"] == primary_freeze["freeze_hash"],
            sensitivity_freeze["freeze_hash"],
        )
    else:
        record("full_label_sensitivity_freeze_parent", False, f"missing {sensitivity_freeze_path}")
    manifest_path = PAPER_ROOT / "manifests/dataset_upload/dataset_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    record(
        "dataset_manifest",
        len(manifest) == len(DATASETS)
        and all(row["sha256"] == primary_freeze["dataset_sha256"][row["name"]] for row in manifest),
        f"{len(manifest)} canonical datasets; sha256={sha256_file(manifest_path)}",
    )

    matrix_specs = [
        (
            "primary_matrix",
            PAPER_ROOT / "experiments/protocol_v1/run_master.csv",
            len(DATASETS) * len(VARIANTS) * len(SEEDS),
            ["dataset", "variant", "seed"],
        ),
        (
            "baseline_matrix",
            PAPER_ROOT / "experiments/baselines_v1/run_master.csv",
            len(DATASETS) * len(BASELINES) * len(SEEDS),
            ["dataset", "method", "seed"],
        ),
        (
            "stress_matrix",
            PAPER_ROOT / "experiments/stress_v1/stress_runs.csv",
            126,
            ["dataset", "variant", "seed", "contamination", "estimator"],
        ),
    ]
    for name, path, expected, keys in matrix_specs:
        if not path.is_file():
            record(name, False, f"missing {path}")
            continue
        frame = pd.read_csv(path)
        unique = len(frame[keys].drop_duplicates())
        record(name, len(frame) == expected and unique == expected, f"rows={len(frame)}, unique={unique}, expected={expected}")

    stress_runtime_audit = PAPER_ROOT / "experiments/stress_v1/runtime_audit.json"
    if stress_runtime_audit.is_file():
        audit = json.loads(stress_runtime_audit.read_text(encoding="utf-8"))
        record(
            "stress_cuda_runtime_audit",
            bool(audit.get("passed")) and int(audit.get("accepted_runs", 0)) == 126,
            f"accepted={audit.get('accepted_runs')}, rejected={audit.get('rejected_runs')}",
        )
    else:
        record("stress_cuda_runtime_audit", False, f"missing {stress_runtime_audit}")

    downstream = PAPER_ROOT / "experiments/downstream_v1/dataset_variant_metrics.csv"
    if downstream.is_file():
        frame = pd.read_csv(downstream)
        expected = len(DATASETS) * len(VARIANTS) * 3
        record("downstream_matrix", len(frame) == expected, f"rows={len(frame)}, expected={expected}")
    else:
        record("downstream_matrix", False, f"missing {downstream}")

    leiden = PAPER_ROOT / "experiments/leiden_fixed_v1/aggregate/run_metrics.csv"
    if leiden.is_file():
        frame = pd.read_csv(leiden)
        expected = len(DATASETS) * len(VARIANTS) * len(SEEDS)
        unique = len(frame[["dataset", "variant", "model_seed"]].drop_duplicates())
        record(
            "fixed_leiden_matrix",
            len(frame) == expected and unique == expected,
            f"rows={len(frame)}, unique={unique}, expected={expected}",
        )
    else:
        record("fixed_leiden_matrix", False, f"missing {leiden}")

    sensitivity = PAPER_ROOT / "experiments/sensitivity_full_labels_v1/aggregate/run_metrics.csv"
    if sensitivity.is_file():
        frame = pd.read_csv(sensitivity)
        expected = len(DATASETS) * 3 * len(SEEDS)
        unique = len(frame[["dataset", "variant", "model_seed"]].drop_duplicates())
        record(
            "full_label_sensitivity_matrix",
            len(frame) == expected and unique == expected,
            f"rows={len(frame)}, unique={unique}, expected={expected}",
        )
    else:
        record("full_label_sensitivity_matrix", False, f"missing {sensitivity}")

    formats = ("svg", "pdf", "png", "tiff")
    missing_figures = [
        f"fig{number}.{suffix}"
        for number in range(1, 7)
        for suffix in formats
        if not any((PAPER_ROOT / "figures/final").glob(f"fig{number}_*.{suffix}"))
    ]
    record("figures_1_to_6", not missing_figures, "missing=" + ",".join(missing_figures))

    panel = json.loads((PAPER_ROOT / "configs/human_pancreas_marker_panel.json").read_text(encoding="utf-8"))
    record("pancreas_marker_panel", panel.get("verification_status") == "verified", panel.get("verification_status", "missing"))

    license_audit = PAPER_ROOT / "planning/data_license_audit.md"
    license_text = license_audit.read_text(encoding="utf-8") if license_audit.is_file() else ""
    record(
        "data_license_audit",
        "Status: **closed**" in license_text and "| Open |" not in license_text,
        "closed" if "Status: **closed**" in license_text else "open submission blocker",
    )

    if args.remote:
        store = RemoteStore()
        data_ok = all(
            store.exists(f"{REMOTE_DATA_ROOT}/datasets/confirmatory_v1/{name}.h5ad")
            for name in DATASETS
        )
        record("remote_datasets_present", data_ok, f"root={REMOTE_DATA_ROOT}")
        counts = {}
        for version, expected in (("protocol_v1", 108), ("baselines_v1", 108), ("stress_v1", 126)):
            result = store.run(
                f"find {REMOTE_RESULT_ROOT}/runs/{version} -type f -name COMPLETED -print",
                check=False,
            )
            count = len([line for line in result.stdout.splitlines() if line.strip()])
            counts[version] = count
            record(f"remote_{version}", count >= expected, f"completed_markers={count}, expected_at_least={expected}")
        result = store.run(
            f"find {REMOTE_RESULT_ROOT}/downstream/leiden_fixed_v1/{PROTOCOL_VERSION} "
            "-type f -name COMPLETED -print",
            check=False,
        )
        count = len([line for line in result.stdout.splitlines() if line.strip()])
        record("remote_fixed_leiden_v1", count >= 108, f"completed_markers={count}, expected_at_least=108")

    payload = {
        "created_utc": utc_now(),
        "passed": all(check["passed"] for check in checks),
        "checks": checks,
    }
    write_json(args.output, payload)
    print(json.dumps(payload, indent=2))
    return 0 if payload["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
