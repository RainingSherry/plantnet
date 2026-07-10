from __future__ import annotations

import argparse
import json
from pathlib import Path

from .config import PAPER_ROOT, SEEDS, sha256_payload
from .io_utils import sha256_file, utc_now, write_json
from .prepare_full_label_sensitivity import VERSION
from .run_full_label_sensitivity import SENSITIVITY_VARIANTS, sensitivity_code_digests
from .run_scVICAR import runtime_identity


def main() -> None:
    parser = argparse.ArgumentParser(description="Freeze full-label sensitivity data, code, runtime, and matrix")
    parser.add_argument(
        "--manifest", type=Path,
        default=PAPER_ROOT / f"manifests/{VERSION}/dataset_manifest.json",
    )
    parser.add_argument(
        "--output", type=Path,
        default=PAPER_ROOT / f"experiments/{VERSION}/source_freeze.json",
    )
    args = parser.parse_args()
    rows = json.loads(args.manifest.read_text(encoding="utf-8"))
    if len(rows) != 6 or len({row["name"] for row in rows}) != 6:
        raise ValueError("Sensitivity freeze requires six unique dataset records")
    primary = json.loads((PAPER_ROOT / "experiments/protocol_v1/source_freeze.json").read_text(encoding="utf-8"))
    planned = [
        {"dataset": row["name"], "variant": variant, "seed": seed}
        for row in rows for variant in SENSITIVITY_VARIANTS for seed in SEEDS
    ]
    payload = {
        "version": VERSION,
        "parent_primary_freeze_hash": primary["freeze_hash"],
        "dataset_manifest_sha256": sha256_file(args.manifest),
        "dataset_sha256": {row["name"]: row["sha256"] for row in rows},
        "dataset_n_cells": {row["name"]: int(row["n_cells"]) for row in rows},
        "dataset_n_labels": {row["name"]: int(row["n_labels"]) for row in rows},
        "variants": list(SENSITIVITY_VARIANTS),
        "seeds": list(SEEDS),
        "planned_matrix": planned,
        "expected_runs": 54,
        "code_sha256": sensitivity_code_digests(),
        "runtime": runtime_identity(),
    }
    freeze = {"created_utc": utc_now(), "freeze_hash": sha256_payload(payload), **payload}
    if args.output.exists():
        existing = json.loads(args.output.read_text(encoding="utf-8"))
        if existing["freeze_hash"] != freeze["freeze_hash"]:
            raise FileExistsError(f"Refusing to overwrite different sensitivity freeze: {args.output}")
    else:
        write_json(args.output, freeze)
    planned_path = args.output.parent / "planned_matrix.json"
    if planned_path.exists():
        if json.loads(planned_path.read_text(encoding="utf-8")) != planned:
            raise FileExistsError(f"Existing sensitivity plan differs: {planned_path}")
    else:
        write_json(planned_path, planned)
    print(freeze["freeze_hash"])


if __name__ == "__main__":
    main()
