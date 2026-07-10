from __future__ import annotations

import argparse
import importlib.metadata
import json
import sys
from pathlib import Path

from .config import PAPER_ROOT, PROTOCOL_VERSION, sha256_payload
from .io_utils import sha256_file, utc_now, write_json
from .secondary_evaluation import SECONDARY_VERSION


PACKAGES = ("anndata", "scanpy", "leidenalg", "igraph", "numpy", "scikit-learn", "scipy")


def main() -> None:
    parser = argparse.ArgumentParser(description="Freeze fixed-Leiden runner and execution environment")
    parser.add_argument(
        "--output", type=Path,
        default=PAPER_ROOT / f"experiments/{SECONDARY_VERSION}/source_freeze.json",
    )
    args = parser.parse_args()
    primary_path = PAPER_ROOT / f"experiments/{PROTOCOL_VERSION}/source_freeze.json"
    master_path = PAPER_ROOT / f"experiments/{PROTOCOL_VERSION}/run_master.csv"
    runner_path = PAPER_ROOT / "code/secondary_evaluation.py"
    primary = json.loads(primary_path.read_text(encoding="utf-8"))
    versions = {}
    for name in PACKAGES:
        try:
            versions[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            versions[name] = None
    payload = {
        "secondary_version": SECONDARY_VERSION,
        "parent_primary_freeze_hash": primary["freeze_hash"],
        "parent_primary_freeze_sha256": sha256_file(primary_path),
        "formal_run_master_sha256": sha256_file(master_path),
        "runner_sha256": sha256_file(runner_path),
        "environment": {
            "python": str(Path(sys.executable).resolve()),
            "python_sha256": sha256_file(Path(sys.executable).resolve()),
            "packages": versions,
        },
        "protocol": {
            "n_neighbors": 15,
            "resolution": 1.0,
            "oracle_sweep": False,
            "known_k": False,
        },
        "expected_runs": 108,
    }
    freeze = {"created_utc": utc_now(), "freeze_hash": sha256_payload(payload), **payload}
    if args.output.exists():
        existing = json.loads(args.output.read_text(encoding="utf-8"))
        if existing["freeze_hash"] != freeze["freeze_hash"]:
            raise FileExistsError(f"Refusing to overwrite a different secondary freeze: {args.output}")
    else:
        write_json(args.output, freeze)
    print(freeze["freeze_hash"])


if __name__ == "__main__":
    main()
