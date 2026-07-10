from __future__ import annotations

import argparse
import json
from pathlib import Path

from .config import PAPER_ROOT, sha256_payload
from .io_utils import utc_now, write_json
from .run_baseline import BASELINES, BASELINE_VERSION, baseline_identity


def build_freeze() -> dict:
    primary = json.loads(
        (PAPER_ROOT / "experiments/protocol_v1/source_freeze.json").read_text(encoding="utf-8")
    )
    payload = {
        "baseline_version": BASELINE_VERSION,
        "primary_freeze_hash": primary["freeze_hash"],
        "dataset_sha256": primary["dataset_sha256"],
        "methods": {method: baseline_identity(method) for method in sorted(BASELINES)},
    }
    return {"created_utc": utc_now(), "freeze_hash": sha256_payload(payload), **payload}


def main() -> None:
    parser = argparse.ArgumentParser(description="Freeze external baseline code and environments")
    parser.add_argument(
        "--output", type=Path,
        default=PAPER_ROOT / "experiments/baselines_v1/source_freeze.json",
    )
    args = parser.parse_args()
    freeze = build_freeze()
    if args.output.exists():
        existing = json.loads(args.output.read_text(encoding="utf-8"))
        if existing.get("freeze_hash") != freeze["freeze_hash"]:
            raise FileExistsError(f"Refusing to overwrite a different baseline freeze: {args.output}")
        print(existing["freeze_hash"])
        return
    write_json(args.output, freeze)
    print(freeze["freeze_hash"])


if __name__ == "__main__":
    main()
