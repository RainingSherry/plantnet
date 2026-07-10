from __future__ import annotations

import argparse
import json
from pathlib import Path

from .config import PAPER_ROOT, sha256_payload
from .io_utils import sha256_file, utc_now, write_json
from .run_scVICAR import runtime_identity
from .run_stress import STRESS_VERSION, stress_code_digests


def main() -> None:
    parser = argparse.ArgumentParser(description="Freeze stress_v1 code, environment, plan, and primary parent")
    parser.add_argument("--output", type=Path, default=PAPER_ROOT / "experiments/stress_v1/source_freeze.json")
    args = parser.parse_args()
    primary = json.loads((PAPER_ROOT / "experiments/protocol_v1/source_freeze.json").read_text(encoding="utf-8"))
    plan = PAPER_ROOT / "experiments/stress_v1/planned_matrix.json"
    payload = {
        "stress_version": STRESS_VERSION,
        "parent_primary_freeze_hash": primary["freeze_hash"],
        "code_sha256": stress_code_digests(),
        "runtime": runtime_identity(),
        "planned_matrix_sha256": sha256_file(plan),
        "expected_runs": 126,
    }
    freeze = {"created_utc": utc_now(), "freeze_hash": sha256_payload(payload), **payload}
    if args.output.exists():
        existing = json.loads(args.output.read_text(encoding="utf-8"))
        if existing["freeze_hash"] != freeze["freeze_hash"]:
            raise FileExistsError(f"Refusing to overwrite a different stress freeze: {args.output}")
    else:
        write_json(args.output, freeze)
    print(freeze["freeze_hash"])


if __name__ == "__main__":
    main()
