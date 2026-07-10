from __future__ import annotations

import argparse
import json
from pathlib import Path

from .config import PAPER_ROOT, protocol_snapshot, sha256_payload
from .io_utils import git_revision, sha256_file, utc_now, write_json
from .run_scVICAR import code_digests, runtime_identity


def build_freeze(dataset_manifest: Path) -> dict:
    manifest = json.loads(dataset_manifest.read_text(encoding="utf-8"))
    data_sha256 = {str(row["name"]): str(row["sha256"]) for row in manifest}
    payload = {
        "protocol": protocol_snapshot(),
        "code_sha256": code_digests(),
        "runtime": runtime_identity(),
        "dataset_sha256": data_sha256,
        "dataset_manifest_sha256": sha256_file(dataset_manifest),
    }
    return {
        "created_utc": utc_now(),
        "freeze_hash": sha256_payload(payload),
        "git": git_revision(PAPER_ROOT.parents[1]),
        **payload,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Freeze the exact code/data/protocol identity for protocol_v1")
    parser.add_argument(
        "--dataset-manifest",
        type=Path,
        default=PAPER_ROOT / "manifests" / "dataset_upload" / "dataset_manifest.json",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=PAPER_ROOT / "experiments" / "protocol_v1" / "source_freeze.json",
    )
    args = parser.parse_args()
    freeze = build_freeze(args.dataset_manifest)
    if args.output.exists():
        existing = json.loads(args.output.read_text(encoding="utf-8"))
        if existing.get("freeze_hash") != freeze["freeze_hash"]:
            raise FileExistsError(f"Refusing to overwrite a different freeze: {args.output}")
        print(existing["freeze_hash"])
        return
    write_json(args.output, freeze)
    print(freeze["freeze_hash"])


if __name__ == "__main__":
    main()
