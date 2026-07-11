from __future__ import annotations

import argparse
import json
import shutil
from pathlib import Path

from .config import PAPER_ROOT, PROJECT_ROOT, REMOTE_RESULT_ROOT, sha256_payload
from .io_utils import (
    git_revision,
    sha256_file,
    utc_now,
    verify_checksum_manifest,
    write_checksum_manifest,
    write_json,
)
from .remote_store import RemoteStore


def selected_files() -> list[Path]:
    files: set[Path] = set()
    manuscript = PAPER_ROOT / "manuscript"
    for pattern in ("*.tex", "*.bib", "sections/*.tex", "generated/*.tex"):
        files.update(path for path in manuscript.glob(pattern) if path.is_file())
    files.update(
        path for path in (PAPER_ROOT / "manuscript_zh").glob("*.md")
        if path.is_file()
    )
    for name in (
        "full_manuscript_argument_audit.md",
        "terminology_ledger.md",
        "terminology_provenance_audit.md",
    ):
        path = PAPER_ROOT / "planning" / name
        if path.is_file():
            files.add(path)
    files.update(
        path for path in (
            manuscript / "scVICAR_draft.pdf",
            manuscript / "scVICAR_supplement_draft.pdf",
            PAPER_ROOT / "manifests/release_audit.json",
            PAPER_ROOT / "manifests/release_audit_remote.json",
        ) if path.is_file()
    )
    files.update((PAPER_ROOT / "figures/final").glob("fig[1-6]_*.pdf"))
    files.update((PAPER_ROOT / "figures/final").glob("fig[1-6]_*.svg"))
    files.update((PAPER_ROOT / "tables/generated").glob("*.tex"))
    files.update((PAPER_ROOT / "tables/development_valid15").glob("*.tex"))
    files.update((PAPER_ROOT / "tables/development_valid15").glob("*.csv"))
    for path in (PAPER_ROOT / "figures/source_data").glob("*/SHA256SUMS"):
        files.add(path)
        files.update(item for item in path.parent.glob("*.csv") if item.is_file())
    aggregate_files = (
        "experiments/protocol_v1/run_master.csv",
        "experiments/baselines_v1/run_master.csv",
        "experiments/baselines_v1/superseded_or_invalid_runs.csv",
        "experiments/stress_v1/stress_runs.csv",
        "experiments/downstream_v1/dataset_variant_metrics.csv",
        "experiments/downstream_v1/downstream_contrasts.csv",
        "experiments/leiden_fixed_v1/aggregate/variant_overall_metrics.csv",
        "experiments/leiden_fixed_v1/aggregate/contrasts.csv",
        "experiments/sensitivity_full_labels_v1/aggregate/variant_overall_metrics.csv",
        "experiments/sensitivity_full_labels_v1/aggregate/contrasts.csv",
    )
    files.update(path for rel in aggregate_files if (path := PAPER_ROOT / rel).is_file())
    return sorted(files)


def main() -> None:
    parser = argparse.ArgumentParser(description="Create a content-addressed scVICAR paper snapshot")
    parser.add_argument("--staging-root", type=Path, default=PAPER_ROOT / ".staging/paper_snapshots")
    parser.add_argument("--upload", action="store_true")
    args = parser.parse_args()
    sources = selected_files()
    if not sources:
        raise RuntimeError("Paper snapshot selection is empty")
    for name in ("release_audit.json", "release_audit_remote.json"):
        audit = json.loads((PAPER_ROOT / "manifests" / name).read_text(encoding="utf-8"))
        if not audit.get("passed"):
            raise RuntimeError(f"Refusing snapshot with failed audit: {name}")
    source_hashes = {
        path.relative_to(PAPER_ROOT).as_posix(): sha256_file(path) for path in sources
    }
    snapshot_id = sha256_payload(source_hashes)
    target = args.staging_root / snapshot_id
    if target.exists():
        verify_checksum_manifest(target)
    else:
        for source in sources:
            destination = target / source.relative_to(PAPER_ROOT)
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source, destination)
        write_json(target / "snapshot_metadata.json", {
            "snapshot_id": snapshot_id,
            "created_utc": utc_now(),
            "git": git_revision(PROJECT_ROOT),
            "source_files": len(source_hashes),
            "source_sha256": source_hashes,
            "release_audit_passed": True,
        })
        (target / "COMPLETED").write_text("paper snapshot complete\n", encoding="utf-8")
        write_checksum_manifest(target)
        verify_checksum_manifest(target)
    remote = f"{REMOTE_RESULT_ROOT}/paper_snapshots/{snapshot_id}"
    if args.upload:
        store = RemoteStore(); store.ensure_layout(); store.upload_directory_atomic(target, remote)
        result = store.run(f"cd {remote} && sha256sum -c SHA256SUMS", check=False)
        if result.returncode != 0:
            raise RuntimeError(f"Remote paper snapshot checksum failed: {remote}")
    print(f"snapshot_id={snapshot_id} files={len(source_hashes)} remote={remote}")


if __name__ == "__main__":
    main()
