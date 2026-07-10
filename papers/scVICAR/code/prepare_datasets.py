from __future__ import annotations

import argparse
import csv
import json
import tempfile
from dataclasses import asdict
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

from .config import DATASETS, PROJECT_ROOT, REMOTE_DATA_ROOT
from .io_utils import sha256_file, utc_now, write_json
from .remote_store import RemoteStore


def prepare_dataset(name: str, output_path: Path) -> dict:
    spec = DATASETS[name]
    source = (PROJECT_ROOT / spec.source_path).resolve()
    if not source.is_file():
        raise FileNotFoundError(source)
    adata = ad.read_h5ad(source)
    if spec.source_label_key not in adata.obs:
        raise KeyError(f"{spec.source_label_key!r} absent from {name}: {list(adata.obs.columns)}")

    labels = adata.obs[spec.source_label_key].astype(str).str.strip()
    excluded = {item.casefold() for item in spec.exclude_labels}
    eligible = ~labels.str.casefold().isin(excluded)
    counts = labels[eligible].value_counts()
    retained_types = set(counts[counts >= spec.min_class_size].index.astype(str))
    keep = eligible & labels.isin(retained_types)
    filtered = adata[keep.to_numpy()].copy()
    filtered_labels = labels[keep].astype(str).to_numpy()
    filtered.obs["resolved_label"] = pd.Categorical(filtered_labels)
    filtered.uns["scvicar_dataset_protocol"] = {
        "dataset": name,
        "source_label_key": spec.source_label_key,
        "excluded_labels": list(spec.exclude_labels),
        "minimum_class_size": int(spec.min_class_size),
        "source_cells": int(adata.n_obs),
        "retained_cells": int(filtered.n_obs),
        "retained_clusters": int(pd.Series(filtered_labels).nunique()),
        "created_utc": utc_now(),
    }
    if filtered.n_obs == 0:
        raise ValueError(f"Filtering removed every cell from {name}")
    actual_clusters = int(filtered.obs["resolved_label"].nunique())
    if actual_clusters != spec.expected_clusters:
        raise ValueError(f"{name}: expected {spec.expected_clusters} clusters, found {actual_clusters}")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    filtered.write_h5ad(output_path, compression="gzip")
    return {
        **asdict(spec),
        "source_path_resolved": str(source),
        "canonical_filename": output_path.name,
        "created_utc": filtered.uns["scvicar_dataset_protocol"]["created_utc"],
        "source_cells": int(adata.n_obs),
        "retained_cells": int(filtered.n_obs),
        "n_genes": int(filtered.n_vars),
        "retained_clusters": actual_clusters,
        "file_size_bytes": int(output_path.stat().st_size),
        "sha256": sha256_file(output_path),
        "remote_path": f"{REMOTE_DATA_ROOT}/datasets/confirmatory_v1/{output_path.name}",
    }


def write_manifest(rows: list[dict], output_dir: Path) -> None:
    write_json(output_dir / "dataset_manifest.json", rows)
    columns = sorted({key for row in rows for key in row})
    with (output_dir / "dataset_manifest.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: json.dumps(value, ensure_ascii=False) if isinstance(value, (list, dict, tuple)) else value for key, value in row.items()})


def upload_one(name: str, local_path: Path, row: dict, store: RemoteStore) -> None:
    remote_path = row["remote_path"]
    if store.exists(remote_path):
        remote_digest = store.run(f"sha256sum {remote_path}").stdout.split()[0]
        if remote_digest != row["sha256"]:
            raise FileExistsError(f"Immutable remote dataset differs: {remote_path}")
        return
    store.upload_file(local_path, remote_path, immutable=True)


def main() -> None:
    parser = argparse.ArgumentParser(description="Prepare label-filtered confirmatory scVICAR datasets")
    parser.add_argument("--datasets", nargs="+", default=list(DATASETS), choices=sorted(DATASETS))
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--upload", action="store_true")
    parser.add_argument("--keep-local", action="store_true")
    args = parser.parse_args()

    if args.output_dir is None:
        output_dir = Path(tempfile.mkdtemp(prefix="scvicar-datasets-"))
    else:
        output_dir = args.output_dir.resolve()
        output_dir.mkdir(parents=True, exist_ok=True)
    store = RemoteStore() if args.upload else None
    if store:
        store.ensure_layout()

    rows: list[dict] = []
    for name in args.datasets:
        output_path = output_dir / f"{name}.h5ad"
        row = prepare_dataset(name, output_path)
        rows.append(row)
        if store:
            upload_one(name, output_path, row, store)
        if args.upload and not args.keep_local:
            output_path.unlink()

    write_manifest(rows, output_dir)
    if store:
        for filename in ("dataset_manifest.json", "dataset_manifest.csv"):
            local = output_dir / filename
            remote = f"{REMOTE_DATA_ROOT}/manifests/{filename}"
            if store.exists(remote):
                existing = store.run(f"sha256sum {remote}").stdout.split()[0]
                if existing != sha256_file(local):
                    versioned = f"{REMOTE_DATA_ROOT}/manifests/{Path(filename).stem}-{utc_now().replace(':', '')}{Path(filename).suffix}"
                    store.upload_file(local, versioned, immutable=True)
            else:
                store.upload_file(local, remote, immutable=True)
    print(output_dir)


if __name__ == "__main__":
    main()

