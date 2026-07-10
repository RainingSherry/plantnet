from __future__ import annotations

import argparse
import json
import shutil
import tempfile
from pathlib import Path

import anndata as ad
import pandas as pd

from .config import DATASETS, PAPER_ROOT, PROJECT_ROOT, REMOTE_DATA_ROOT
from .io_utils import sha256_file, utc_now, write_json
from .remote_store import RemoteStore


VERSION = "sensitivity_full_labels_v1"


def prepare_one(name: str, output: Path) -> dict:
    spec = DATASETS[name]
    source = (PROJECT_ROOT / spec.source_path).resolve()
    if not source.is_file():
        raise FileNotFoundError(source)
    adata = ad.read_h5ad(source)
    if spec.source_label_key not in adata.obs:
        raise KeyError(f"{spec.source_label_key!r} absent from {name}")
    labels = adata.obs[spec.source_label_key].astype(str).str.strip()
    adata.obs["resolved_label"] = pd.Categorical(labels.to_numpy())
    counts = labels.value_counts(dropna=False).sort_index()
    created = utc_now()
    adata.uns["scvicar_full_label_sensitivity"] = {
        "version": VERSION,
        "dataset": name,
        "source_label_key": spec.source_label_key,
        "filtering": "none; every source cell and string-valued source label retained",
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "n_labels": int(labels.nunique(dropna=False)),
        "created_utc": created,
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(output, compression="gzip")
    return {
        "name": name,
        "version": VERSION,
        "source_path": str(source),
        "source_label_key": spec.source_label_key,
        "canonical_filename": output.name,
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "n_labels": int(labels.nunique(dropna=False)),
        "label_counts": {str(key): int(value) for key, value in counts.items()},
        "filtering": "none",
        "created_utc": created,
        "file_size_bytes": int(output.stat().st_size),
        "sha256": sha256_file(output),
        "remote_path": f"{REMOTE_DATA_ROOT}/datasets/{VERSION}/{output.name}",
    }


def upload_immutable(store: RemoteStore, local: Path, remote: str, digest: str) -> None:
    if store.exists(remote):
        observed = store.run(f"sha256sum {remote}").stdout.split()[0]
        if observed != digest:
            raise FileExistsError(f"Remote sensitivity dataset differs: {remote}")
        return
    store.upload_file(local, remote, immutable=True)


def main() -> None:
    parser = argparse.ArgumentParser(description="Prepare unfiltered full-label sensitivity H5AD files")
    parser.add_argument("--datasets", nargs="+", default=list(DATASETS), choices=sorted(DATASETS))
    parser.add_argument(
        "--output-dir", type=Path,
        default=PAPER_ROOT / f"manifests/{VERSION}",
    )
    parser.add_argument("--upload", action="store_true")
    parser.add_argument("--keep-local", action="store_true")
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    store = RemoteStore() if args.upload else None
    rows = []
    for name in args.datasets:
        with tempfile.TemporaryDirectory(prefix=f"scvicar-{VERSION}-{name}-") as temp:
            local = Path(temp) / f"{name}.h5ad"
            row = prepare_one(name, local)
            rows.append(row)
            if store:
                upload_immutable(store, local, row["remote_path"], row["sha256"])
            if args.keep_local:
                target = args.output_dir / local.name
                shutil.copy2(local, target)
    manifest = args.output_dir / "dataset_manifest.json"
    write_json(manifest, rows)
    if store:
        remote = f"{REMOTE_DATA_ROOT}/manifests/{VERSION}_dataset_manifest_{sha256_file(manifest)[:16]}.json"
        upload_immutable(store, manifest, remote, sha256_file(manifest))
    print(json.dumps({"version": VERSION, "datasets": len(rows), "manifest": str(manifest)}, indent=2))


if __name__ == "__main__":
    main()
