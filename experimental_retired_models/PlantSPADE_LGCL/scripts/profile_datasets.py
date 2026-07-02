#!/usr/bin/env python
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd
import yaml

SCRIPT_DIR = Path(__file__).resolve().parent
ROOT = next(parent for parent in [SCRIPT_DIR, *SCRIPT_DIR.parents] if (parent / "methods" / "DeepLearning" / "scMAE_family.py").exists())
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))


def ensure_dir(path):
    Path(path).mkdir(parents=True, exist_ok=True)
    return str(path)


def save_json(payload, path):
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, default=str)


def parse_args():
    parser = argparse.ArgumentParser(description="Profile configured plant h5ad datasets.")
    parser.add_argument(
        "--datasets_config",
        default=str(SCRIPT_DIR.parent / "configs" / "datasets_8plant.yaml"),
    )
    parser.add_argument(
        "--output_dir",
        default=str(ROOT / "results" / "PlantSPADE_LGCL_protocol" / "dataset_profiles"),
    )
    parser.add_argument("--strict", action="store_true", help="Fail on missing/unreadable datasets.")
    return parser.parse_args()


def main():
    args = parse_args()
    with open(args.datasets_config, "r", encoding="utf-8") as handle:
        cfg = yaml.safe_load(handle)

    output_dir = Path(ensure_dir(args.output_dir))
    rows = []
    for entry in cfg.get("datasets", []):
        dataset_name = entry["dataset_name"]
        file_path = Path(entry["file_path"])
        dataset_dir = Path(ensure_dir(output_dir / dataset_name))
        if not file_path.exists():
            payload = {
                **entry,
                "status": "missing",
                "error": f"File not found: {file_path}",
            }
            save_json(payload, str(dataset_dir / "dataset_profile.json"))
            rows.append(payload)
            print(f"[missing] {dataset_name}: {file_path}", flush=True)
            if args.strict:
                raise FileNotFoundError(file_path)
            continue

        try:
            print(f"[read] {dataset_name}: {file_path}", flush=True)
            import scanpy as sc

            from experimental_retired_models.PlantSPADE_LGCL.data import profile_anndata

            adata = sc.read_h5ad(file_path)
            profile = profile_anndata(
                adata,
                dataset_name=dataset_name,
                label_key=entry.get("label_key"),
                input_mode="auto",
            )
            payload = {**entry, **profile, "status": "ok"}
            save_json(payload, str(dataset_dir / "dataset_profile.json"))
            rows.append(payload)
            print(f"[ok] {dataset_name}: {profile['n_cells']} cells, {profile['n_genes']} genes", flush=True)
        except Exception as exc:
            payload = {**entry, "status": "error", "error": str(exc)}
            save_json(payload, str(dataset_dir / "dataset_profile.json"))
            rows.append(payload)
            print(f"[error] {dataset_name}: {exc}", flush=True)
            if args.strict:
                raise

    summary = pd.json_normalize(rows)
    summary.to_csv(output_dir / "dataset_profiles_summary.csv", index=False)
    print(f"Saved profiles to {output_dir}", flush=True)


if __name__ == "__main__":
    main()
