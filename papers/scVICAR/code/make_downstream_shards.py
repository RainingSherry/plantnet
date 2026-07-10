from __future__ import annotations

import argparse
import json
from pathlib import Path, PurePosixPath

import pandas as pd

from .config import PAPER_ROOT, PROTOCOL_VERSION, REMOTE_RESULT_ROOT
from .remote_store import RemoteStore


def main() -> None:
    parser = argparse.ArgumentParser(description="Shard only incomplete frozen downstream run IDs")
    parser.add_argument(
        "--run-master", type=Path,
        default=PAPER_ROOT / f"experiments/{PROTOCOL_VERSION}/run_master.csv",
    )
    parser.add_argument(
        "--output-dir", type=Path,
        default=PAPER_ROOT / f"experiments/downstream_v1/shards",
    )
    parser.add_argument("--n-shards", type=int, default=4)
    args = parser.parse_args()
    if args.n_shards < 1:
        raise ValueError("--n-shards must be positive")
    frame = pd.read_csv(args.run_master)
    formal = frame[frame["execution_mode"] == "formal"].copy()
    if len(formal) != 108 or formal["run_id"].nunique() != 108:
        raise RuntimeError("Downstream sharding requires the exact 108-run formal master")
    store = RemoteStore()
    root = f"{REMOTE_RESULT_ROOT}/downstream/{PROTOCOL_VERSION}"
    found = store.run(f"find {root} -type f -name COMPLETED -print", check=False)
    completed = {
        PurePosixPath(line.strip()).parent.name
        for line in found.stdout.splitlines() if line.strip()
    }
    known = set(formal["run_id"])
    if not completed.issubset(known):
        raise ValueError(f"Remote downstream namespace contains {len(completed - known)} unknown run IDs")
    pending = formal[~formal["run_id"].isin(completed)].sort_values(
        ["dataset", "variant", "seed", "run_id"]
    ).reset_index(drop=True)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    manifest = {
        "formal_runs": 108, "remote_complete": len(completed), "pending": len(pending),
        "n_shards": args.n_shards, "shards": [],
    }
    for index in range(args.n_shards):
        shard = pending.iloc[index::args.n_shards].copy()
        path = args.output_dir / f"shard_{index:02d}.csv"
        shard.to_csv(path, index=False)
        manifest["shards"].append({
            "index": index, "path": str(path), "runs": len(shard),
            "run_ids": shard["run_id"].tolist(),
        })
    (args.output_dir / "shard_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps({key: manifest[key] for key in ("remote_complete", "pending", "n_shards")}, indent=2))


if __name__ == "__main__":
    main()
