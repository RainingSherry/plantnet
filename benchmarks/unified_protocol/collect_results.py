#!/usr/bin/env python
import argparse
import glob
import os

import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(description="Collect unified benchmark CSV files.")
    parser.add_argument("--root", required=True)
    parser.add_argument("--out_csv", required=True)
    return parser.parse_args()


def main():
    args = parse_args()
    frames = []
    for path in glob.glob(os.path.join(args.root, "**", "*.csv"), recursive=True):
        if os.path.basename(path) == os.path.basename(args.out_csv):
            continue
        try:
            df = pd.read_csv(path)
        except Exception:
            continue
        required = {"dataset", "method", "cluster_method"}
        if required.issubset(df.columns):
            frames.append(df)
    if not frames:
        raise SystemExit("No result CSV files found.")
    out = pd.concat(frames, ignore_index=True)
    out = out.drop_duplicates(subset=["dataset", "method", "cluster_method"], keep="last")
    out = out.sort_values(["dataset", "method", "cluster_method"])
    os.makedirs(os.path.dirname(args.out_csv), exist_ok=True)
    out.to_csv(args.out_csv, index=False)
    print(out.to_string(index=False))


if __name__ == "__main__":
    main()
