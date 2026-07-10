from __future__ import annotations

import argparse
from itertools import product
from pathlib import Path

import numpy as np
import pandas as pd

from .config import PAPER_ROOT


CONTRASTS = [
    ("fixed", "nomix", "F_vs_NoMix"),
    ("topology_full", "nomix", "T_vs_NoMix"),
    ("topology_full", "fixed", "T_vs_F"),
]


def bootstrap_ci(values: np.ndarray, iterations: int, rng: np.random.Generator) -> tuple[float, float]:
    if values.size == 0:
        return np.nan, np.nan
    draws = rng.choice(values, size=(iterations, values.size), replace=True).mean(axis=1)
    return tuple(np.percentile(draws, [2.5, 97.5]).astype(float))


def sign_flip_pvalue(values: np.ndarray, iterations: int, rng: np.random.Generator) -> float:
    if values.size == 0:
        return np.nan
    observed = abs(float(values.mean()))
    if values.size <= 16:
        signs = np.asarray(list(product([-1.0, 1.0], repeat=values.size)))
    else:
        signs = rng.choice([-1.0, 1.0], size=(iterations, values.size))
    null = np.abs((signs * values).mean(axis=1))
    return float((np.sum(null >= observed) + 1) / (len(null) + 1))


def holm_adjust(pvalues: list[float]) -> list[float]:
    result = [np.nan] * len(pvalues)
    valid = [(index, value) for index, value in enumerate(pvalues) if np.isfinite(value)]
    ordered = sorted(valid, key=lambda item: item[1])
    running = 0.0
    m = len(ordered)
    for rank, (index, value) in enumerate(ordered):
        running = max(running, min(1.0, (m - rank) * value))
        result[index] = running
    return result


def contrast_table(frame: pd.DataFrame, metric: str, iterations: int = 10000) -> pd.DataFrame:
    formal = frame[frame["execution_mode"] == "formal"].copy()
    seed_mean = formal.groupby(["dataset", "variant"], as_index=False)[metric].mean()
    pivot = seed_mean.pivot(index="dataset", columns="variant", values=metric)
    rows = []
    rng = np.random.default_rng(20260710)
    for left, right, name in CONTRASTS:
        available = pivot[[left, right]].dropna() if {left, right}.issubset(pivot.columns) else pd.DataFrame()
        delta = available[left].to_numpy() - available[right].to_numpy() if not available.empty else np.array([])
        lo, hi = bootstrap_ci(delta, iterations, rng)
        rows.append(
            {
                "contrast": name,
                "metric": metric,
                "n_datasets": len(delta),
                "mean_delta": float(np.mean(delta)) if len(delta) else np.nan,
                "median_delta": float(np.median(delta)) if len(delta) else np.nan,
                "ci95_low": lo,
                "ci95_high": hi,
                "wins": int(np.sum(delta > 0)),
                "ties": int(np.sum(delta == 0)),
                "losses": int(np.sum(delta < 0)),
                "permutation_p": sign_flip_pvalue(delta, iterations, rng),
            }
        )
    adjusted = holm_adjust([row["permutation_p"] for row in rows])
    for row, value in zip(rows, adjusted):
        row["holm_p"] = value
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-master", type=Path, default=PAPER_ROOT / "experiments" / "protocol_v1" / "run_master.csv")
    parser.add_argument("--output-dir", type=Path, default=PAPER_ROOT / "tables")
    parser.add_argument("--iterations", type=int, default=10000)
    args = parser.parse_args()
    frame = pd.read_csv(args.run_master)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for metric in ["ari", "nmi", "acc", "f1_macro"]:
        contrast_table(frame, metric, args.iterations).to_csv(args.output_dir / f"confirmatory_contrasts_{metric}.csv", index=False)


if __name__ == "__main__":
    main()

