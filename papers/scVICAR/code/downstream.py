from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.optimize import linear_sum_assignment
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, balanced_accuracy_score, f1_score, recall_score
from sklearn.model_selection import train_test_split
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

from .config import SPLIT_SEEDS
from .io_utils import utc_now, write_json


def dense_vector(values: Any) -> np.ndarray:
    return np.asarray(values.toarray() if sp.issparse(values) else values).ravel()


def expression_for_de(adata: ad.AnnData) -> ad.AnnData:
    work = adata.copy()
    matrix = work.layers["counts"] if "counts" in work.layers else (work.raw.X if work.raw is not None else work.X)
    var = work.raw.var.copy() if work.raw is not None and matrix.shape[1] == work.raw.n_vars else work.var.copy()
    work = ad.AnnData(X=matrix.copy(), obs=work.obs.copy(), var=var)
    sample = work.X[: min(256, work.n_obs)]
    vals = sample.data if sp.issparse(sample) else np.asarray(sample).ravel()
    vals = vals[np.isfinite(vals)]
    raw_like = bool(vals.size and np.all(vals >= 0) and np.allclose(vals, np.round(vals), atol=1e-4))
    if raw_like:
        totals = dense_vector(work.X.sum(axis=1)).astype(np.float64)
        scale = np.divide(10000.0, totals, out=np.zeros_like(totals), where=totals > 0)
        if sp.issparse(work.X):
            work.X = work.X.multiply(scale[:, None]).tocsr()
            work.X.data = np.log1p(work.X.data)
        else:
            work.X = np.log1p(np.asarray(work.X) * scale[:, None])
    expressed = dense_vector((work.X > 0).sum(axis=0)) >= 3
    return work[:, expressed].copy()


def rank_markers(work: ad.AnnData, indices: np.ndarray, groups: np.ndarray, n_genes: int = 100) -> dict[str, pd.DataFrame]:
    import scanpy as sc

    subset = work[indices].copy()
    subset.obs["_group"] = pd.Categorical(groups.astype(str))
    output: dict[str, pd.DataFrame] = {}
    counts = subset.obs["_group"].value_counts()
    valid = counts[counts >= 3].index.astype(str).tolist()
    if len(valid) < 2:
        return output
    subset = subset[subset.obs["_group"].astype(str).isin(valid)].copy()
    sc.tl.rank_genes_groups(subset, "_group", method="wilcoxon", n_genes=min(n_genes, subset.n_vars), tie_correct=True)
    for group in valid:
        frame = sc.get.rank_genes_groups_df(subset, group=group)
        frame = frame.replace([np.inf, -np.inf], np.nan).dropna(subset=["names"])
        frame["names"] = frame["names"].astype(str)
        output[group] = frame.drop_duplicates("names").head(n_genes).reset_index(drop=True)
    return output


def marker_sets(markers: dict[str, pd.DataFrame], n: int) -> dict[str, list[str]]:
    return {key: frame["names"].astype(str).head(n).tolist() for key, frame in markers.items()}


def overlap_matrix(reference: dict[str, pd.DataFrame], predicted: dict[str, pd.DataFrame], n: int = 100) -> pd.DataFrame:
    ref = marker_sets(reference, n)
    pred = marker_sets(predicted, n)
    return pd.DataFrame(
        [[len(set(ref[r]) & set(pred[p])) for p in pred] for r in ref],
        index=list(ref), columns=list(pred), dtype=float,
    )


def matched_marker_metrics(reference: dict[str, pd.DataFrame], predicted: dict[str, pd.DataFrame]) -> dict[str, float]:
    matrix100 = overlap_matrix(reference, predicted, 100)
    if matrix100.empty:
        return {"recovery_at_20": 0.0, "recovery_at_50": 0.0, "recovery_at_100": 0.0, "jaccard_at_100": 0.0, "overlap_coefficient_at_100": 0.0}
    rows, cols = linear_sum_assignment(-matrix100.to_numpy())
    pairs = [(matrix100.index[i], matrix100.columns[j]) for i, j in zip(rows, cols)]
    # Missing reference types represent failed recovery and must contribute zero
    # rather than silently disappearing from the mean when fewer clusters have
    # valid markers.
    denominator = max(1, len(reference))
    metrics: dict[str, float] = {}
    for n in (20, 50, 100):
        ref = marker_sets(reference, n)
        pred = marker_sets(predicted, n)
        values = [len(set(ref[r]) & set(pred[p])) / float(n) for r, p in pairs]
        metrics[f"recovery_at_{n}"] = float(np.sum(values) / denominator) if values else 0.0
    ref100 = marker_sets(reference, 100)
    pred100 = marker_sets(predicted, 100)
    jaccard = []
    overlap = []
    for r, p in pairs:
        a, b = set(ref100[r]), set(pred100[p])
        inter = len(a & b)
        jaccard.append(inter / max(1, len(a | b)))
        overlap.append(inter / max(1, min(len(a), len(b))))
    metrics["jaccard_at_100"] = float(np.sum(jaccard) / denominator)
    metrics["overlap_coefficient_at_100"] = float(np.sum(overlap) / denominator)
    return metrics


def annotate_clusters(reference: dict[str, pd.DataFrame], predicted: dict[str, pd.DataFrame]) -> dict[str, str]:
    matrix = overlap_matrix(reference, predicted, 100)
    mapping: dict[str, str] = {}
    for cluster in matrix.columns:
        frame = predicted.get(str(cluster))
        if frame is None or len(frame) < 3 or matrix[cluster].max() <= 0:
            mapping[str(cluster)] = "unassigned"
            continue
        candidates = matrix.index[matrix[cluster] == matrix[cluster].max()].astype(str).tolist()
        if len(candidates) == 1:
            mapping[str(cluster)] = candidates[0]
            continue
        cluster_scores = dict(zip(frame["names"].astype(str), frame.get("logfoldchanges", pd.Series(np.zeros(len(frame))))))
        tie_scores = {}
        for cell_type in candidates:
            shared = set(reference[cell_type]["names"].astype(str)) & set(cluster_scores)
            values = np.asarray([cluster_scores[g] for g in shared], dtype=float)
            values = values[np.isfinite(values)]
            tie_scores[cell_type] = float(values.mean()) if values.size else -np.inf
        mapping[str(cluster)] = sorted(candidates, key=lambda item: (-tie_scores[item], item))[0]
    return mapping


def oracle_mapping(y_true: np.ndarray, clusters: np.ndarray) -> dict[str, str]:
    true_values = np.unique(y_true.astype(str))
    cluster_values = np.unique(clusters.astype(str))
    counts = np.zeros((len(true_values), len(cluster_values)), dtype=int)
    for i, label in enumerate(true_values):
        for j, cluster in enumerate(cluster_values):
            counts[i, j] = int(np.sum((y_true.astype(str) == label) & (clusters.astype(str) == cluster)))
    rows, cols = linear_sum_assignment(-counts)
    return {str(cluster_values[j]): str(true_values[i]) for i, j in zip(rows, cols)}


def annotation_metrics(y_true: np.ndarray, y_pred: np.ndarray) -> dict[str, Any]:
    labels = sorted(np.unique(y_true.astype(str)).tolist())
    return {
        "accuracy": float(accuracy_score(y_true, y_pred)),
        "balanced_accuracy": float(balanced_accuracy_score(y_true, y_pred)),
        "macro_f1": float(f1_score(y_true, y_pred, labels=labels, average="macro", zero_division=0)),
        "coverage": float(np.mean(y_pred.astype(str) != "unassigned")),
        "per_class_recall": {
            label: float(value)
            for label, value in zip(labels, recall_score(y_true, y_pred, labels=labels, average=None, zero_division=0))
        },
    }


def run_marker_tasks(adata: ad.AnnData, clusters: np.ndarray, split_seed: int, output_dir: Path) -> dict:
    labels = adata.obs["resolved_label"].astype(str).to_numpy()
    all_indices = np.arange(adata.n_obs)
    reference_idx, evaluation_idx = train_test_split(
        all_indices, train_size=0.5, stratify=labels, random_state=split_seed,
    )
    if np.intersect1d(reference_idx, evaluation_idx).size or len(reference_idx) + len(evaluation_idx) != adata.n_obs:
        raise AssertionError("Reference/evaluation split is not a disjoint partition")
    work = expression_for_de(adata)
    reference = rank_markers(work, reference_idx, labels[reference_idx])
    predicted = rank_markers(work, evaluation_idx, clusters[evaluation_idx].astype(str))
    matrix = overlap_matrix(reference, predicted, 100)
    output_dir.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        output_dir / "split_indices.npz",
        reference=np.asarray(reference_idx, dtype=np.int64),
        evaluation=np.asarray(evaluation_idx, dtype=np.int64),
    )
    matrix.to_csv(output_dir / "marker_overlap.csv")
    for key, frame in reference.items():
        frame.to_csv(output_dir / f"reference_markers__{key}.csv", index=False)
    for key, frame in predicted.items():
        frame.to_csv(output_dir / f"cluster_markers__{key}.csv", index=False)

    recovery = matched_marker_metrics(reference, predicted)
    evaluation_clusters = sorted(np.unique(clusters[evaluation_idx].astype(str)).tolist())
    invalid = [cluster not in predicted or len(predicted[cluster]) < 3 for cluster in evaluation_clusters]
    recovery["invalid_cluster_fraction"] = float(np.mean(invalid)) if invalid else 1.0
    mapping = annotate_clusters(reference, predicted)
    predicted_types = np.asarray([mapping.get(str(cluster), "unassigned") for cluster in clusters[evaluation_idx]], dtype=object)
    annotation = annotation_metrics(labels[evaluation_idx], predicted_types)
    oracle = oracle_mapping(labels[evaluation_idx], clusters[evaluation_idx])
    oracle_types = np.asarray([oracle.get(str(cluster), "unassigned") for cluster in clusters[evaluation_idx]], dtype=object)
    oracle_scores = annotation_metrics(labels[evaluation_idx], oracle_types)
    payload = {
        "split_seed": split_seed,
        "reference_cells": int(len(reference_idx)),
        "evaluation_cells": int(len(evaluation_idx)),
        "split_disjoint": True,
        "split_complete": True,
        "marker_recovery": recovery,
        "marker_annotation": annotation,
        "oracle_upper_bound": oracle_scores,
        "cluster_annotation": mapping,
    }
    write_json(output_dir / "marker_results.json", payload)
    return payload


def run_linear_probe(embedding: np.ndarray, labels: np.ndarray, fraction: float, split_seed: int) -> dict:
    indices = np.arange(len(labels))
    train_idx, test_idx = train_test_split(indices, train_size=fraction, stratify=labels, random_state=split_seed)
    classifier = make_pipeline(
        StandardScaler(),
        # lbfgs optimizes multinomial loss for multiclass data in current
        # scikit-learn releases; the removed multi_class argument is omitted so
        # the fixed protocol works across supported versions.
        LogisticRegression(C=1.0, class_weight="balanced", max_iter=2000, solver="lbfgs"),
    )
    classifier.fit(embedding[train_idx], labels[train_idx])
    pred = classifier.predict(embedding[test_idx])
    values = annotation_metrics(labels[test_idx], pred)
    values.update({"label_fraction": fraction, "split_seed": split_seed, "train_cells": len(train_idx), "test_cells": len(test_idx)})
    return values


def main() -> None:
    parser = argparse.ArgumentParser(description="Run preregistered scVICAR downstream tasks")
    parser.add_argument("--data-path", required=True, type=Path)
    parser.add_argument("--run-dir", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--task", nargs="+", default=["marker_recovery", "marker_annotation", "linear_probe"], choices=["marker_recovery", "marker_annotation", "linear_probe"])
    parser.add_argument("--reference-fraction", type=float, default=0.5)
    parser.add_argument("--label-fractions", default="0.1,0.3")
    parser.add_argument("--split-seeds", default=",".join(map(str, SPLIT_SEEDS)))
    args = parser.parse_args()
    if args.reference_fraction != 0.5:
        raise ValueError("protocol_v1 fixes --reference-fraction at 0.5")
    split_seeds = tuple(int(item) for item in args.split_seeds.split(","))
    if split_seeds != SPLIT_SEEDS:
        raise ValueError(f"protocol_v1 split seeds are {SPLIT_SEEDS}")
    fractions = tuple(float(item) for item in args.label_fractions.split(","))
    if fractions != (0.1, 0.3):
        raise ValueError("protocol_v1 label fractions are 0.1 and 0.3")

    adata = ad.read_h5ad(args.data_path)
    if "resolved_label" not in adata.obs:
        raise KeyError("Canonical dataset must contain obs['resolved_label']")
    embedding = np.load(args.run_dir / "embedding_float32.npz")["embedding"].astype(np.float32)
    clusters = np.load(args.run_dir / "clusters.npz")["predicted"].astype(str)
    if embedding.shape[0] != adata.n_obs or len(clusters) != adata.n_obs:
        raise ValueError("Run artifacts and dataset contain different cells")
    labels = adata.obs["resolved_label"].astype(str).to_numpy()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    summary: dict[str, Any] = {"created_utc": utc_now(), "tasks": args.task, "marker_splits": [], "linear_probe": []}
    if {"marker_recovery", "marker_annotation"}.intersection(args.task):
        for seed in split_seeds:
            summary["marker_splits"].append(run_marker_tasks(adata, clusters, seed, args.output_dir / f"split_{seed}"))
    if "linear_probe" in args.task:
        for fraction in fractions:
            for seed in split_seeds:
                summary["linear_probe"].append(run_linear_probe(embedding, labels, fraction, seed))
        pd.json_normalize(summary["linear_probe"]).to_csv(args.output_dir / "linear_probe.csv", index=False)
    write_json(args.output_dir / "downstream_summary.json", summary)


if __name__ == "__main__":
    main()
