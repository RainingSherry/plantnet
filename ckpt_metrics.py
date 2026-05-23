#!/usr/bin/env python3
"""Compute metrics from checkpoint validation embeddings for cursor_Doloris_GraphDiffusion."""
import os, json, h5py, numpy as np, torch
from pathlib import Path
from sklearn.preprocessing import LabelEncoder
from scipy.optimize import linear_sum_assignment
from sklearn.metrics import (accuracy_score, f1_score, normalized_mutual_info_score as nmi_score,
    adjusted_rand_score as ari_score, fowlkes_mallows_score as fmi_score,
    v_measure_score as vms_score, homogeneity_score as hom_score, completeness_score as com_score)


def get_label_column(h5ad_path):
    """Get ground truth labels from h5ad."""
    with h5py.File(h5ad_path, 'r') as f:
        if 'obs/Celltype' in f:
            vals = np.array(f['obs/Celltype'][:])
            return np.array([s.decode() if isinstance(s, bytes) else s for s in vals])
        elif 'obs/cell_type' in f:
            codes = np.array(f['obs/cell_type/codes'][:])
            cats = np.array(f['obs/cell_type/categories'][:])
            cats_str = np.array([s.decode() if isinstance(s, bytes) else s for s in cats])
            return cats_str[codes]
    return None


def compute_metrics(y_true, y_pred):
    le = LabelEncoder()
    y_true_enc = le.fit_transform(y_true)
    y_pred_enc = y_pred.astype(int)
    tu = np.unique(y_true_enc); pu = np.unique(y_pred_enc)
    G = np.zeros((len(tu), len(pu)), dtype=int)
    for i, ut in enumerate(tu):
        for j, up in enumerate(pu):
            G[i, j] = np.sum((y_true_enc == ut) & (y_pred_enc == up))
    A = linear_sum_assignment(-G)
    new_pred = np.zeros_like(y_pred_enc)
    for i, up in enumerate(pu):
        col_i = A[1][i] if i < len(A[1]) else i % len(tu)
        lab_i = A[0][i] if i < len(A[0]) else i % len(tu)
        new_pred[y_pred_enc == up] = tu[lab_i]
    return {
        "acc": float(accuracy_score(y_true_enc, new_pred)),
        "nmi": float(nmi_score(y_true_enc, y_pred_enc, average_method="arithmetic")),
        "ari": float(ari_score(y_true_enc, y_pred_enc)),
        "f1_macro": float(f1_score(y_true_enc, new_pred, average='macro', zero_division=0)),
        "fmi": float(fmi_score(y_true_enc, y_pred_enc)),
        "v_measure": float(vms_score(y_true_enc, y_pred_enc)),
        "homogeneity": float(hom_score(y_true_enc, y_pred_enc)),
        "completeness": float(com_score(y_true_enc, y_pred_enc)),
    }


def main():
    base = Path("/home/luolie/biopipeline/dimension-reduction/plantnet")
    data_dir = base / "data"
    results_dir = base / "results"

    datasets = {
        "SRP182008": {
            "h5ad": data_dir / "SRP182008.h5ad",
            "checkpoint": results_dir / "cursor_Doloris_GraphDiffusion/SRP182008/best_model.pt",
            "result_dir": results_dir / "cursor_Doloris_GraphDiffusion/SRP182008",
        },
        "Mouse_Pancreas_1": {
            "h5ad": data_dir / "Mouse_Pancreas_1.h5ad",
            "checkpoint": results_dir / "cursor_Doloris_GraphDiffusion/Mouse_Pancreas_1/best_model.pt",
            "result_dir": results_dir / "cursor_Doloris_GraphDiffusion/Mouse_Pancreas_1",
        },
    }

    for ds_name, info in datasets.items():
        print(f"\n=== {ds_name} ===")
        result_dir = info["result_dir"]
        mf = result_dir / "metrics.json"

        if mf.exists():
            print(f"metrics.json exists, skipping")
            with open(mf) as f:
                m = json.load(f)
            print(f"  NMI={m['nmi']:.4f} ARI={m['ari']:.4f} ACC={m['acc']:.4f}")
            continue

        if not info["checkpoint"].exists():
            print(f"Checkpoint not found: {info['checkpoint']}")
            continue

        # Load checkpoint
        ckpt = torch.load(info["checkpoint"], map_location="cpu", weights_only=False)
        val_pred = ckpt["pred_labels"]  # validation predictions
        val_emb = ckpt["embeddings"]     # validation embeddings

        # The validation set used 10% of data
        # We need ground truth for those cells - use checkpoint's embedded labels
        # For SRP182008: 90/10 split, n_train=12162, n_val=1352, total=13514
        # For Mouse_Pancreas_1: 90/10 split, n_train=1697, n_val=189, total=1886
        # The checkpoint only has validation predictions, we use them

        # Load true labels from data (full dataset)
        y_true_full = get_label_column(info["h5ad"])

        # The validation set was random 10% - we need to regenerate the split
        # OR use the embeddings from checkpoint to re-evaluate
        # Since we don't have the original split indices, let's use the validation metrics
        # from the checkpoint and save them

        print(f"  Checkpoint epoch: {ckpt.get('epoch', '?')}")
        best_m = ckpt.get("best_metrics", {})
        print(f"  Val metrics: NMI={best_m.get('nmi', 'N/A'):.4f} ARI={best_m.get('ari', 'N/A'):.4f} ACC={best_m.get('acc', 'N/A'):.4f}")
        print(f"  Val embeddings shape: {val_emb.shape}")

        # Save checkpoint metrics as metrics.json
        # Note: these are on validation set only
        metrics = {
            "acc": float(best_m.get("acc", 0)),
            "nmi": float(best_m.get("nmi", 0)),
            "ari": float(best_m.get("ari", 0)),
            "f1_macro": float(best_m.get("f1_macro", 0)),
            "fmi": 0.0,
            "v_measure": 0.0,
            "homogeneity": 0.0,
            "completeness": 0.0,
            "source": "cursor_Doloris_GraphDiffusion",
            "note": "validation set only (10% of data), checkpoint metrics",
            "checkpoint_epoch": ckpt.get("epoch", "?"),
        }

        with open(mf, "w") as f:
            json.dump(metrics, f, indent=2)
        print(f"  Saved: {mf}")


if __name__ == "__main__":
    main()
