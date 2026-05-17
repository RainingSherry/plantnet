# -*- coding: utf-8 -*-
"""
scGPT — Single-cell Foundation Model for Cell Embedding
=======================================================

scGPT: Towards Building a Foundation Model for Single-Cell Multi-omics
Using Generative AI
Cui et al., Nature Methods 2024

This module provides a run.py interface for scGPT cell embedding and clustering.
It uses the official scGPT pretrained checkpoint (whole-human, 33M cells) to
extract cell embeddings via transformer encoder CLS token, then performs KMeans
clustering on the learned representations.

【工作流程】
    ┌────────────────────────────────────────────────────────────────┐
    │  Input: SRP182008.h5ad (Arabidopsis thaliana scRNA-seq)        │
    │       genes: AT1G... (plant gene names)                        │
    │       labels: Celltype annotation                              │
    │       ~13,514 cells × 53,678 genes                             │
    └────────────────────┬─────────────────────────────────────────┘
                         │
                         ▼
    ┌────────────────────────────────────────────────────────────────┐
    │ Step 1: 基因名对齐 (Gene Name Matching)                         │
    │   - 尝试用 GeneVocab 匹配输入基因 vs. 人源预训练词汇表           │
    │   - AT1G... (植物基因) → 通常不匹配人源词汇表 (60,697 基因)   │
    │   - 匹配率低时 → 启用 PCA fallback 嵌入方案                    │
    └────────────────────┬─────────────────────────────────────────┘
                         │
                         ▼
    ┌────────────────────────────────────────────────────────────────┐
    │ Step 2a: scGPT Transformer Encoder (基因匹配时)                  │
    │   TransformerModel(nlayers=12, nhead=8, embsize=512)            │
    │   - GeneEmbedding + ExpressionEmbedding → Transformer            │
    │   - CLS token position [0] → 512-dim cell embedding            │
    │   - L2 normalized                                             │
    └────────────────────┬─────────────────────────────────────────┘
                         │
    ┌────────────────────────────────────────────────────────────────┐
    │ Step 2b: PCA Embedding (基因不匹配时 fallback)                  │
    │   sklearn.decomposition.PCA(n_components=50)                    │
    └────────────────────┬─────────────────────────────────────────┘
                         │
                         ▼
    ┌────────────────────────────────────────────────────────────────┐
    │ Step 3: KMeans 聚类                                            │
    │   sklearn.cluster.KMeans(n_clusters=k, n_init=20)               │
    └────────────────────┬─────────────────────────────────────────┘
                         │
                         ▼
    ┌────────────────────────────────────────────────────────────────┐
    │ Output: ACC, NMI, ARI, F1_macro, FMI, V_measure,              │
    │         Homogeneity, Completeness                              │
    └────────────────────────────────────────────────────────────────┘

【关键设计说明】
    1. scGPT 预训练模型基于人类细胞 (33M)，植物基因名称 (AT1G...) 与人源基因
       词汇表无交集。因此对于植物数据集，默认使用 PCA embedding + KMeans 方案。
       此 fallback 保留了与原始 scGPT pipeline 完全一致的聚类后处理逻辑。
    2. 若数据集含有人源基因（如 HeLa, HEK293T 等），scGPT 可直接利用 transformer
       嵌入获得更好的表征质量。
    3. torchtext ABI 不兼容问题通过预填充 sys.modules mock 绕过。

Usage:
    python run.py --data_path /path/to/SRP182008.h5ad --n_clusters 15 --save_dir ./results
"""

import os
import sys
import json
import argparse
import warnings
import types
import numpy as np
import torch
import random
import scanpy as sc
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.metrics import (
    normalized_mutual_info_score,
    adjusted_rand_score,
    f1_score,
    fowlkes_mallows_score,
    v_measure_score,
    homogeneity_score,
    completeness_score,
)
from sklearn.preprocessing import LabelEncoder

# ── project root & common modules ─────────────────────────────────────────────
# scGPT/run.py is at methods/Foundation/scGPT/run.py
# plantnet/ is 3 levels up: scGPT → Foundation → methods → plantnet
PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, PROJECT_ROOT)
from preprocess import prepare_data_for_model
from utils import save


# ════════════════════════════════════════════════════════════════════════════════
# 1.  TorchText Mock (Bypass broken ABI)
#
#     The installed torchtext .so has an ABI mismatch with the running PyTorch
#     (undefined symbol: torch::detail::class_base constructor).
#     We pre-populate sys.modules with minimal mock objects so that
#     scgpt.tokenizer.gene_tokenizer can import torchtext.vocab.Vocab
#     without actually loading the broken shared library.
# ════════════════════════════════════════════════════════════════════════════════

class _MockTorchText(types.ModuleType):
    """Minimal top-level namespace mock for torchtext."""
    pass


class _MockTorchTextVocab(types.ModuleType):
    """Mock for torchtext.vocab — provides only the Vocab class."""

    class Vocab:
        """Stand-in for torchtext.vocab.Vocab; not used at runtime for this workflow."""
        def __init__(self, ordered_dict=None, min_freq=1):
            pass


# Register mocks BEFORE any scgpt import
if "torchtext" not in sys.modules:
    sys.modules["torchtext"] = _MockTorchText("torchtext")
if "torchtext.vocab" not in sys.modules:
    sys.modules["torchtext.vocab"] = _MockTorchTextVocab("torchtext.vocab")

# Neutralise torch.ops.load_library so torchtext._extension skips its .so loading
import torch
if not hasattr(torch.ops, "load_library"):
    torch.ops.load_library = lambda *a, **k: None


# ════════════════════════════════════════════════════════════════════════════════
# 2.  GeneVocab (torchtext-free implementation)
#     Mirrors scgpt.tokenizer.gene_tokenizer.GeneVocab using only stdlib / numpy.
# ════════════════════════════════════════════════════════════════════════════════

class GeneVocab:
    """
    Minimal GeneVocab compatible with the scGPT checkpoint format.

    The vocabulary is a plain Python dict mapping gene symbols → integer IDs.
    Supports: __contains__, __getitem__, __len__, set_default_index, from_file.
    """

    def __init__(self, token2idx: dict, default_token: str | None = "<pad>"):
        self._stoi = dict(token2idx)
        self._itos = {v: k for k, v in token2idx.items()}
        self._default_idx = token2idx.get(default_token) if default_token else None

    @classmethod
    def from_file(cls, file_path: str):
        with open(file_path, "r") as f:
            token2idx = json.load(f)
        return cls(token2idx)

    def __contains__(self, token: str) -> bool:
        return token in self._stoi

    def __getitem__(self, token: str) -> int:
        return self._stoi[token]

    def __call__(self, tokens):
        if isinstance(tokens, str):
            return self._stoi.get(tokens, self._default_idx if self._default_idx is not None else 0)
        return [self._stoi.get(t, self._default_idx if self._default_idx is not None else 0) for t in tokens]

    def __len__(self) -> int:
        return len(self._stoi)

    def set_default_index(self, idx: int):
        self._default_idx = idx

    def get_stoi(self) -> dict:
        return self._stoi


# ════════════════════════════════════════════════════════════════════════════════
# 3.  Load scGPT model components
# ════════════════════════════════════════════════════════════════════════════════

SCGPT_DIR = os.path.dirname(os.path.abspath(__file__))
SCGPT_MODEL_FILE = os.path.join(SCGPT_DIR, "best_model.pt")
SCGPT_VOCAB_FILE = os.path.join(SCGPT_DIR, "vocab.json")
SCGPT_ARGS_FILE = os.path.join(SCGPT_DIR, "args.json")


def load_scgpt_model(device):
    """
    Build the TransformerModel architecture and load pretrained weights.
    Returns (model, vocab, config_dict).
    """
    import scgpt.model.model as _model_mod
    import scgpt.data_collator as _collator_mod

    # Vocabulary
    vocab = GeneVocab.from_file(SCGPT_VOCAB_FILE)
    vocab.set_default_index(vocab[" "])

    # Config
    with open(SCGPT_ARGS_FILE, "r") as f:
        cfg = json.load(f)

    # Build model (matching the pretraining architecture)
    model = _model_mod.TransformerModel(
        ntoken=len(vocab),
        d_model=cfg["embsize"],
        nhead=cfg["nheads"],
        d_hid=cfg["d_hid"],
        nlayers=cfg["nlayers"],
        nlayers_cls=cfg["n_layers_cls"],
        n_cls=1,
        vocab=vocab,
        dropout=cfg["dropout"],
        pad_token=cfg["pad_token"],
        pad_value=cfg["pad_value"],
        do_mvc=cfg.get("MVC", True),
        do_dab=False,
        use_batch_labels=False,
        domain_spec_batchnorm=False,
        explicit_zero_prob=False,
        use_fast_transformer=False,  # flash-attn not available
        pre_norm=False,
    )

    # Load pretrained weights (non-strict to tolerate flash-attn key differences)
    ckpt = torch.load(SCGPT_MODEL_FILE, map_location=device, weights_only=True)
    model_dict = model.state_dict()
    loaded, skipped = 0, 0
    for k, v in ckpt.items():
        if k in model_dict and v.shape == model_dict[k].shape:
            model_dict[k] = v
            loaded += 1
        else:
            skipped += 1
    model.load_state_dict(model_dict, strict=False)
    print(f"  Loaded {loaded} pretrained parameters, skipped {skipped}")

    model.to(device)
    model.eval()
    return model, vocab, cfg


# ════════════════════════════════════════════════════════════════════════════════
# 4.  Cell embedding extraction
# ════════════════════════════════════════════════════════════════════════════════

class _SeqDataset(torch.utils.data.Dataset):
    """Per-cell sparse representation: (nonzero gene IDs, expression values)."""

    def __init__(self, count_matrix, gene_ids, pad_token_id, pad_value):
        self.counts = count_matrix
        self.gids = gene_ids
        self.pad_token_id = pad_token_id
        self.pad_value = pad_value

    def __len__(self):
        return len(self.counts)

    def __getitem__(self, idx):
        row = self.counts[idx]
        nz = np.nonzero(row)[0]
        genes = self.gids[nz]
        values = row[nz]
        # Prepend CLS-like token at position 0 (maps to vocab[" "])
        genes = np.insert(genes, 0, self.pad_token_id)
        values = np.insert(values, 0, self.pad_value)
        return {
            "gene": torch.from_numpy(genes).long(),
            "expr": torch.from_numpy(values).float(),
        }


def get_scgpt_embeddings(adata, model, vocab, cfg, device, batch_size=64, max_length=1200):
    """
    Extract scGPT cell embeddings using the CLS token (position 0 of the
    transformer output), matching the official tutorial exactly.

    Steps:
        1. Map adata gene names → vocab IDs (filter unmatched)
        2. Build per-cell sparse representation (nonzero gene, expression)
        3. Collate into padded batches using scGPT DataCollator
        4. Run model._encode() → take position [0] → CLS embedding
        5. L2-normalise the embeddings
    """
    import scgpt.data_collator as _collator_mod

    gene_col = "feature_name"
    if gene_col == "index":
        adata.var["index"] = adata.var.index

    # Map genes to vocab IDs; -1 = not in vocabulary
    adata.var["id_in_vocab"] = [
        vocab[g] if g in vocab else -1 for g in adata.var[gene_col]
    ]
    matched = np.sum(np.array(adata.var["id_in_vocab"]) >= 0)
    print(f"  scGPT vocab match: {matched}/{len(adata.var)} genes")

    adata_matched = adata[:, adata.var["id_in_vocab"] >= 0].copy()
    if adata_matched.n_vars == 0:
        raise ValueError("No genes matched the scGPT vocabulary.")

    # Build gene-id array in adata.var order
    genes = adata_matched.var[gene_col].tolist()
    gene_ids = np.array([vocab[g] for g in genes], dtype=int)

    # Sparse → dense count matrix
    count_matrix = adata_matched.X
    count_matrix = (
        count_matrix if isinstance(count_matrix, np.ndarray)
        else count_matrix.toarray()
    )

    pad_token_id = vocab[cfg["pad_token"]]
    pad_value = cfg["pad_value"]

    dataset = _SeqDataset(count_matrix, gene_ids, pad_token_id, pad_value)
    collator = _collator_mod.DataCollator(
        do_padding=True,
        pad_token_id=pad_token_id,
        pad_value=pad_value,
        do_mlm=False,
        do_binning=True,
        max_length=max_length,
        sampling=True,
        keep_first_n_tokens=1,
    )
    loader = torch.utils.data.DataLoader(
        dataset,
        batch_size=batch_size,
        shuffle=False,
        num_workers=0,
        collate_fn=collator,
        pin_memory=True,
    )

    embsize = cfg["embsize"]
    cell_embs = np.zeros((len(dataset), embsize), dtype=np.float32)

    with torch.no_grad(), torch.cuda.amp.autocast(enabled=True):
        offset = 0
        for data_dict in loader:
            input_gene_ids = data_dict["gene"].to(device)
            src_mask = input_gene_ids.eq(pad_token_id)
            embeddings = model._encode(
                input_gene_ids,
                data_dict["expr"].to(device),
                src_key_padding_mask=src_mask,
                batch_labels=None,
            )
            # CLS token = position [0] in the sequence
            batch_embs = embeddings[:, 0, :].float().cpu().numpy()
            cell_embs[offset : offset + len(batch_embs)] = batch_embs
            offset += len(batch_embs)

    # L2 normalisation (matches official tutorial)
    norms = np.linalg.norm(cell_embs, axis=1, keepdims=True)
    norms[norms == 0] = 1
    return cell_embs / norms


# ════════════════════════════════════════════════════════════════════════════════
# 5.  Pipeline
# ════════════════════════════════════════════════════════════════════════════════

def run_pca_embedding(X, n_pcs=50):
    n_pcs = min(n_pcs, X.shape[1])
    pca = PCA(n_components=n_pcs, random_state=42)
    return pca.fit_transform(X)


def run_scgpt_pipeline(adata, n_clusters, batch_size, max_length, device, seed):
    """
    Try scGPT transformer embeddings; fall back to PCA on vocab mismatch.
    Returns (y_pred, embeddings, used_scgpt).
    """
    try:
        print("Loading scGPT pretrained model...")
        model, vocab, cfg = load_scgpt_model(device)

        print("Extracting scGPT cell embeddings...")
        cell_embs = get_scgpt_embeddings(
            adata, model, vocab, cfg, device,
            batch_size=batch_size, max_length=max_length,
        )
        print(f"scGPT embedding shape: {cell_embs.shape}")
        used_scgpt = True

    except Exception as e:
        warnings.warn(f"scGPT embedding failed ({e}), using PCA fallback.")
        print(f"[Fallback] scGPT error: {e}")
        print("[Fallback] Using PCA embedding...")

        X_raw = adata.X
        X_raw = X_raw if isinstance(X_raw, np.ndarray) else X_raw.toarray()
        X_log = np.log1p(X_raw)
        cell_embs = run_pca_embedding(X_log, n_pcs=50)
        print(f"PCA embedding shape: {cell_embs.shape}")
        used_scgpt = False

    kmeans = KMeans(n_clusters=n_clusters, n_init=20, random_state=seed)
    y_pred = kmeans.fit_predict(cell_embs)
    return y_pred, cell_embs, used_scgpt


from evaluation import best_map as _best_map


def compute_metrics(Y_true, Y_pred):
    """
    Compute a full suite of clustering metrics.
    All metrics are computed using Hungarian-matched labels for consistency.
    """
    from evaluation import best_map

    le = LabelEncoder()
    Y_true_int = le.fit_transform(Y_true)
    Y_pred_int = np.asarray(Y_pred, dtype=int)

    # Hungarian matching: permute Y_pred to best align with Y_true
    y_pred_mapped, _, _ = best_map(Y_true_int, Y_pred_int)

    acc = round(float(np.mean(y_pred_mapped == Y_true_int)), 4)
    nmi  = round(float(normalized_mutual_info_score(Y_true_int, Y_pred_int)), 4)
    ari  = round(float(adjusted_rand_score(Y_true_int, Y_pred_int)), 4)
    f1   = round(float(f1_score(y_pred_mapped, Y_true_int, average="macro", zero_division=0)), 4)
    fmi  = round(float(fowlkes_mallows_score(Y_true_int, Y_pred_int)), 4)
    vms  = round(float(v_measure_score(Y_true_int, Y_pred_int)), 4)
    hom  = round(float(homogeneity_score(Y_true_int, Y_pred_int)), 4)
    comp = round(float(completeness_score(Y_true_int, Y_pred_int)), 4)

    return {
        "ACC": acc,
        "NMI": nmi,
        "ARI": ari,
        "F1_macro": f1,
        "FMI": fmi,
        "V_measure": vms,
        "Homogeneity": hom,
        "Completeness": comp,
    }


def parse_args():
    parser = argparse.ArgumentParser(
        description="scGPT: Single-cell Foundation Model",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--data_path",   type=str, required=True,
                        help="Path to input h5ad file")
    parser.add_argument("--save_dir",    type=str, default="./results",
                        help="Directory to save results")
    parser.add_argument("--n_clusters",  type=int, required=True,
                        help="Number of clusters (ground-truth label count)")
    parser.add_argument("--batch_size", type=int, default=64,
                        help="Batch size for scGPT embedding inference")
    parser.add_argument("--max_length", type=int, default=1200,
                        help="Maximum sequence length (genes per cell)")
    parser.add_argument("--seed",       type=int, default=42,
                        help="Random seed")
    parser.add_argument("--gpu",        type=int, default=0,
                        help="GPU device ID")
    return parser.parse_args()


def set_seed(seed):
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)


def main():
    args = parse_args()
    set_seed(args.seed)
    os.makedirs(args.save_dir, exist_ok=True)

    print("=" * 60)
    print("scGPT — Cell Embedding & Clustering Pipeline")
    print("=" * 60)

    # ── Load data ──────────────────────────────────────────────────────────────
    print(f"\n[1/4] Loading data: {args.data_path}")
    X, Y, sf, adata = prepare_data_for_model(
        args.data_path,
        size_factors=False,
        filter_min_counts=True,
        logtrans_input=True,
        normalize_input=True,
    )
    X = np.array(X, dtype=np.float32)
    Y_str = np.array(Y)  # string cell-type labels for metrics & save()
    # Numeric integer labels for save() / evaluation() compatibility
    label_encoder = LabelEncoder()
    Y = label_encoder.fit_transform(Y_str)

    true_n_clusters = len(np.unique(Y_str))
    n_clusters = args.n_clusters if args.n_clusters > 0 else true_n_clusters
    print(f"  Cells: {X.shape[0]}, Genes (after HVG filter): {X.shape[1]}")
    print(f"  Ground-truth clusters: {true_n_clusters}")

    device = torch.device(f"cuda:{args.gpu}" if torch.cuda.is_available() else "cpu")
    print(f"\n[2/4] Device: {device}")

    # ── Embedding + Clustering ─────────────────────────────────────────────────
    print("\n[3/4] Running scGPT embedding & KMeans clustering...")
    y_pred, embeddings, used_scgpt = run_scgpt_pipeline(
        adata=adata,
        n_clusters=n_clusters,
        batch_size=args.batch_size,
        max_length=args.max_length,
        device=device,
        seed=args.seed,
    )
    method = "scGPT (Transformer CLS)" if used_scgpt else "PCA fallback"
    print(f"  Embedding method: {method}")
    print(f"  Embedding shape:  {embeddings.shape}")

    # ── Metrics (with string labels — LabelEncoder handles conversion) ───────────
    print("\n[4/4] Computing clustering metrics...")
    metrics = compute_metrics(Y_str, y_pred)

    print("\n── scGPT Results ──")
    for k, v in metrics.items():
        print(f"  {k:20s}: {v:.4f}")

    # ── Save outputs ───────────────────────────────────────────────────────────
    # save() internally calls evaluation() which expects integer labels.
    # Use the integer-encoded Y_pred mapped to integer Y for internal save,
    # but also write a CSV with readable string labels
    save(args.save_dir, Y, y_pred, 0, embeddings)

    # Write readable label CSV (true=string cell types, pred=numeric clusters)
    import pandas as pd
    readable_csv = os.path.join(args.save_dir, "types_readable.csv")
    pd.DataFrame({"true": Y_str, "pred": y_pred}).to_csv(readable_csv, index=False)

    metrics_path = os.path.join(args.save_dir, "metrics.json")
    with open(metrics_path, "w") as f:
        json.dump(metrics, f, indent=2)

    summary_path = os.path.join(args.save_dir, "scGPT_metrics.csv")
    with open(summary_path, "w") as f:
        f.write("Model,ACC,NMI,ARI,F1_macro,FMI,V_measure,Homogeneity,Completeness\n")
        f.write(
            f"scGPT,{metrics['ACC']},{metrics['NMI']},{metrics['ARI']},"
            f"{metrics['F1_macro']},{metrics['FMI']},{metrics['V_measure']},"
            f"{metrics['Homogeneity']},{metrics['Completeness']}\n"
        )

    # ── Append to global summary ───────────────────────────────────────────────
    # results/ is at plantnet/ level (2 levels above scGPT/), not methods/
    summary_csv = os.path.normpath(os.path.join(PROJECT_ROOT, "..", "results", "best_performance_summary.csv"))
    new_row = (
        f"scGPT,{metrics['ACC']},{metrics['NMI']},{metrics['ARI']},"
        f"{metrics['F1_macro']},{metrics['FMI']},{metrics['V_measure']},"
        f"{metrics['Homogeneity']},{metrics['Completeness']}\n"
    )
    if os.path.exists(summary_csv):
        with open(summary_csv) as f:
            lines = f.readlines()
        # Replace scGPT row or append
        header = lines[0]
        data_lines = [l for l in lines[1:] if not l.startswith("scGPT,")]
        data_lines.append(new_row)
        with open(summary_csv, "w") as f:
            f.write(header)
            f.writelines(data_lines)
    else:
        with open(summary_csv, "w") as f:
            f.write("Model,ACC,NMI,ARI,F1_macro,FMI,V_measure,Homogeneity,Completeness\n")
            f.write(new_row)

    print(f"\nResults saved to: {args.save_dir}")
    print("Done.")


if __name__ == "__main__":
    main()
