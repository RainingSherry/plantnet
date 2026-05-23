import argparse
import json
import os
import random
import sys
from pathlib import Path
from typing import List

import numpy as np
import torch
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score

if not hasattr(np, "string_"):
    np.string_ = np.bytes_
if not hasattr(np, "unicode_"):
    np.unicode_ = np.str_

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT / "methods"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from evaluation import evaluation
from utils import save
from graphs.build_cell_gene_graph import build_cell_support_bundle
from graphs.build_gene_graph import build_gene_graph_bundle, save_gene_graph_bundle
from model import GraphDiffusionModel


def set_seed(seed: int) -> None:
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed(seed)
        torch.cuda.manual_seed_all(seed)


def parse_args():
    parser = argparse.ArgumentParser(description="Plant GraphDiffusion for scRNA-seq clustering")
    parser.add_argument("--data_path", type=str, required=True)
    parser.add_argument("--save_dir", type=str, required=True)
    parser.add_argument("--n_clusters", type=int, default=0)
    parser.add_argument("--label_col", type=str, default="Celltype")
    parser.add_argument("--n_top_genes", type=int, default=2000)
    parser.add_argument("--gene_feature_dim", type=int, default=64)
    parser.add_argument("--gene_hidden_dim", type=int, default=64)
    parser.add_argument("--gene_output_dim", type=int, default=64)
    parser.add_argument("--cell_dim", type=int, default=64)
    parser.add_argument("--epochs", type=int, default=80)
    parser.add_argument("--lr", type=float, default=1e-3)
    parser.add_argument("--weight_decay", type=float, default=1e-4)
    parser.add_argument("--corr_threshold", type=float, default=0.35)
    parser.add_argument("--coexpr_top_k", type=int, default=20)
    parser.add_argument("--support_weight_mode", type=str, default="log1p_count", choices=["log1p_count", "normalized_count", "rank_weight", "tfidf"])
    parser.add_argument("--save_graph", action="store_true")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--gpu", type=int, default=0)
    parser.add_argument("--no_cuda", action="store_true")
    return parser.parse_args()


def tensors_for_support(support_list: List[np.ndarray], device: torch.device) -> List[torch.Tensor]:
    return [torch.tensor(x, dtype=torch.long, device=device) for x in support_list]


def weight_tensors_for_support(weight_list: List[np.ndarray], device: torch.device) -> List[torch.Tensor]:
    return [torch.tensor(x, dtype=torch.float32, device=device) for x in weight_list]


def main():
    args = parse_args()
    set_seed(args.seed)
    device = torch.device(f"cuda:{args.gpu}" if torch.cuda.is_available() and not args.no_cuda else "cpu")
    os.makedirs(args.save_dir, exist_ok=True)

    print("Building plant gene graph...")
    gene_graph = build_gene_graph_bundle(
        data_path=args.data_path,
        label_col=args.label_col,
        n_top_genes=args.n_top_genes,
        corr_threshold=args.corr_threshold,
        coexpr_top_k=args.coexpr_top_k,
        feature_dim=args.gene_feature_dim,
        seed=args.seed,
    )
    save_gene_graph_bundle(gene_graph, os.path.join(args.save_dir, "graph"))

    print("Building cell support graph...")
    support_bundle = build_cell_support_bundle(
        data_path=args.data_path,
        gene_graph=gene_graph,
        label_col=args.label_col,
        support_weight_mode=args.support_weight_mode,
    )

    y_true = support_bundle.labels
    n_clusters = args.n_clusters if args.n_clusters > 0 else len(np.unique(y_true))
    module_count = max(8, min(64, len(gene_graph.marker_modules) if len(gene_graph.marker_modules) > 0 else n_clusters * 2))

    model = GraphDiffusionModel(
        gene_input_dim=gene_graph.gene_features.shape[1],
        gene_hidden_dim=args.gene_hidden_dim,
        gene_output_dim=args.gene_output_dim,
        cell_dim=args.cell_dim,
        module_count=module_count,
        n_clusters=n_clusters,
        gat_layers=2,
        dropout=0.1,
    ).to(device)

    optimizer = torch.optim.AdamW(model.parameters(), lr=args.lr, weight_decay=args.weight_decay)

    gene_features = torch.tensor(gene_graph.gene_features, dtype=torch.float32, device=device)
    edge_index = torch.tensor(gene_graph.edge_index, dtype=torch.long, device=device)
    edge_weight = torch.tensor(gene_graph.edge_weight, dtype=torch.float32, device=device)
    support_indices = tensors_for_support(support_bundle.support_indices, device)
    support_weights = weight_tensors_for_support(support_bundle.support_weights, device)
    x_dense = torch.tensor(support_bundle.x_dense, dtype=torch.float32, device=device)

    best_ari = -1.0
    best_embedding = None
    best_pred = None
    best_module_activation = None

    print(f"Training on {device} with {support_bundle.x_dense.shape[0]} cells and {support_bundle.x_dense.shape[1]} genes")
    for epoch in range(args.epochs):
        model.train()
        optimizer.zero_grad()
        gene_embeddings, cell_embeddings, refined_embeddings, logits, probs, pooled_modules, module_activation = model(
            gene_features=gene_features,
            edge_index=edge_index,
            edge_weight=edge_weight,
            support_indices=support_indices,
            support_weights=support_weights,
        )
        losses = model.compute_losses(
            refined_embeddings=refined_embeddings,
            logits=logits,
            module_activation=module_activation,
            pooled_modules=pooled_modules,
            cluster_head=model.cluster_head,
        )
        losses["total"].backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), 5.0)
        optimizer.step()

        with torch.no_grad():
            embedding_np = refined_embeddings.detach().cpu().numpy()
            pred_labels = KMeans(n_clusters=n_clusters, random_state=args.seed, n_init=20).fit_predict(embedding_np)
            acc, nmi, ari, f1_macro, fmi, v_measure, hom, com, mapped_pred = evaluation(y_true, pred_labels)
        if ari > best_ari:
            best_ari = ari
            best_embedding = embedding_np.copy()
            best_pred = mapped_pred.copy()
            best_module_activation = module_activation.detach().cpu().numpy().copy()
            save(args.save_dir, y_true, mapped_pred, epoch, best_embedding)
            np.save(os.path.join(args.save_dir, "module_activation.npy"), best_module_activation)
            np.save(os.path.join(args.save_dir, "gene_embeddings.npy"), gene_embeddings.detach().cpu().numpy())
            with open(os.path.join(args.save_dir, "training_state.json"), "w") as handle:
                json.dump({
                    "epoch": epoch,
                    "acc": acc,
                    "nmi": nmi,
                    "ari": ari,
                    "f1_macro": f1_macro,
                    "v_measure": v_measure,
                    "homogeneity": hom,
                    "completeness": com,
                }, handle, indent=2)
        if (epoch + 1) % 10 == 0 or epoch == args.epochs - 1:
            sil = silhouette_score(embedding_np, pred_labels) if len(np.unique(pred_labels)) > 1 else -1.0
            print(
                f"Epoch {epoch + 1}/{args.epochs} loss={losses['total'].item():.4f} "
                f"ARI={ari:.4f} NMI={nmi:.4f} ACC={acc:.4f} SIL={sil:.4f}"
            )

    if best_embedding is not None:
        np.save(os.path.join(args.save_dir, "best_embedding.npy"), best_embedding)
        np.save(os.path.join(args.save_dir, "best_pred.npy"), best_pred)
        torch.save(model.state_dict(), os.path.join(args.save_dir, "graphdiffusion_model.pt"))
        top_gene_scores = np.linalg.norm(np.load(os.path.join(args.save_dir, "gene_embeddings.npy")), axis=1)
        top_gene_idx = np.argsort(top_gene_scores)[::-1][:200]
        top_gene_records = [
            {"gene": gene_graph.gene_names[int(idx)], "score": float(top_gene_scores[int(idx)])}
            for idx in top_gene_idx
        ]
        with open(os.path.join(args.save_dir, "top_contributing_genes.json"), "w") as handle:
            json.dump(top_gene_records, handle, indent=2)

    print(f"Finished. Best ARI={best_ari:.4f}. Results saved to {args.save_dir}")


if __name__ == "__main__":
    main()
