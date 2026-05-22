"""
train/train_plantdiffcluster.py
==================================
PlantDiffCluster 训练脚本。

支持三种图结构消融实验：
  1. HVG Co-expression Graph  — 高变基因共表达图（推荐）
  2. Marker Graph             — Marker-boosted 图
  3. Random Graph             — 随机图（负对照）

训练流程：
  1. 加载 SRP182008 数据集
  2. 构建基因图（可配置）
  3. 初始化 PlantDiffCluster 模型
  4. 联合优化：
     - 重构损失 MSE(decoder(X), X)
     - 零值掩码损失 BCE(mask_pred, M)
     - 聚类损失 GMM变分 / DEC KL / Contrastive
  5. 周期性评估（ARI, NMI, ACC）
  6. 保存最优模型
  7. 可视化（UMAP/t-SNE）

使用方式：
  python train/train_plantdiffcluster.py --config configs/srp182008.yaml

  或直接运行：
  python train/train_plantdiffcluster.py --graph_type coexpression --n_clusters 15 --epochs 100
"""

from __future__ import annotations

import os
import sys
import json
import time
import random
import argparse
import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from pathlib import Path
from datetime import datetime
from typing import Optional, Tuple, List, Dict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scanpy as sc

# 添加项目路径
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from datasets.srp182008_dataset import SRP182008Dataset, create_dataloader
from models import PlantDiffCluster, save_checkpoint


# ---------------------------------------------------------------------------
# 评估指标（来自 scMAE/scCluBench）
# ---------------------------------------------------------------------------

def compute_clustering_metrics(y_true: np.ndarray, y_pred: np.ndarray) -> dict:
    """计算聚类指标（ACC, NMI, ARI, F1-macro 等）。"""
    from sklearn.metrics import (
        accuracy_score, f1_score,
        normalized_mutual_info_score as nmi_score,
        adjusted_rand_score as ari_score,
    )
    from scipy.optimize import linear_sum_assignment

    # Hungarian algorithm 重排
    y_true_unique = np.unique(y_true)
    y_pred_unique = np.unique(y_pred)
    n_class = len(y_true_unique)
    n_pred = len(y_pred_unique)

    # 构建匹配矩阵
    G = np.zeros((n_class, n_pred), dtype=int)
    for i, ut in enumerate(y_true_unique):
        for j, up in enumerate(y_pred_unique):
            G[i, j] = np.sum((y_true == ut) & (y_pred == up))

    # Hungarian
    A = linear_sum_assignment(-G)
    new_pred = np.zeros_like(y_pred)
    for i, up in enumerate(y_pred_unique):
        col_idx = A[1][i] if i < len(A[1]) else i % n_class
        label_idx = A[0][i] if i < len(A[0]) else i % n_class
        new_pred[y_pred == up] = y_true_unique[label_idx]

    acc = accuracy_score(y_true, new_pred)
    f1 = f1_score(y_true, new_pred, average="macro")
    nmi = nmi_score(y_true, y_pred, average_method="arithmetic")
    ari = ari_score(y_true, y_pred)

    return {"acc": acc, "f1_macro": f1, "nmi": nmi, "ari": ari}


def plot_umap(embeddings: np.ndarray, labels: np.ndarray, save_path: str, title: str = "UMAP"):
    """绘制 UMAP 可视化。"""
    try:
        adata = sc.AnnData(embeddings)
        sc.pp.neighbors(adata, n_neighbors=15)
        sc.tl.umap(adata)
        sc.pl.umap(adata, color=np.arange(len(np.unique(labels))), save=save_path, show=False)
    except Exception as e:
        print(f"[UMAP] Failed: {e}")


# ---------------------------------------------------------------------------
# 训练器
# ---------------------------------------------------------------------------

class Trainer:
    """PlantDiffCluster 训练器。"""

    def __init__(
        self,
        model: nn.Module,
        train_loader,
        val_loader,
        optimizer,
        scheduler,
        device: str,
        save_dir: str,
        config: dict,
    ):
        self.model = model
        self.train_loader = train_loader
        self.val_loader = val_loader
        self.optimizer = optimizer
        self.scheduler = scheduler
        self.device = device
        self.save_dir = Path(save_dir)
        self.config = config
        self.best_metrics = {"nmi": 0.0, "ari": 0.0, "acc": 0.0}
        self.best_epoch = 0
        self.loss_history = []

        self.save_dir.mkdir(parents=True, exist_ok=True)

        # 保存配置
        with open(self.save_dir / "config.json", "w") as f:
            json.dump(config, f, indent=2, default=str)

    def train_epoch(self, epoch: int) -> dict:
        """训练一个 epoch。"""
        self.model.train()
        total_loss = 0.0
        total_mse = 0.0
        total_mask = 0.0
        total_cluster = 0.0
        n_batches = 0

        for batch in self.train_loader:
            X = batch["X"].to(self.device)
            labels = batch["label"].to(self.device)
            support_idx = batch["support_idx"].to(self.device)
            support_mask = batch["support_mask"].to(self.device)
            support_weight = batch["support_weight"].to(self.device)

            # 采样扩散时间步
            t = torch.randint(
                0, self.model.config.get("num_timesteps", 500), (X.size(0),), device=self.device
            )

            self.optimizer.zero_grad()

            # 前向传播
            output = self.model(
                X=X,
                cell_type=labels,
                support_weight=support_weight,
                support_mask=support_mask,
                support_idx=support_idx,
                t=t,
            )

            losses = output["losses"]
            loss = losses.get("loss", losses.get("mse", torch.tensor(0.0, device=self.device)))

            # 反向传播
            loss.backward()
            torch.nn.utils.clip_grad_norm_(self.model.parameters(), max_norm=1.0)
            self.optimizer.step()

            total_loss += loss.item()
            total_mse += losses.get("mse", torch.tensor(0.0)).item()
            total_mask += losses.get("mask", torch.tensor(0.0)).item()
            total_cluster += losses.get("cluster", torch.tensor(0.0)).item()
            n_batches += 1

        avg_losses = {
            "loss": total_loss / n_batches,
            "mse": total_mse / n_batches,
            "mask": total_mask / n_batches,
            "cluster": total_cluster / n_batches,
        }
        self.loss_history.append(avg_losses)

        return avg_losses

    @torch.no_grad()
    def evaluate(self, epoch: int) -> dict:
        """在验证集上评估。"""
        self.model.eval()

        all_embeddings = []
        all_labels = []

        for batch in self.val_loader:
            X = batch["X"].to(self.device)
            labels = batch["label"].to(self.device)
            support_idx = batch["support_idx"]
            support_mask = batch["support_mask"]
            support_weight = batch["support_weight"]

            output = self.model(
                X=X,
                cell_type=labels,
                support_weight=support_weight,
                support_mask=support_mask,
                support_idx=support_idx,
                t=None,  # 推理模式
            )

            emb = output["cell_z"].cpu().numpy()
            all_embeddings.append(emb)
            all_labels.append(labels.cpu().numpy())

        embeddings = np.vstack(all_embeddings)
        labels = np.concatenate(all_labels)

        # K-Means 聚类
        from sklearn.cluster import KMeans
        kmeans = KMeans(
            n_clusters=self.model.n_clusters,
            n_init=20,
            random_state=42,
        )
        pred_labels = kmeans.fit_predict(embeddings)

        metrics = compute_clustering_metrics(labels, pred_labels)

        # UMAP 可视化（每 10 个 epoch）
        if epoch % 10 == 0:
            try:
                save_path = str(self.save_dir / f"umap_epoch{epoch}.pdf")
                plot_umap(embeddings, labels, save_path, title=f"Epoch {epoch}")
            except Exception as e:
                print(f"[UMAP] Failed: {e}")

        return metrics, embeddings, pred_labels

    def run(self, epochs: int, eval_interval: int = 5):
        """运行完整训练流程。"""
        print(f"\n{'='*60}")
        print(f"Training PlantDiffCluster for {epochs} epochs")
        print(f"Device: {self.device}")
        print(f"Save dir: {self.save_dir}")
        print(f"{'='*60}\n")

        for epoch in range(1, epochs + 1):
            start_time = time.time()

            # 训练
            train_losses = self.train_epoch(epoch)
            epoch_time = time.time() - start_time

            # 日志
            lr = self.optimizer.param_groups[0]["lr"]
            print(
                f"[Epoch {epoch}/{epochs}] "
                f"loss={train_losses['loss']:.4f} "
                f"mse={train_losses['mse']:.4f} "
                f"mask={train_losses['mask']:.4f} "
                f"cluster={train_losses['cluster']:.4f} "
                f"lr={lr:.6f} "
                f"time={epoch_time:.1f}s"
            )

            # 学习率调度
            if self.scheduler is not None:
                self.scheduler.step()

            # 评估
            if epoch % eval_interval == 0 or epoch == epochs:
                metrics, embeddings, pred_labels = self.evaluate(epoch)
                print(
                    f"[Eval] "
                    f"ACC={metrics['acc']:.4f} "
                    f"NMI={metrics['nmi']:.4f} "
                    f"ARI={metrics['ari']:.4f} "
                    f"F1={metrics['f1_macro']:.4f}"
                )

                # 保存最优模型
                if metrics["nmi"] > self.best_metrics["nmi"]:
                    self.best_metrics = metrics
                    self.best_epoch = epoch
                    save_path = self.save_dir / "best_model.pt"
                    save_checkpoint(
                        self.model, self.optimizer, epoch, str(save_path),
                        best_metrics=metrics,
                        embeddings=embeddings,
                        pred_labels=pred_labels,
                    )
                    print(f"  → Best model saved! (NMI={metrics['nmi']:.4f})")

            # 定期保存
            if epoch % 20 == 0:
                save_checkpoint(
                    self.model, self.optimizer, epoch,
                    str(self.save_dir / f"checkpoint_epoch{epoch}.pt")
                )

        print(f"\n{'='*60}")
        print(f"Training complete! Best epoch={self.best_epoch}, NMI={self.best_metrics['nmi']:.4f}")
        print(f"{'='*60}")

        # 保存损失历史
        with open(self.save_dir / "loss_history.json", "w") as f:
            json.dump(self.loss_history, f, indent=2)

        # 绘制损失曲线
        self._plot_losses()

    def _plot_losses(self):
        """绘制训练损失曲线。"""
        if not self.loss_history:
            return

        epochs = list(range(1, len(self.loss_history) + 1))
        loss_vals = [h["loss"] for h in self.loss_history]
        mse_vals = [h["mse"] for h in self.loss_history]
        mask_vals = [h["mask"] for h in self.loss_history]
        cluster_vals = [h["cluster"] for h in self.loss_history]

        fig, axes = plt.subplots(2, 2, figsize=(12, 8))
        axes[0, 0].plot(epochs, loss_vals, label="Total Loss")
        axes[0, 0].set_title("Total Loss")
        axes[0, 0].set_xlabel("Epoch")
        axes[0, 0].grid(True)

        axes[0, 1].plot(epochs, mse_vals, label="MSE", color="tab:blue")
        axes[0, 1].set_title("Reconstruction MSE")
        axes[0, 1].set_xlabel("Epoch")
        axes[0, 1].grid(True)

        axes[1, 0].plot(epochs, mask_vals, label="Mask Loss", color="tab:orange")
        axes[1, 0].set_title("Sparsity Mask Loss")
        axes[1, 0].set_xlabel("Epoch")
        axes[1, 0].grid(True)

        axes[1, 1].plot(epochs, cluster_vals, label="Cluster Loss", color="tab:green")
        axes[1, 1].set_title("Cluster Loss")
        axes[1, 1].set_xlabel("Epoch")
        axes[1, 1].grid(True)

        plt.tight_layout()
        plt.savefig(self.save_dir / "loss_curves.pdf", dpi=150)
        plt.close()
        print(f"[Loss curves] Saved to {self.save_dir / 'loss_curves.pdf'}")


# ---------------------------------------------------------------------------
# 入口
# ---------------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(
        description="Train PlantDiffCluster on SRP182008",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # 数据
    parser.add_argument("--data_path", type=str,
                        default="/home/luolie/biopipeline/dimension-reduction/plantnet/data/SRP182008.h5ad",
                        help="Path to SRP182008.h5ad")
    parser.add_argument("--save_dir", type=str,
                        default="./results/srp182008",
                        help="Results directory")
    parser.add_argument("--n_hvg", type=int, default=1500, help="Number of HVG")
    parser.add_argument("--graph_type", type=str, default="coexpression",
                        choices=["coexpression", "marker", "random"],
                        help="Gene graph type")
    parser.add_argument("--support_strategy", type=str, default="log1p",
                        choices=["log1p", "rank", "tfidf", "norm"],
                        help="Support set weighting strategy")
    parser.add_argument("--dropout_rate", type=float, default=0.0,
                        help="Support set dropout rate")

    # 模型
    parser.add_argument("--gene_dim", type=int, default=64, help="Gene embedding dim")
    parser.add_argument("--hidden_dim", type=int, default=256, help="GAT hidden dim")
    parser.add_argument("--embed_dim", type=int, default=128, help="Cell embedding dim")
    parser.add_argument("--time_embed_dim", type=int, default=128, help="Time embedding dim")
    parser.add_argument("--n_layers", type=int, default=2, help="Number of GAT layers")
    parser.add_argument("--pooling_strategy", type=str, default="attention",
                        choices=["attention", "mean", "weighted_sum", "topk"],
                        help="Pooling strategy")
    parser.add_argument("--pooling_topk", type=int, default=50, help="Top-K for pooling")
    parser.add_argument("--n_clusters", type=int, default=15, help="Number of clusters")
    parser.add_argument("--cluster_strategy", type=str, default="gmm",
                        choices=["gmm", "contrastive", "dec"],
                        help="Cluster strategy")

    # 扩散
    parser.add_argument("--use_diffusion", action="store_true", default=True,
                        help="Use diffusion refiner")
    parser.add_argument("--use_mask_predictor", action="store_true", default=True,
                        help="Use sparsity mask predictor")
    parser.add_argument("--num_timesteps", type=int, default=500, help="Diffusion timesteps")
    parser.add_argument("--ddim_steps", type=int, default=20, help="DDIM sampling steps")
    parser.add_argument("--beta_schedule", type=str, default="cosine",
                        choices=["linear", "cosine", "sqrt"])

    # 训练
    parser.add_argument("--epochs", type=int, default=100, help="Training epochs")
    parser.add_argument("--batch_size", type=int, default=64, help="Batch size")
    parser.add_argument("--lr", type=float, default=1e-3, help="Learning rate")
    parser.add_argument("--weight_decay", type=float, default=1e-4, help="Weight decay")
    parser.add_argument("--eval_interval", type=int, default=5, help="Eval every N epochs")
    parser.add_argument("--lambda_cluster", type=float, default=0.1, help="Cluster loss weight")
    parser.add_argument("--cell_type_num", type=int, default=15, help="Number of cell types")
    parser.add_argument("--seed", type=int, default=42, help="Random seed")
    parser.add_argument("--num_workers", type=int, default=0, help="DataLoader workers")

    # 其他
    parser.add_argument("--device", type=str, default="cuda", help="Device (cuda/cpu)")
    parser.add_argument("--resume", type=str, default=None, help="Resume from checkpoint")

    return parser.parse_args()


def set_seed(seed: int):
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)


def main():
    args = parse_args()
    set_seed(args.seed)

    device = args.device if torch.cuda.is_available() else "cpu"
    print(f"Device: {device}")

    # ---- 数据集 ----
    print(f"\nLoading dataset from {args.data_path}...")
    dataset = SRP182008Dataset(
        h5ad_path=args.data_path,
        n_hvg=args.n_hvg,
        graph_type=args.graph_type,
        support_strategy=args.support_strategy,
        dropout_rate=args.dropout_rate,
        random_seed=args.seed,
    )

    # 分割训练/验证集
    n_train = int(0.9 * len(dataset))
    n_val = len(dataset) - n_train
    train_dataset, val_dataset = torch.utils.data.random_split(
        dataset, [n_train, n_val],
        generator=torch.Generator().manual_seed(args.seed),
    )

    train_loader = create_dataloader(
        train_dataset, batch_size=args.batch_size, shuffle=True,
        drop_last=True, num_workers=args.num_workers,
    )
    val_loader = create_dataloader(
        val_dataset, batch_size=args.batch_size * 2, shuffle=False,
        drop_last=False, num_workers=args.num_workers,
    )

    print(f"Train: {len(train_dataset)} cells, Val: {len(val_dataset)} cells")

    # ---- 模型配置 ----
    model_config = {
        "gene_dim": args.gene_dim,
        "hidden_dim": args.hidden_dim,
        "embed_dim": args.embed_dim,
        "time_embed_dim": args.time_embed_dim,
        "n_layers": args.n_layers,
        "heads": [4] * args.n_layers,
        "pooling_strategy": args.pooling_strategy,
        "pooling_topk": args.pooling_topk,
        "n_clusters": args.n_clusters,
        "cluster_strategy": args.cluster_strategy,
        "use_diffusion": args.use_diffusion,
        "use_mask_predictor": args.use_mask_predictor,
        "num_timesteps": args.num_timesteps,
        "ddim_steps": args.ddim_steps,
        "beta_schedule": args.beta_schedule,
        "refiner_depth": 3,
        "refiner_hidden_dim": args.hidden_dim,
        "lambda_cluster": args.lambda_cluster,
        "cell_type_num": args.cell_type_num,
        "use_decoder": True,
        "decoder_hidden_dim": args.hidden_dim,
        "dropout": 0.1,
    }

    # ---- 模型 ----
    print(f"\nBuilding PlantDiffCluster...")
    model = PlantDiffCluster(
        n_genes=dataset.n_hvg_actual,
        gene_names=list(dataset.gene_names),
        graph_dict=dataset.graph_dict,
        config=model_config,
    ).to(device)

    # 统计参数量
    n_params = sum(p.numel() for p in model.parameters() if p.requires_grad)
    print(f"Model parameters: {n_params:,}")

    # ---- 优化器 ----
    optimizer = torch.optim.AdamW(
        model.parameters(), lr=args.lr, weight_decay=args.weight_decay,
    )
    scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(
        optimizer, T_max=args.epochs, eta_min=1e-6,
    )

    # ---- 训练 ----
    trainer = Trainer(
        model=model,
        train_loader=train_loader,
        val_loader=val_loader,
        optimizer=optimizer,
        scheduler=scheduler,
        device=device,
        save_dir=args.save_dir,
        config={
            "args": vars(args),
            "model_config": model_config,
        },
    )

    trainer.run(epochs=args.epochs, eval_interval=args.eval_interval)

    print(f"\nAll results saved to: {args.save_dir}")


if __name__ == "__main__":
    main()
