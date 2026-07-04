#!/usr/bin/env python3
"""Phase 2 —— 条件化 / nuisance 匹配 swap-noise 消融（骨干 = scMAE + DEC + 每维 std-floor 赢家）。

科学假设（见 ROADMAP.md Section 5 Phase 2）：
    原始 scMAE 的 swap-noise 从「全体细胞」里为某个基因抽取供体值，因此「这个值是否
    被替换过」常常能靠技术轴（文库大小 library size / 检测基因数 n_detected / 零率
    zero-rate）判断出来 —— 这是一条捷径。若把供体池限制在「同一 nuisance 分箱」内，
    模型就无法再靠技术轴识别被换的位置，只能去学更细的 gene-gene 条件结构。

五个 corruption arm（--corruption）：
    zero          : 零填充 mask（= 现有赢家 AdaptiveSwitch.random_mask，作为「复现赢家」对照）
    swap_global   : S0 全局 swap，供体从全体细胞抽（复现原始 scMAE swap-noise，swap 基线）
    swap_lib      : S1 供体只从「同一 library-size 分箱」抽
    swap_ndet     : S2 供体只从「同一 n_detected 分箱」抽
    swap_zerolib  : S3 供体从「(zero-rate x library) 联合分箱」抽

判据：若 S1/S2/S3 的 ARI 相对 S0 提升 >= 0.02 且多种子稳定，说明真正从 scMAE 内部改对了。

严格控制变量：除 corruption 外，mask_prob / DEC / std-floor / force_gate=1 /
variance_weight=0.02 等一切与赢家一致，保证干净归因。

防泄露：不用测试 ARI 选任何超参；label 仅在最终评测用一次（见 Section 4）。
"""
from __future__ import annotations

import argparse
import os
import random
import sys
import time
from pathlib import Path

for _var in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_var, "1")
os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")
os.environ.setdefault("NUMBA_CACHE_DIR", "/tmp/numba-cache")

import numpy as np
import scanpy as sc
import torch
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import normalize
from torch.utils.data import DataLoader, Dataset

ROOT = Path(__file__).resolve().parents[2]
ADAPTIVE_SWITCH_DIR = ROOT / "experimental_retired_models" / "Granularity_scMAE_experiments" / "AdaptiveSwitch_scMAE"
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(ADAPTIVE_SWITCH_DIR) not in sys.path:
    sys.path.insert(0, str(ADAPTIVE_SWITCH_DIR))
from clusterability import compute_clusterability
from loss import AdaptiveSwitchLoss, compute_gate
from model import AdaptiveSwitchScMAE
from methods.DeepLearning import scMAE_family as family
from methods.DeepLearning.scMAE_family import ensure_csr, select_count_source
from methods.shared_utils import ensure_dir, save_json


def _register_null_h5ad_reader() -> None:
    """兼容部分 h5ad 里 null 编码的字段（否则 anndata 读取会报错）。"""
    try:
        import h5py
        from anndata._io.specs.registry import _REGISTRY, IOSpec

        def _read_null(*args, **kwargs):
            return None

        for typ in (h5py.Dataset, h5py.Group):
            try:
                _REGISTRY.register_read(typ, IOSpec("null", "0.1.0"))(_read_null)
            except Exception:
                pass
    except Exception:
        pass


_register_null_h5ad_reader()


class ExprDataset(Dataset):
    """返回 (全局索引, 编码器输入 enc, 未缩放 log 表达 log, 标签)。

    - enc  : scale 后的编码器输入（corruption 施加在这个空间，与 model.random_mask 一致）
    - log  : 未缩放 log1p 表达，作为 scMAE 重构 target（与赢家一致）
    """

    def __init__(self, enc: np.ndarray, log: np.ndarray, labels: np.ndarray):
        self.enc = torch.as_tensor(enc, dtype=torch.float32)
        self.log = torch.as_tensor(log, dtype=torch.float32)
        self.labels = torch.as_tensor(labels, dtype=torch.long)

    def __len__(self) -> int:
        return int(self.enc.shape[0])

    def __getitem__(self, idx: int):
        return int(idx), self.enc[idx], self.log[idx], self.labels[idx]


def compute_nuisance(data_path: str, input_mode: str, n_cells_expected: int,
                     log_expr: np.ndarray) -> dict:
    """从 UNSCALED 原始 counts 计算 nuisance 变量（library size / n_detected / zero-rate）。

    关键：nuisance 必须从「未 scale」的原始数据算，绝不能从编码器输入（scale 后）算 ——
    否则 library-size 这条技术轴已被归一化抹掉，Phase 2 的整个假设就不成立。

    - library_size : 每个细胞的原始总 counts（全基因，未做 HVG 截断，语义最干净）
    - n_detected   : 每个细胞检测到的基因数（原始 counts > 0 的个数，全基因）
    - zero_rate    : 每个细胞在 HVG 空间里的零率（= 编码器实际看到的稀疏度）

    若数据本身没有可用的原始 counts（只有 log1p），则用 log_expr 退化计算并在 profile 标注。
    返回各 nuisance 向量 + 元信息 dict。
    """
    n_cells = int(log_expr.shape[0])
    source_desc = "raw-counts"
    inferred = "raw"
    lib = ndet = None
    try:
        adata = sc.read_h5ad(data_path)
        source_x, _gene_names, _var, source_desc, inferred = select_count_source(adata, input_mode)
        counts = ensure_csr(source_x)  # [n_cells, n_all_genes]，稀疏，不 densify
        if int(counts.shape[0]) == n_cells and inferred == "raw":
            lib = np.asarray(counts.sum(axis=1)).ravel().astype(np.float64)
            ndet = np.asarray((counts > 0).sum(axis=1)).ravel().astype(np.float64)
    except Exception as exc:  # noqa: BLE001 —— 读原始 counts 失败则安全退化
        source_desc = f"raw-read-failed:{type(exc).__name__}"
        inferred = "fallback"

    fallback = lib is None or ndet is None
    if fallback:
        # 退化路径：从 log1p HVG 表达近似（library=对数表达行和；n_detected=非零基因数）
        lib = log_expr.astype(np.float64).sum(axis=1)
        ndet = (log_expr > 0.0).sum(axis=1).astype(np.float64)
        source_desc = source_desc if source_desc.startswith("raw-read-failed") else "log1p-fallback"
        inferred = "fallback"

    zero_rate = (log_expr <= 0.0).mean(axis=1).astype(np.float64)  # HVG 空间稀疏度
    meta = {
        "nuisance_source": source_desc,
        "inferred_mode": inferred,
        "fallback_used": bool(fallback),
        "n_cells": int(n_cells),
        "library_min": float(lib.min()), "library_median": float(np.median(lib)), "library_max": float(lib.max()),
        "n_detected_min": float(ndet.min()), "n_detected_median": float(np.median(ndet)), "n_detected_max": float(ndet.max()),
        "zero_rate_min": float(zero_rate.min()), "zero_rate_median": float(np.median(zero_rate)), "zero_rate_max": float(zero_rate.max()),
    }
    return {"library": lib, "n_detected": ndet, "zero_rate": zero_rate, "meta": meta}


def quantile_bins(values: np.ndarray, n_bins: int) -> np.ndarray:
    """等频分箱：按秩切成 n_bins 个近似等大小的箱，返回每个细胞的箱标签 [0, n_bins)。

    用秩（argsort of argsort）而非 np.quantile 边界，天然对大量并列值（scRNA 常见）稳健，
    保证每个箱的样本数尽量均衡，避免某些箱空掉导致供体池退化。
    """
    n = int(values.shape[0])
    n_bins = max(1, min(int(n_bins), n))
    order = np.argsort(values, kind="stable")
    ranks = np.empty(n, dtype=np.int64)
    ranks[order] = np.arange(n, dtype=np.int64)
    # 秩 -> 箱号，等频切分
    bins = (ranks * n_bins) // n
    return np.clip(bins, 0, n_bins - 1).astype(np.int64)


def build_bin_labels(corruption: str, nuisance: dict, n_nuisance_bins: int, n_joint_bins: int) -> np.ndarray:
    """按 corruption arm 生成每个细胞的「供体池分箱标签」。

    - swap_global : 全部细胞同一个箱（供体池 = 全体）
    - swap_lib    : 按 library-size 等频分箱
    - swap_ndet   : 按 n_detected 等频分箱
    - swap_zerolib: (zero-rate 分箱) x (library 分箱) 的联合箱（各 n_joint_bins 个）
    """
    n = int(nuisance["library"].shape[0])
    if corruption in ("zero", "swap_global"):
        return np.zeros(n, dtype=np.int64)
    if corruption == "swap_lib":
        return quantile_bins(nuisance["library"], n_nuisance_bins)
    if corruption == "swap_ndet":
        return quantile_bins(nuisance["n_detected"], n_nuisance_bins)
    if corruption == "swap_zerolib":
        zb = quantile_bins(nuisance["zero_rate"], n_joint_bins)
        lb = quantile_bins(nuisance["library"], n_joint_bins)
        return (zb * int(n_joint_bins) + lb).astype(np.int64)
    raise ValueError(f"unknown corruption {corruption}")


class SwapCorruptor:
    """条件化 swap-noise 供体采样器（施加在编码器输入空间）。

    对 batch 里每个细胞、每个被选中的基因位置，独立地从「同一 nuisance 箱」里随机抽
    一个别的细胞的该基因值来替换（per-gene independent donor —— 这才是真正的 swap-noise，
    区别于 scMAE_family.apply_scmae_noise 那种「整行换同一个供体细胞」的实现）。

    与 model.random_mask 的区别只在「替换值来源」这一个变量：
      - random_mask（zero）: 被选位置填 0
      - SwapCorruptor       : 被选位置填 同箱内另一细胞的同基因值
    返回的 mask 指示矩阵语义完全一致（1=该位置被选中 corrupt），因此 loss.py 的 BCE
    mask-discriminator 和加权重构无需任何改动，保证干净归因。

    注意：供体从「全体 enc 矩阵」的同箱细胞里抽（不是仅 batch 内），避免小 batch 把供体
    池压得过小。zero 箱内自身也可能被抽中（近似恒等），靠 donor-pool 大小统计来监控。
    """

    def __init__(self, enc: np.ndarray, bin_labels: np.ndarray, seed: int):
        self.enc = np.ascontiguousarray(enc, dtype=np.float32)  # [N, G] 全体编码器输入
        self.bin_labels = np.asarray(bin_labels, dtype=np.int64)
        self.rng = np.random.default_rng(int(seed) + 4242)
        # 预建：箱号 -> 该箱所有细胞的全局索引
        self.bin_to_members: dict[int, np.ndarray] = {}
        for b in np.unique(self.bin_labels):
            self.bin_to_members[int(b)] = np.where(self.bin_labels == b)[0].astype(np.int64)

    def pool_stats(self) -> dict:
        """各箱供体池大小统计 —— 监控「箱太窄 -> swap 退化成近似恒等」这一失败模式。"""
        sizes = np.array([m.shape[0] for m in self.bin_to_members.values()], dtype=np.float64)
        n_total = float(self.bin_labels.shape[0])
        # cell-weighted 平均：随机一个细胞期望遇到的供体池大小
        cell_weighted = float((sizes * sizes).sum() / max(n_total, 1.0))
        return {
            "n_bins": int(sizes.shape[0]),
            "pool_mean": float(sizes.mean()),
            "pool_min": float(sizes.min()),
            "pool_max": float(sizes.max()),
            "pool_median": float(np.median(sizes)),
            "pool_cell_weighted_mean": cell_weighted,
            "bin_sizes": [int(s) for s in sizes.tolist()],
        }

    def corrupt(self, global_idx: np.ndarray, mask_prob: float, device):
        """对一个 batch 生成 (corrupted_enc[tensor], mask[tensor], effective_change_rate)。

        - global_idx : batch 内每个样本的全局行号（用于查箱 & 取原值）
        - mask_prob  : 与 zero-mask 完全相同的选中概率
        效果：先 bernoulli 选位置（并保证每行至少一位被选，和 random_mask 一致），
        再对每个被选位置从同箱抽供体值替换。
        """
        idx = np.asarray(global_idx, dtype=np.int64)
        x = self.enc[idx]                      # [B, G] 原始（batch 顺序）
        B, G = x.shape
        # 1) 选中矩阵（设计 mask，语义与 model.random_mask 一致）
        sel = self.rng.random((B, G)) < float(mask_prob)
        empty = ~sel.any(axis=1)
        if empty.any():
            cols = self.rng.integers(0, G, size=int(empty.sum()))
            sel[np.where(empty)[0], cols] = True
        # 2) 逐箱、逐列向量化抽供体
        corrupted = x.copy()
        for b in np.unique(self.bin_labels[idx]):
            members = self.bin_to_members[int(b)]           # 该箱全体细胞（全局索引）
            rows_local = np.where(self.bin_labels[idx] == b)[0]  # batch 内属于该箱的行
            if members.shape[0] == 0:
                continue
            sub_sel = sel[rows_local]                        # [b_rows, G]
            if not sub_sel.any():
                continue
            rr, cc = np.where(sub_sel)                       # 该子块内被选位置
            # 为每个 (行, 列) 独立抽一个供体细胞（在池内）。直接按全局行号索引，
            # 避免每 batch 物化 donor_bank=self.enc[members]（swap_global 时那是 44808x1000
            # ≈179MB 的拷贝，是 CPU 瓶颈）。members[donor_pick] 给出全局行号，只 gather 需要的位置。
            donor_pick = self.rng.integers(0, members.shape[0], size=rr.shape[0])
            corrupted[rows_local[rr], cc] = self.enc[members[donor_pick], cc]
        mask = sel.astype(np.float32)
        # 3) 实际改动比例（诊断：与 mask_prob 比较，确认 swap 真在换、没退化成恒等）
        changed = (corrupted != x)
        eff = float(changed.mean())
        corrupted_t = torch.as_tensor(corrupted, dtype=torch.float32, device=device)
        mask_t = torch.as_tensor(mask, dtype=torch.float32, device=device)
        return corrupted_t, mask_t, eff


def parse_args():
    p = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--data_path", required=True)
    p.add_argument("--save_dir", required=True)
    p.add_argument("--dataset_name", default=None)
    p.add_argument("--label_key", default="auto")
    p.add_argument("--input_mode", default="auto", choices=["auto", "raw", "log1p"])
    p.add_argument("--n_top_genes", type=int, default=1000)
    p.add_argument("--target_sum", type=float, default=10000.0)
    p.add_argument("--scale_input", type=family.str2bool, default=True)
    p.add_argument("--n_clusters", type=int, required=True)
    p.add_argument("--epochs", type=int, default=80)
    p.add_argument("--batch_size", type=int, default=256)
    p.add_argument("--lr", type=float, default=1e-3)
    p.add_argument("--weight_decay", type=float, default=1e-5)
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--gpu", type=int, default=1)
    p.add_argument("--no_cuda", action="store_true")
    p.add_argument("--skip_eval", action="store_true")
    p.add_argument("--hidden_size", type=int, default=128)
    p.add_argument("--dropout", type=float, default=0.05)
    p.add_argument("--mask_prob", type=float, default=0.4)
    p.add_argument("--masked_data_weight", type=float, default=0.75)
    p.add_argument("--mask_weight", type=float, default=0.65)
    p.add_argument("--cluster_weight", type=float, default=0.35)
    p.add_argument("--consistency_weight", type=float, default=0.05)
    p.add_argument("--variance_weight", type=float, default=0.02)
    p.add_argument("--var_mode", default="hinge", choices=["hinge", "cov", "koleo", "both"])
    p.add_argument("--entropy_weight", type=float, default=0.10)
    p.add_argument("--confidence_threshold", type=float, default=0.35)
    p.add_argument("--warmup_epochs", type=int, default=20)
    p.add_argument("--target_update_interval", type=int, default=5)
    p.add_argument("--neighbor_k", type=int, default=15)
    p.add_argument("--knn_pca_dim", type=int, default=50)
    p.add_argument("--gate_kappa", type=float, default=0.15)
    p.add_argument("--force_gate", type=float, default=1.0)
    p.add_argument("--gate_ema", type=float, default=0.5)
    # Phase 2 conditional-shuffle 专属旋钮
    p.add_argument("--corruption", default="zero",
                   choices=["zero", "swap_global", "swap_lib", "swap_ndet", "swap_zerolib"],
                   help="corruption 方式；zero=复现赢家，swap_*=Phase 2 的四个 arm")
    p.add_argument("--n_nuisance_bins", type=int, default=10, help="S1/S2 单轴等频箱数")
    p.add_argument("--n_joint_bins", type=int, default=5, help="S3 联合分箱每轴箱数（总箱=平方）")
    p.add_argument("--smoke", action="store_true")
    return p.parse_args()


def set_seed(seed: int) -> None:
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False


def get_device(gpu: int, no_cuda: bool) -> torch.device:
    if no_cuda or not torch.cuda.is_available():
        return torch.device("cpu")
    if gpu in {0, 7}:
        raise ValueError("Physical GPU 0 and GPU 7 are intentionally avoided.")
    return torch.device(f"cuda:{gpu}")


def build_neighbor_indices(data: np.ndarray, k: int, pca_dim: int, seed: int):
    n_cells, n_genes = data.shape
    k = max(1, min(int(k), n_cells - 1))
    dim = max(2, min(int(pca_dim), n_cells - 1, n_genes - 1))
    emb = PCA(n_components=dim, random_state=seed).fit_transform(data.astype(np.float64))
    emb = normalize(emb, norm="l2", axis=1)
    nn_ = NearestNeighbors(n_neighbors=k + 1, metric="cosine").fit(emb)
    _, neighbors = nn_.kneighbors(emb, return_distance=True)
    return neighbors[:, 1:].astype(np.int64)


@torch.no_grad()
def extract_all(model, loader, device):
    """干净（未 corrupt）前向，导出 KMeans 聚类用的 embedding 与 soft-assignment q。"""
    model.eval()
    emb, q_all, labels = [], [], []
    for _, x, _, y in loader:
        out = model(x.to(device))
        emb.append(out["latent"].detach().cpu().numpy())
        q_all.append(out["cluster_q"].detach().cpu().numpy())
        labels.append(y.numpy())
    emb = np.nan_to_num(np.concatenate(emb).astype(np.float32))
    q_all = np.concatenate(q_all).astype(np.float32)
    labels = np.concatenate(labels).astype(np.int64)
    return emb, q_all, labels


def effective_dimensionality(std: np.ndarray) -> dict:
    var = np.square(std.astype(np.float64))
    pr = float((var.sum() ** 2) / max(float(np.square(var).sum()), 1e-12))
    return {
        "std_min": float(std.min()),
        "std_median": float(np.median(std)),
        "std_max": float(std.max()),
        "effective_dim_pr": pr,
        "dims_std_gt_0p1": int((std > 0.1).sum()),
        "dims_std_gt_1p0": int((std > 1.0).sum()),
    }


def cluster_aligned_eff_dim(emb: np.ndarray, labels: np.ndarray) -> dict:
    """between-class scatter 特征谱的 participation ratio：多少维真正携带细胞型判别信号。"""
    emb = emb.astype(np.float64)
    grand = emb.mean(axis=0, keepdims=True)
    classes = np.unique(labels)
    d = emb.shape[1]
    Sb = np.zeros((d, d), dtype=np.float64)
    for c in classes:
        m = labels == c
        n_c = int(m.sum())
        if n_c == 0:
            continue
        diff = emb[m].mean(axis=0, keepdims=True) - grand
        Sb += n_c * (diff.T @ diff)
    ev = np.linalg.eigvalsh(Sb).clip(min=0.0)
    pr = float((ev.sum() ** 2) / max(float(np.square(ev).sum()), 1e-12))
    return {"cluster_aligned_eff_dim": pr, "between_class_scatter_trace": float(np.trace(Sb))}


def main() -> int:
    args = parse_args()
    if args.smoke:
        args.epochs = min(args.epochs, 3)
        args.warmup_epochs = min(args.warmup_epochs, 1)
    set_seed(args.seed)
    save_dir = Path(ensure_dir(args.save_dir))
    save_json(vars(args), str(save_dir / "args.json"))
    device = get_device(args.gpu, args.no_cuda)
    dataset_name = args.dataset_name or Path(args.data_path).stem

    # target（未 scale log1p，重构 target）与 encoder（scale 后，corruption 空间）分开加载
    target_bundle = family.load_scmae_dataset(
        args.data_path, args.input_mode, args.n_top_genes, args.target_sum, False, args.label_key, args.seed
    )
    if args.scale_input:
        encoder_bundle = family.load_scmae_dataset(
            args.data_path, args.input_mode, args.n_top_genes, args.target_sum, True, args.label_key, args.seed
        )
        encoder_data = np.asarray(encoder_bundle.data, dtype=np.float32)
    else:
        encoder_data = np.asarray(target_bundle.data, dtype=np.float32)
    log_expr = np.asarray(target_bundle.data, dtype=np.float32)
    labels = np.asarray(target_bundle.labels, dtype=np.int64)
    n_clusters = int(args.n_clusters if args.n_clusters > 0 else len(np.unique(labels)))
    save_json(target_bundle.profile, str(save_dir / "dataset_profile.json"))

    # nuisance 从 UNSCALED 原始 counts 算；再据此建供体池分箱
    nuisance = compute_nuisance(args.data_path, args.input_mode, encoder_data.shape[0], log_expr)
    bin_labels = build_bin_labels(args.corruption, nuisance, args.n_nuisance_bins, args.n_joint_bins)
    corruptor = None
    pool_stats = {"n_bins": 1, "pool_mean": float(encoder_data.shape[0]),
                  "pool_min": float(encoder_data.shape[0]), "pool_max": float(encoder_data.shape[0]),
                  "pool_median": float(encoder_data.shape[0]),
                  "pool_cell_weighted_mean": float(encoder_data.shape[0]), "bin_sizes": [int(encoder_data.shape[0])]}
    if args.corruption != "zero":
        corruptor = SwapCorruptor(encoder_data, bin_labels, args.seed)
        pool_stats = corruptor.pool_stats()
    neighbor_profile = {"corruption": args.corruption, "nuisance_meta": nuisance["meta"],
                        "n_nuisance_bins": int(args.n_nuisance_bins), "n_joint_bins": int(args.n_joint_bins),
                        "donor_pool": pool_stats}
    save_json(neighbor_profile, str(save_dir / "neighbor_profile.json"))

    nb_indices = build_neighbor_indices(encoder_data, args.neighbor_k, args.knn_pca_dim, args.seed)

    dataset = ExprDataset(encoder_data, log_expr, labels)
    generator = torch.Generator().manual_seed(args.seed)
    train_loader = DataLoader(dataset, batch_size=args.batch_size, shuffle=True, drop_last=True, generator=generator)
    full_loader = DataLoader(dataset, batch_size=max(args.batch_size * 4, 512), shuffle=False, drop_last=False)

    model = AdaptiveSwitchScMAE(encoder_data.shape[1], n_clusters, args.hidden_size, args.dropout).to(device)
    criterion = AdaptiveSwitchLoss(
        args.masked_data_weight, args.mask_weight, args.cluster_weight,
        args.consistency_weight, args.variance_weight, args.entropy_weight,
        args.confidence_threshold, args.var_mode,
    )
    optimizer = torch.optim.Adam(model.parameters(), lr=args.lr, weight_decay=args.weight_decay)

    p_targets = None
    clusterab = np.ones(encoder_data.shape[0], dtype=np.float32)
    gate = float(args.force_gate)
    kl_ref = 0.0
    centers_initialized = False
    eff_change_accum: list[float] = []  # 累计实际改动率，验证 swap 真在换
    history = {k: [] for k in ["loss", "scmae_loss", "sharp_loss", "variance_loss", "gate", "kl_ref", "eff_change", "designed_mask"]}
    start = time.time()
    print(
        f"Device={device} dataset={dataset_name} cells={encoder_data.shape[0]} clusters={n_clusters} "
        f"corruption={args.corruption} donor_pool_mean={pool_stats['pool_mean']:.1f} "
        f"n_bins={pool_stats['n_bins']} mask_prob={args.mask_prob} varw={args.variance_weight} "
        f"nuisance_src={nuisance['meta']['nuisance_source']}"
    )

    for epoch in range(1, max(1, args.epochs) + 1):
        if epoch > args.warmup_epochs and (
            (epoch - args.warmup_epochs - 1) % max(1, args.target_update_interval) == 0 or p_targets is None
        ):
            emb, q_full, _ = extract_all(model, full_loader, device)
            if not centers_initialized:
                km = KMeans(n_clusters=n_clusters, n_init=20, random_state=args.seed).fit(emb)
                model.initialize_centers(torch.as_tensor(km.cluster_centers_, dtype=torch.float32, device=device))
                emb, q_full, _ = extract_all(model, full_loader, device)
                centers_initialized = True
            sharp_p = AdaptiveSwitchScMAE.sharpen(torch.as_tensor(q_full)).numpy().astype(np.float32)
            p_targets = sharp_p
            clusterab, _ = compute_clusterability(emb, q_full, nb_indices, k=args.neighbor_k)
            g_new, kl_ref = compute_gate(sharp_p, q_full, args.gate_kappa)
            gate = args.gate_ema * gate + (1.0 - args.gate_ema) * g_new if args.force_gate < 0 else float(args.force_gate)

        cluster_scale = 0.0 if epoch <= args.warmup_epochs else min(1.0, (epoch - args.warmup_epochs) / max(1, args.warmup_epochs))
        model.train()
        sums = {k: 0.0 for k in ["loss", "scmae_loss", "sharp_loss", "variance_loss", "eff_change", "designed_mask"]}
        batches = 0
        for idx, x_cpu, log_cpu, _ in train_loader:
            idx_np = idx.numpy().astype(np.int64, copy=False)
            x = x_cpu.to(device)
            target = log_cpu.to(device)
            c_batch = torch.as_tensor(clusterab[idx_np], dtype=torch.float32, device=device)
            # ==== 唯一被改动的变量：corruption 生成 strong/weak 两个视图 ====
            # 两个诊断量对所有 arm 用同一定义，保证可比：
            #   designed_mask = mask.mean()   设计选中率，应 ~= mask_prob（各 arm 相同 -> 干净归因）
            #   eff_change = (strong!=x)       实际数值改动率（稀疏 scRNA 上 swap 因 0<->0 而偏低）
            if corruptor is None:  # zero：完全等同赢家 model.random_mask（保持随机流一致）
                strong, mask = model.random_mask(x, args.mask_prob)
                weak, _ = model.random_mask(x, max(0.05, args.mask_prob * 0.5))
                eff_change = float((strong != x).float().mean().detach().cpu())
            else:                  # swap_*：同箱供体 swap-noise，mask 语义不变
                strong, mask, eff_change = corruptor.corrupt(idx_np, args.mask_prob, device)
                weak, _, _ = corruptor.corrupt(idx_np, max(0.05, args.mask_prob * 0.5), device)
            designed_mask = float(mask.mean().detach().cpu())
            out = model(strong)
            weak_out = model(weak)
            p_batch = None if p_targets is None else torch.as_tensor(p_targets[idx_np], dtype=torch.float32, device=device)
            loss, parts = criterion(out, weak_out, target, mask, p_batch, c_batch, cluster_scale, gate)
            if not torch.isfinite(loss):
                raise FloatingPointError(f"Non-finite loss at epoch {epoch}: {parts}")
            optimizer.zero_grad(set_to_none=True)
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 5.0)
            optimizer.step()

            sums["loss"] += float(loss.detach().cpu())
            sums["scmae_loss"] += float(parts["scmae_loss"])
            sums["sharp_loss"] += float(parts["sharp_loss"])
            sums["variance_loss"] += float(parts["variance_loss"])
            sums["eff_change"] += float(eff_change)
            sums["designed_mask"] += float(designed_mask)
            batches += 1
        for key in sums:
            history[key].append(sums[key] / max(1, batches))
        history["gate"].append(float(gate))
        history["kl_ref"].append(float(kl_ref))
        eff_change_accum.append(history["eff_change"][-1])
        if epoch == 1 or epoch == args.epochs or epoch % 10 == 0:
            print(
                f"Epoch {epoch:03d}/{args.epochs} loss={history['loss'][-1]:.4f} "
                f"scmae={history['scmae_loss'][-1]:.4f} sharp={history['sharp_loss'][-1]:.4f} "
                f"var={history['variance_loss'][-1]:.4f} eff_change={history['eff_change'][-1]:.3f} gate={gate:.3f}"
            )

    embedding, q_out, labels_out = extract_all(model, full_loader, device)
    np.save(save_dir / "embedding_final.npy", embedding.astype(np.float32))
    np.save(save_dir / "labels.npy", labels_out.astype(np.int64))
    save_json(history, str(save_dir / "training_history.json"))
    std_profile = effective_dimensionality(embedding.std(axis=0))
    aligned = cluster_aligned_eff_dim(embedding, labels_out)
    mean_eff_change = float(np.mean(eff_change_accum)) if eff_change_accum else 0.0

    eval_result = None
    preds = None
    if not args.skip_eval:
        eval_result = family.write_kmeans_known_k_outputs(
            save_dir, dataset_name, "conditional-shuffle ablation", args.seed,
            embedding, labels_out, n_clusters,
            {"corruption": args.corruption, "donor_pool_mean": float(pool_stats["pool_mean"]),
             "n_bins": int(pool_stats["n_bins"])},
        )
        preds = eval_result["preds"]["kmeans_known_k"]
        save_json(eval_result["fixed"], str(save_dir / "metrics.json"))

    counts = np.bincount(preds.astype(np.int64), minlength=n_clusters).astype(np.float64) if preds is not None else np.zeros(n_clusters)
    frac = counts / max(1.0, counts.sum())
    summary = {
        "dataset": dataset_name,
        "seed": int(args.seed),
        "n_clusters": int(n_clusters),
        "runtime_seconds": float(time.time() - start),
        "corruption": args.corruption,
        "mask_prob": float(args.mask_prob),
        "mean_effective_change_rate": mean_eff_change,
        "mean_designed_mask_rate": float(np.mean(history["designed_mask"])) if history["designed_mask"] else 0.0,
        "n_nuisance_bins": int(args.n_nuisance_bins),
        "n_joint_bins": int(args.n_joint_bins),
        "donor_pool": pool_stats,
        "nuisance_meta": nuisance["meta"],
        "variance_weight": float(args.variance_weight),
        "force_gate": float(args.force_gate),
        "fixed_metrics": eval_result["fixed"] if eval_result is not None else {},
        "std_profile": std_profile,
        "cluster_aligned": aligned,
        "final_scmae_loss": float(history["scmae_loss"][-1]) if history["scmae_loss"] else 0.0,
        "cluster_mass_min": float(frac.min()) if frac.size else 0.0,
        "cluster_mass_max": float(frac.max()) if frac.size else 0.0,
        "final_gate": float(gate),
    }
    save_json(summary, str(save_dir / "summary.json"))
    ari = summary["fixed_metrics"].get("kmeans_known_k", {}).get("ari")
    nmi = summary["fixed_metrics"].get("kmeans_known_k", {}).get("nmi")
    print(
        f"[RESULT] {dataset_name} corruption={args.corruption} ARI={ari} NMI={nmi} "
        f"eff_change={mean_eff_change:.3f} donor_pool_mean={pool_stats['pool_mean']:.1f} "
        f"aligned_eff_dim={aligned['cluster_aligned_eff_dim']:.1f}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

