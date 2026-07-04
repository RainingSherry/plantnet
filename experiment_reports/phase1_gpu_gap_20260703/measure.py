#!/usr/bin/env python3
"""Phase 1 —— CPU↔GPU 缺口实测。

对一个代表性深度单细胞方法(默认 AdaptiveSwitch=scMAE+DEC+std-floor 旗舰)的完整流水线
分阶段计时，量化「哪些阶段是 CPU-bound、GPU 化能省多少、随细胞数如何变化」。

分阶段(现有 pipeline 的真实阶段)：
  load        : 读 h5ad + 归一化 + HVG (+scale)   —— CPU (scanpy/numpy)
  neighbor    : PCA(50) + KNN 图                    —— CPU (sklearn)   [经典步骤,cuML可加速]
  kmeans_init : DEC 中心初始化 KMeans(n_init=20)    —— CPU (sklearn)   [经典步骤,cuML可加速]
  train       : scMAE+DEC 训练循环                  —— GPU (torch)     [已GPU]
  eval_kmeans : 最终 KMeans(known-k, n_init=20)     —— CPU (sklearn)   [经典步骤,cuML可加速]
输出每阶段 wall-time + GPU 利用率(后台采样) + 峰值显存 + 峰值 CPU 内存。
device=cpu 时训练也在 CPU，用来对比"若无 GPU"的训练代价。

用法示例:
  python measure.py --data_path <h5ad> --dataset_name X --n_clusters K --device cuda --gpu 1
  python measure.py ... --device cpu --cpu_threads 16   (CPU baseline, 绑定线程数避免超订)
"""
from __future__ import annotations
import argparse, json, os, sys, time, threading
from pathlib import Path


def _pin_threads(n: int):
    for v in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
        os.environ[v] = str(n)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data_path", required=True)
    ap.add_argument("--save_dir", required=True)
    ap.add_argument("--dataset_name", default=None)
    ap.add_argument("--n_clusters", type=int, required=True)
    ap.add_argument("--label_key", default="auto")
    ap.add_argument("--n_top_genes", type=int, default=1000)
    ap.add_argument("--epochs", type=int, default=30)  # 计时用途,不需80
    ap.add_argument("--warmup_epochs", type=int, default=10)
    ap.add_argument("--batch_size", type=int, default=256)
    ap.add_argument("--device", choices=["cuda", "cpu"], default="cuda")
    ap.add_argument("--gpu", type=int, default=1)
    ap.add_argument("--cpu_threads", type=int, default=16)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    _pin_threads(args.cpu_threads)
    os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")
    if args.device == "cuda":
        os.environ["CUDA_VISIBLE_DEVICES"] = str(args.gpu)

    import numpy as np
    import torch
    import torch.nn.functional as F
    from sklearn.cluster import KMeans
    from sklearn.decomposition import PCA
    from sklearn.neighbors import NearestNeighbors
    from sklearn.preprocessing import normalize
    from torch.utils.data import DataLoader, Dataset

    ROOT = Path(__file__).resolve().parents[2]
    ADSW = ROOT / "experimental_retired_models" / "Granularity_scMAE_experiments" / "AdaptiveSwitch_scMAE"
    for p in (str(ROOT), str(ADSW)):
        if p not in sys.path:
            sys.path.insert(0, p)
    from loss import AdaptiveSwitchLoss, compute_gate
    from model import AdaptiveSwitchScMAE
    from clusterability import compute_clusterability
    from methods.DeepLearning import scMAE_family as family
    from methods.shared_utils import ensure_dir, save_json

    # null-h5ad shim
    try:
        import h5py
        from anndata._io.specs.registry import _REGISTRY, IOSpec
        for typ in (h5py.Dataset, h5py.Group):
            try:
                _REGISTRY.register_read(typ, IOSpec("null", "0.1.0"))(lambda *a, **k: None)
            except Exception:
                pass
    except Exception:
        pass

    dev = torch.device("cuda:0" if args.device == "cuda" and torch.cuda.is_available() else "cpu")
    save_dir = Path(ensure_dir(args.save_dir))
    ds_name = args.dataset_name or Path(args.data_path).stem
    timings = {}
    T = lambda: time.perf_counter()

    # ---- GPU 利用率后台采样器 ----
    gpu_util, stop = [], threading.Event()
    def sampler():
        import subprocess
        while not stop.is_set():
            try:
                out = subprocess.run(["nvidia-smi", "--query-gpu=utilization.gpu,memory.used",
                                      "--format=csv,noheader,nounits", "-i", str(args.gpu)],
                                     capture_output=True, text=True, timeout=2).stdout.strip()
                u, m = out.split(",")
                gpu_util.append((float(u), float(m)))
            except Exception:
                pass
            stop.wait(0.5)
    if args.device == "cuda":
        threading.Thread(target=sampler, daemon=True).start()

    # ===== STAGE: load =====
    t0 = T()
    tgt = family.load_scmae_dataset(args.data_path, "auto", args.n_top_genes, 10000.0, False, args.label_key, args.seed)
    enc_b = family.load_scmae_dataset(args.data_path, "auto", args.n_top_genes, 10000.0, True, args.label_key, args.seed)
    enc = np.asarray(enc_b.data, dtype=np.float32)
    log_expr = np.asarray(tgt.data, dtype=np.float32)
    labels = np.asarray(tgt.labels, dtype=np.int64)
    n_clusters = int(args.n_clusters if args.n_clusters > 0 else len(np.unique(labels)))
    timings["load"] = T() - t0
    n_cells, n_genes = enc.shape

    # ===== STAGE: neighbor (PCA + KNN) —— 经典步骤 =====
    t0 = T()
    dim = max(2, min(50, n_cells - 1, n_genes - 1))
    emb_pca = PCA(n_components=dim, random_state=args.seed).fit_transform(enc.astype(np.float64))
    emb_pca = normalize(emb_pca, norm="l2", axis=1)
    nn_ = NearestNeighbors(n_neighbors=16, metric="cosine").fit(emb_pca)
    _, nb_idx = nn_.kneighbors(emb_pca, return_distance=True)
    nb_idx = nb_idx[:, 1:].astype(np.int64)
    timings["neighbor_pca_knn"] = T() - t0

    # ===== build model =====
    class DS(Dataset):
        def __init__(s, e, l, y): s.e=torch.as_tensor(e); s.l=torch.as_tensor(l); s.y=torch.as_tensor(y)
        def __len__(s): return int(s.e.shape[0])
        def __getitem__(s, i): return int(i), s.e[i], s.l[i], s.y[i]
    ds = DS(enc, log_expr, labels)
    gen = torch.Generator().manual_seed(args.seed)
    train_loader = DataLoader(ds, batch_size=args.batch_size, shuffle=True, drop_last=True, generator=gen)
    full_loader = DataLoader(ds, batch_size=max(args.batch_size*4, 512), shuffle=False)
    model = AdaptiveSwitchScMAE(n_genes, n_clusters, 128, 0.05).to(dev)
    crit = AdaptiveSwitchLoss(0.75, 0.65, 0.35, 0.05, 0.02, 0.10, 0.35, "hinge")
    opt = torch.optim.Adam(model.parameters(), lr=1e-3, weight_decay=1e-5)

    @torch.no_grad()
    def extract():
        model.eval(); E=[]; Q=[]
        for _, x, _, _ in full_loader:
            o = model(x.to(dev)); E.append(o["latent"].cpu().numpy()); Q.append(o["cluster_q"].cpu().numpy())
        return np.concatenate(E).astype(np.float32), np.concatenate(Q).astype(np.float32)

    # ===== STAGE: train (含首个 kmeans_init 单独计时) =====
    p_targets=None; clusterab=np.ones(n_cells,dtype=np.float32); gate=1.0; inited=False
    kmeans_init_time=0.0; train_epoch_times=[]
    for epoch in range(1, args.epochs+1):
        if epoch > args.warmup_epochs and (p_targets is None or (epoch-args.warmup_epochs-1) % 5 == 0):
            emb, q = extract()
            if not inited:
                tk=T()
                km=KMeans(n_clusters=n_clusters, n_init=20, random_state=args.seed).fit(emb)
                kmeans_init_time=T()-tk
                model.initialize_centers(torch.as_tensor(km.cluster_centers_, dtype=torch.float32, device=dev))
                emb, q = extract(); inited=True
            sp=AdaptiveSwitchScMAE.sharpen(torch.as_tensor(q)).numpy().astype(np.float32); p_targets=sp
            clusterab,_=compute_clusterability(emb, q, nb_idx, k=15)
            gate,_=compute_gate(sp, q, 0.15)
        cs = 0.0 if epoch<=args.warmup_epochs else min(1.0,(epoch-args.warmup_epochs)/max(1,args.warmup_epochs))
        model.train(); te=T()
        for idx, x, lg, _ in train_loader:
            idx=idx.numpy(); x=x.to(dev); tgt_b=lg.to(dev)
            cb=torch.as_tensor(clusterab[idx], dtype=torch.float32, device=dev)
            strong, mask = model.random_mask(x, 0.4); weak,_=model.random_mask(x, 0.2)
            out=model(strong); wout=model(weak)
            pb=None if p_targets is None else torch.as_tensor(p_targets[idx], dtype=torch.float32, device=dev)
            loss,_=crit(out, wout, tgt_b, mask, pb, cb, cs, gate)
            opt.zero_grad(set_to_none=True); loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(),5.0); opt.step()
        train_epoch_times.append(T()-te)
    timings["kmeans_init"]=kmeans_init_time
    timings["train_total"]=float(np.sum(train_epoch_times))
    timings["train_per_epoch"]=float(np.mean(train_epoch_times))

    # ===== STAGE: eval_kmeans —— 经典步骤 =====
    emb,_=extract()
    t0=T()
    KMeans(n_clusters=n_clusters, n_init=20, random_state=args.seed).fit_predict(emb)
    timings["eval_kmeans"]=T()-t0

    stop.set(); time.sleep(0.6)
    peak_mem_mb = float(torch.cuda.max_memory_allocated(0)/1e6) if dev.type=="cuda" else 0.0
    utils=[u for u,_ in gpu_util]; mems=[m for _,m in gpu_util]
    total=sum(v for k,v in timings.items() if k in
              ["load","neighbor_pca_knn","kmeans_init","train_total","eval_kmeans"])
    classic=timings["neighbor_pca_knn"]+timings["kmeans_init"]+timings["eval_kmeans"]+timings["load"]
    summary={
        "dataset":ds_name,"device":args.device,"n_cells":int(n_cells),"n_genes":int(n_genes),
        "n_clusters":n_clusters,"epochs":args.epochs,"cpu_threads":args.cpu_threads,
        "timings_sec":timings,
        "total_measured_sec":total,
        "cpu_classic_sec":classic,
        "cpu_classic_frac":classic/max(total,1e-9),
        "gpu_util_mean":float(np.mean(utils)) if utils else None,
        "gpu_util_p90":float(np.percentile(utils,90)) if utils else None,
        "gpu_mem_peak_mb_smi":float(np.max(mems)) if mems else None,
        "torch_peak_mem_mb":peak_mem_mb,
    }
    save_json(summary, str(save_dir/"gap.json"))
    print(f"[GAP] {ds_name} dev={args.device} n={n_cells} | load={timings['load']:.1f} "
          f"nbr={timings['neighbor_pca_knn']:.1f} kmInit={timings['kmeans_init']:.1f} "
          f"train={timings['train_total']:.1f}({timings['train_per_epoch']:.2f}/ep) "
          f"evalKM={timings['eval_kmeans']:.1f} | classic_frac={summary['cpu_classic_frac']:.0%} "
          f"gpu_util={summary['gpu_util_mean']} peakMemMB={summary['gpu_mem_peak_mb_smi']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
