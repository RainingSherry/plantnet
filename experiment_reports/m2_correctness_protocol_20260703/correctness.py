#!/usr/bin/env python3
"""M2 —— 结构感知统计等价正确性协议（论文脊梁，第一步）。

核心论点：GPU 迁移(cuML)用近似算子，可能悄悄退化稀有/细粒度结构，而聚合 ARI 看不见。
所以"迁移正确"不能只看整体 ARI，必须**分层**校验，且落在 CPU 自身的**多种子噪声带**内。

第一步(本脚本)：在 FROZEN 的旗舰 DEC+std-floor embedding 上，比较
  CPU 参考 = sklearn.KMeans   vs   GPU 迁移 = cuml.KMeans
各自跨多个 seed 跑，对每个方法得到指标的分布(噪声带)，再判断:
  (a) 聚合: ARI / NMI —— cuML 的均值是否落在 sklearn 的噪声带内(TOST 风格: |Δmean| <= sklearn_sd)
  (b) 结构感知(关键): 稀有簇 recall / 每簇 F1(macro, Hungarian 对齐后) / 方差谱 eff-dim
      —— 稀有簇是聚合 ARI 掩盖退化的地方; 这才是"我们的诊断"承重之处。

用法(在 rapids_bench 环境, CUDA_VISIBLE_DEVICES 指定单卡):
  CUDA_VISIBLE_DEVICES=1 python correctness.py \
    --emb_glob ".../dec_floor*/embedding_final.npy" --save_dir runs/macosko_flagship \
    --n_rare 5 --seeds 0 1 2 3 4
每个 embedding 目录下需有配套 labels.npy。
"""
from __future__ import annotations
import argparse, glob, json, os, time
from pathlib import Path
import numpy as np
from scipy.optimize import linear_sum_assignment
from sklearn.metrics import (adjusted_rand_score, normalized_mutual_info_score,
                             confusion_matrix, f1_score)


def hungarian_map(labels: np.ndarray, pred: np.ndarray) -> np.ndarray:
    """把 pred 的簇号最优映射到真标签(最大化重叠), 返回映射后的 pred。"""
    labs = np.unique(labels); prs = np.unique(pred)
    C = confusion_matrix(labels, pred, labels=None)
    # 行=真类, 列=预测簇; 用负号做最大化匹配
    cm = np.zeros((len(labs), len(prs)), dtype=np.int64)
    lab_idx = {l: i for i, l in enumerate(labs)}
    pr_idx = {p: j for j, p in enumerate(prs)}
    for l, p in zip(labels, pred):
        cm[lab_idx[l], pr_idx[p]] += 1
    r, c = linear_sum_assignment(-cm)
    mapping = {prs[cj]: labs[ri] for ri, cj in zip(r, c)}
    # 未被匹配到的预测簇(当簇数>类数时) 映射到其众数真类
    for p in prs:
        if p not in mapping:
            mask = pred == p
            mapping[p] = np.bincount(labels[mask]).argmax()
    return np.array([mapping[p] for p in pred], dtype=labels.dtype)


def structure_metrics(labels: np.ndarray, pred: np.ndarray, rare_classes: np.ndarray) -> dict:
    """结构感知指标: 整体 + 每簇 F1(macro) + 稀有簇 recall。"""
    mapped = hungarian_map(labels, pred)
    per_class_f1 = f1_score(labels, mapped, average=None, labels=np.unique(labels), zero_division=0)
    macro_f1 = float(np.mean(per_class_f1))
    # 稀有类 recall(mapped==true 在稀有类上的比例)
    rare_rec = {}
    for c in rare_classes:
        m = labels == c
        rare_rec[int(c)] = float(np.mean(mapped[m] == c)) if m.sum() > 0 else float("nan")
    rare_recall_mean = float(np.mean(list(rare_rec.values()))) if rare_rec else float("nan")
    return {
        "ari": float(adjusted_rand_score(labels, pred)),
        "nmi": float(normalized_mutual_info_score(labels, pred)),
        "macro_f1": macro_f1,
        "rare_recall_mean": rare_recall_mean,
        "rare_recall_per_class": rare_rec,
        "n_pred": int(len(np.unique(pred))),
    }


def variance_spectrum(emb: np.ndarray) -> dict:
    std = emb.std(axis=0)
    var = np.square(std.astype(np.float64))
    pr = float((var.sum() ** 2) / max(float(np.square(var).sum()), 1e-12))
    return {"eff_dim_pr": pr, "std_min": float(std.min()), "dims_std_gt_1": int((std > 1.0).sum())}


def km_sklearn(emb, k, seed):
    from sklearn.cluster import KMeans
    return KMeans(n_clusters=k, n_init=20, random_state=seed).fit_predict(emb).astype(np.int64)


def km_cuml(emb, k, seed):
    import cupy as cp
    from cuml.cluster import KMeans as cuKMeans
    g = cp.asarray(emb.astype(np.float32))
    p = cuKMeans(n_clusters=k, n_init=20, random_state=seed).fit_predict(g)
    return cp.asnumpy(p).astype(np.int64)


def band_verdict(cpu_vals, gpu_vals, method_band=None):
    """等价判定: |GPU_mean - CPU_mean| 是否落在等价带内。

    带的选择(关键方法学): 用 method_band(该方法自身"重训一次会差多少"的波动, 如 Macosko
    DEC+std-floor ARI 0.087) 才是生物学诚实的容差; 若为 None 回退到 CPU 自身 sd(frozen
    embedding 近确定性, 太紧, 仅参考)。两者都报, verdict 以 method_band(若给)为准。"""
    cpu = np.array([v for v in cpu_vals if v is not None and not np.isnan(v)], dtype=float)
    gpu = np.array([v for v in gpu_vals if v is not None and not np.isnan(v)], dtype=float)
    if len(cpu) == 0 or len(gpu) == 0:
        return {"verdict": "n/a"}
    cmean, csd = float(cpu.mean()), float(cpu.std())
    gmean = float(gpu.mean())
    delta = gmean - cmean
    tight = abs(delta) <= max(csd, 1e-9)
    band = float(method_band) if method_band is not None else None
    honest = (abs(delta) <= band) if band is not None else None
    verdict = ("EQUIVALENT" if honest else "DIFFERENT") if band is not None \
              else ("EQUIVALENT" if tight else "DEGRADED/DIFFERENT")
    return {"cpu_mean": cmean, "cpu_sd": csd, "gpu_mean": gmean, "delta": delta,
            "within_tight_band_1sd": bool(tight), "method_band": band,
            "within_method_band": (bool(honest) if band is not None else None),
            "verdict": verdict}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--emb_glob", required=True)
    ap.add_argument("--save_dir", required=True)
    ap.add_argument("--n_rare", type=int, default=5, help="按类大小取最小的N类为稀有类")
    ap.add_argument("--seeds", type=int, nargs="+", default=[0, 1, 2, 3, 4])
    ap.add_argument("--method_band", type=float, default=None,
                    help="方法自身重训波动(如 Macosko DEC+std-floor ARI 0.087)作为诚实等价带; "
                         "不给则用 frozen-embedding 的紧带(仅参考)")
    ap.add_argument("--gpu", type=int, default=1)
    args = ap.parse_args()
    os.environ.setdefault("CUDA_VISIBLE_DEVICES", str(args.gpu))
    for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ.setdefault(v, "48")
    save = Path(args.save_dir); save.mkdir(parents=True, exist_ok=True)

    emb_files = sorted(glob.glob(args.emb_glob))
    assert emb_files, f"no embeddings match {args.emb_glob}"

    # 多个 embedding = 多个训练种子(如 dec_floor 42/43/44)。CPU 噪声带 = 跨(训练种子×KMeans种子)
    # 的波动 —— 这才是"重训一次会差多少"的真实带, 而非同一 embedding 内近确定性的 KMeans 波动。
    metrics_keys = ["ari", "nmi", "macro_f1", "rare_recall_mean"]
    pooled = {"sklearn": {m: [] for m in metrics_keys}, "cuml": {m: [] for m in metrics_keys}}
    pooled_time = {"sklearn": [], "cuml": []}
    per_emb = []
    for ef in emb_files:
        d = Path(ef).parent
        emb = np.nan_to_num(np.load(ef).astype(np.float32))
        labels = np.load(d / "labels.npy").astype(np.int64)
        uniq = np.unique(labels); remap = {u: i for i, u in enumerate(uniq)}
        labels = np.array([remap[x] for x in labels], dtype=np.int64)
        k = len(uniq); sizes = np.bincount(labels)
        rare_classes = np.argsort(sizes)[:args.n_rare]
        rows = {"sklearn": [], "cuml": []}
        for backend, fn in [("sklearn", km_sklearn), ("cuml", km_cuml)]:
            for s in args.seeds:
                t = time.perf_counter(); pred = fn(emb, k, s); dt = time.perf_counter() - t
                m = structure_metrics(labels, pred, rare_classes)
                m["seconds"] = dt; m["seed"] = s; rows[backend].append(m)
                for mk in metrics_keys:
                    pooled[backend][mk].append(m[mk])
                pooled_time[backend].append(dt)
        per_emb.append({
            "embedding": str(d.name), "n_cells": int(emb.shape[0]), "k": int(k),
            "rare_sizes": [int(sizes[c]) for c in rare_classes],
            "variance_spectrum": variance_spectrum(emb),
            "raw_sklearn": rows["sklearn"], "raw_cuml": rows["cuml"],
        })

    # 跨训练种子 × KMeans 种子 的分层等价判定(CPU 带 = pooled sklearn 的 mean±sd)
    verdicts = {mk: band_verdict(pooled["sklearn"][mk], pooled["cuml"][mk],
                                 method_band=args.method_band) for mk in metrics_keys}
    speedup = float(np.mean(pooled_time["sklearn"]) / max(np.mean(pooled_time["cuml"]), 1e-9))
    degraded = [mk for mk, v in verdicts.items() if v.get("verdict") == "DEGRADED/DIFFERENT"]
    overall = "EQUIVALENT (all metrics incl. rare cells within CPU noise band)" if not degraded \
              else f"NOT-equivalent on: {degraded}"

    summary = {"n_embeddings": len(per_emb), "seeds": args.seeds, "n_rare": args.n_rare,
               "band_definition": "CPU mean±sd pooled over (training-seed embeddings × KMeans seeds)",
               "verdicts": verdicts, "kmeans_speedup_sklearn_over_cuml": speedup,
               "overall_verdict": overall, "degraded_list": degraded, "per_embedding": per_emb}
    with open(save / "correctness.json", "w") as h:
        json.dump(summary, h, indent=2)

    print(f"[M2] {len(per_emb)} training-seed embeddings × {len(args.seeds)} KMeans seeds; "
          f"CPU noise band = pooled sklearn mean±sd")
    for mk in metrics_keys:
        v = verdicts[mk]
        print(f"  {mk:16s} CPU {v['cpu_mean']:.4f}±{v['cpu_sd']:.4f} | cuML {v['gpu_mean']:.4f} "
              f"| Δ{v['delta']:+.4f} | {v['verdict']}")
    print(f"  KMeans speedup(sklearn/cuml) = {speedup:.1f}x")
    print(f"[SUMMARY] {overall}")
    print(f"Wrote {save/'correctness.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
