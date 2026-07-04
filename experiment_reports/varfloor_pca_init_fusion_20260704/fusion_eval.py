#!/usr/bin/env python3
"""Step A: 事后融合评测。对每个 embedding-init 的 run(存了 VF latent + PCA emb),
比较 KMeans-known-K 在 3 种输入上的表现:
  VF   = VarFloor 128d latent
  PCA  = PCA 128d
  FUSE = z-score 后 [PCA ⊕ VF] 拼接(各标准化, 避免尺度失衡)
分层指标: ARI/NMI/macro-F1/稀有簇recall。判断融合是"取最大"还是"插值/稀释"。
只读 center_init=embedding 的 run(原版 VF), 跨种子聚合。"""
from __future__ import annotations
import glob, json
from pathlib import Path
import numpy as np
from scipy.optimize import linear_sum_assignment
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score, f1_score, confusion_matrix


def hmap(labels, pred):
    labs = np.unique(labels); prs = np.unique(pred)
    cm = np.zeros((len(labs), len(prs)), dtype=np.int64)
    li = {l: i for i, l in enumerate(labs)}; pi = {p: j for j, p in enumerate(prs)}
    for l, p in zip(labels, pred): cm[li[l], pi[p]] += 1
    r, c = linear_sum_assignment(-cm)
    mp = {prs[cj]: labs[ri] for ri, cj in zip(r, c)}
    for p in prs:
        if p not in mp: mp[p] = np.bincount(labels[pred == p]).argmax()
    return np.array([mp[p] for p in pred])


def metrics(labels, pred, rare):
    mapped = hmap(labels, pred)
    f1 = f1_score(labels, mapped, average=None, labels=np.unique(labels), zero_division=0)
    rr = [float(np.mean(mapped[labels == c] == c)) for c in rare if (labels == c).any()]
    return {"ari": float(adjusted_rand_score(labels, pred)),
            "nmi": float(normalized_mutual_info_score(labels, pred)),
            "macro_f1": float(np.mean(f1)),
            "rare_recall": float(np.mean(rr)) if rr else float("nan")}


def zscore(x):
    return (x - x.mean(0, keepdims=True)) / (x.std(0, keepdims=True) + 1e-8)


def km(emb, k, seed):
    return KMeans(n_clusters=k, n_init=20, random_state=seed).fit_predict(emb)


def main():
    base = Path(__file__).resolve().parent
    out = base / "fusion_summary.json"
    from collections import defaultdict
    agg = defaultdict(lambda: defaultdict(list))  # dataset -> variant -> [metric dicts]
    n_rare = 5
    # 只取 center_init=embedding 的 run(原版 VF)
    for run in sorted((base / "runs").glob("*__embedding__seed*")):
        vf_p = run / "embedding_final.npy"; pca_p = run / "pca_embedding.npy"; lab_p = run / "labels.npy"
        if not (vf_p.exists() and pca_p.exists() and lab_p.exists()):
            continue
        ds = run.name.split("__")[0]; seed = int(run.name.split("seed")[-1])
        vf = np.nan_to_num(np.load(vf_p).astype(np.float32))
        pca = np.nan_to_num(np.load(pca_p).astype(np.float32))
        labels = np.load(lab_p).astype(np.int64)
        u = np.unique(labels); remap = {x: i for i, x in enumerate(u)}
        labels = np.array([remap[x] for x in labels]); k = len(u)
        rare = np.argsort(np.bincount(labels))[:n_rare]
        fuse = np.concatenate([zscore(pca), zscore(vf)], axis=1)
        for name, emb in [("VF", vf), ("PCA", pca), ("FUSE", fuse)]:
            agg[ds][name].append(metrics(labels, km(emb, k, seed), rare))
        print(f"[{ds} seed{seed}] done")

    # 聚合 mean
    summary = {}
    for ds, variants in agg.items():
        summary[ds] = {}
        for name, rows in variants.items():
            summary[ds][name] = {m: float(np.mean([r[m] for r in rows])) for m in ["ari", "nmi", "macro_f1", "rare_recall"]}
    json.dump(summary, open(out, "w"), indent=2)

    # 打印判定表
    print("\n=== Step A 融合结果 (ARI / macro_f1 / rare_recall, 跨种子均值) ===")
    print(f"{'dataset':20s} {'VF_ari':>7s} {'PCA_ari':>7s} {'FUSE_ari':>8s} {'verdict':>28s}")
    for ds in sorted(summary):
        s = summary[ds]; vf = s["VF"]["ari"]; pca = s["PCA"]["ari"]; fu = s["FUSE"]["ari"]
        mx = max(vf, pca); mn = min(vf, pca)
        if fu >= mx - 0.005: v = "取最大/更好(理想)"
        elif fu <= mn + 0.005: v = "退化到较弱(坏)"
        else: v = "插值(中庸)"
        print(f"{ds:20s} {vf:>7.3f} {pca:>7.3f} {fu:>8.3f} {v:>28s}")
    print(f"\nWrote {out}")


if __name__ == "__main__":
    raise SystemExit(main())
