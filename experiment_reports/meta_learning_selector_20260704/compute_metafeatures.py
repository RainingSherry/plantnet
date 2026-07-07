#!/usr/bin/env python3
"""算23数据集的 label-free meta-feature(不用标签, 除了K=已知簇数是给定的)。
输出 metafeatures.csv 供 LOO 算法选择用。"""
import json, csv, os, sys, warnings
import numpy as np
warnings.filterwarnings("ignore")
sys.path.insert(0, "/home/luolie/biopipeline/dimension-reduction/plantnet")
# null-h5ad shim
try:
    import h5py
    from anndata._io.specs.registry import _REGISTRY, IOSpec
    for t in (h5py.Dataset, h5py.Group):
        try: _REGISTRY.register_read(t, IOSpec("null","0.1.0"))(lambda *a,**k: None)
        except Exception: pass
except Exception: pass
import scanpy as sc, scipy.sparse as sp
from sklearn.decomposition import PCA

M=json.load(open("experiment_reports/formal_varfloor_pca_completion_20260704/task_manifest.json"))
seen={}
for e in M:
    if e["dataset"] not in seen: seen[e["dataset"]]=(e["data_path"], e["n_clusters"])

def entropy(p):
    p=p[p>0]; return float(-(p*np.log(p)).sum())

rows=[]
for d,(path,k) in seen.items():
    a=sc.read_h5ad(path)
    X=a.X
    if sp.issparse(X): X=X.tocsr()
    N,G=X.shape
    # 稀疏度
    nnz=X.nnz if sp.issparse(X) else int((X!=0).sum())
    sparsity=1.0-nnz/(N*G)
    # 每细胞检测基因数 & library
    if sp.issparse(X):
        ndet=np.asarray((X>0).sum(1)).ravel(); lib=np.asarray(X.sum(1)).ravel()
    else:
        ndet=(X>0).sum(1); lib=X.sum(1)
    # HVG 1000 -> PCA 谱形状(内在维度代理), 在 log1p+scale 后
    work=a.copy()
    try:
        if float(np.nanmax(X.max() if not sp.issparse(X) else X.data.max()))>30:
            sc.pp.normalize_total(work, target_sum=1e4); sc.pp.log1p(work)
    except Exception: pass
    try:
        sc.pp.highly_variable_genes(work, n_top_genes=min(1000,G-1), flavor="seurat")
        work=work[:, work.var.highly_variable]
    except Exception: pass
    Xd=work.X.toarray() if sp.issparse(work.X) else np.asarray(work.X)
    Xd=np.nan_to_num(Xd.astype(np.float64))
    Xd=(Xd-Xd.mean(0))/(Xd.std(0)+1e-8)
    npc=min(50, Xd.shape[0]-1, Xd.shape[1]-1)
    pca=PCA(n_components=npc, random_state=0).fit(Xd)
    evr=pca.explained_variance_ratio_
    # 谱形状: участие ratio(有效维), 前10 PC 累计方差, 前1 PC占比
    part_ratio=float((evr.sum()**2)/(np.square(evr).sum()+1e-12))
    top10=float(evr[:10].sum()); pc1=float(evr[0])
    rows.append({
        "dataset":d, "n_cells":N, "n_genes":G, "k":k, "log_n_cells":float(np.log10(N)),
        "cells_per_cluster":N/k, "sparsity":round(sparsity,4),
        "ndet_median":float(np.median(ndet)), "ndet_cv":float(np.std(ndet)/(np.mean(ndet)+1e-8)),
        "lib_cv":float(np.std(lib)/(np.mean(lib)+1e-8)),
        "pca_part_ratio":round(part_ratio,3), "pca_top10_var":round(top10,4),
        "pca_pc1_var":round(pc1,4),
    })
    print(f"{d:28s} N={N} G={G} k={k} sparsity={sparsity:.3f} pca_pr={part_ratio:.1f} pc1={pc1:.3f}", flush=True)

out="experiment_reports/meta_learning_selector_20260704/metafeatures.csv"
with open(out,"w",newline="") as h:
    w=csv.DictWriter(h, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
print("Wrote", out)
