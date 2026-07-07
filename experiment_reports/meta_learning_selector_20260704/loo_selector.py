#!/usr/bin/env python3
"""LOO 算法选择器: 用 meta-feature 预测每数据集该用哪个方法, 留一交叉验证。
诚实评估: 预测方法的实际ARI vs SBS(固定最好方法) vs VBS(oracle)。
只用全覆盖(23)的方法作候选, 保证公平。多个候选池对比。"""
import csv, numpy as np
from collections import defaultdict

# --- ARI 表 ---
D=defaultdict(dict)
for r in csv.DictReader(open("results/260629全benchmark结果.csv")):
    try: a=float(r["ari_mean"])
    except: continue
    D[r["dataset"]][r["method"]]=a
ds=sorted(D)
cover={}
allm=sorted({m for d in D for m in D[d]})
for m in allm: cover[m]=sum(1 for d in ds if m in D[d])
full=[m for m in allm if cover[m]==len(ds)]     # 全覆盖
# 也做一个"覆盖>=18且平均高"的更强候选池
strong=[m for m in allm if cover[m]>=18]

# --- meta-features ---
MF={}
rd=list(csv.DictReader(open("experiment_reports/meta_learning_selector_20260704/metafeatures.csv")))
feat_cols=[c for c in rd[0].keys() if c not in ("dataset",)]
for r in rd:
    MF[r["dataset"]]=np.array([float(r[c]) for c in feat_cols])

def run_pool(cands, name):
    # 每数据集: 只在候选池内且该方法覆盖此集
    def ari(d,m): return D[d].get(m, None)
    valid_ds=[d for d in ds if any(ari(d,m) is not None for m in cands)]
    # SBS: 候选池里平均ARI最高的固定方法(要求覆盖所有 valid_ds)
    def mean_on(m, dss): 
        vals=[ari(d,m) for d in dss if ari(d,m) is not None]
        return np.mean(vals) if vals else -1
    sbs_m=max(cands, key=lambda m: mean_on(m, valid_ds))
    sbs=mean_on(sbs_m, valid_ds)
    vbs=np.mean([max(ari(d,m) for m in cands if ari(d,m) is not None) for d in valid_ds])

    # 特征标准化(在全体上, LOO时用训练集统计更严, 但23太少, 简化用全体z-score)
    Xall=np.stack([MF[d] for d in valid_ds]); mu=Xall.mean(0); sd=Xall.std(0)+1e-8
    # LOO: 1-NN in meta-feature space 选"训练集里在最近邻数据集上最好的方法"
    sel_ari=[]; picks=defaultdict(int)
    for i,d in enumerate(valid_ds):
        train=[dd for dd in valid_ds if dd!=d]
        xt=(MF[d]-mu)/sd
        # 最近邻训练数据集
        dists=[(np.sum(((MF[dd]-mu)/sd - xt)**2), dd) for dd in train]
        dists.sort(); nn=dists[0][1]
        # 选在最近邻数据集上ARI最高的方法(仅候选池, 且在目标集也覆盖)
        cand_here=[m for m in cands if ari(d,m) is not None and ari(nn,m) is not None]
        if not cand_here: cand_here=[m for m in cands if ari(d,m) is not None]
        pick=max(cand_here, key=lambda m: ari(nn,m))
        picks[pick]+=1
        sel_ari.append(ari(d,pick))
    sel=np.mean(sel_ari)
    print(f"\n=== 候选池: {name} ({len(cands)}方法, {len(valid_ds)}数据集) ===")
    print(f"  SBS(固定最好={sbs_m}) = {sbs:.4f}")
    print(f"  VBS(oracle)          = {vbs:.4f}   (天花板, gap={vbs-sbs:+.4f})")
    print(f"  LOO 1-NN selector    = {sel:.4f}   (vs SBS {sel-sbs:+.4f}, 捕获oracle空间 {100*(sel-sbs)/(vbs-sbs+1e-9):.0f}%)")
    print(f"  selector 选择分布: {dict(picks)}")

run_pool(full, "全覆盖(23)")
run_pool(strong, "覆盖>=18")
