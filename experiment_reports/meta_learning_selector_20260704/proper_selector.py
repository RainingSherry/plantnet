#!/usr/bin/env python3
"""正式 selector 套件, 严格矩形 Common-18 强池, LODO, 修掉泄露。
分 known-K / deployment(去掉真K特征) 两版。多个 selector + baseline + paired bootstrap CI。"""
import csv, numpy as np
from collections import defaultdict
from sklearn.linear_model import Ridge
from sklearn.tree import DecisionTreeClassifier

# --- ARI 矩阵 ---
D=defaultdict(dict)
for r in csv.DictReader(open("results/260629全benchmark结果.csv")):
    try: a=float(r["ari_mean"])
    except: continue
    D[r["dataset"]][r["method"]]=a
ds_all=sorted(D); methods=sorted({m for d in D for m in D[d]})
covers={m:set(d for d in ds_all if m in D[d]) for m in methods}
strong=[m for m in methods if len(covers[m])>=18]
uni=set(ds_all)
for m in strong: uni&=covers[m]
uni=sorted(uni)                       # 严格矩形18集
M=strong
A=np.array([[D[d][m] for m in M] for d in uni])   # [18, 10]
sbs_j=int(np.argmax(A.mean(0))); SBS=A.mean(0)[sbs_j]
VBS=A.max(1).mean()
print(f"矩形: {len(M)}方法 × {len(uni)}数据集 | SBS={M[sbs_j]}={SBS:.4f} VBS={VBS:.4f} gap={VBS-SBS:+.4f}")

# --- meta-features ---
rd={r["dataset"]:r for r in csv.DictReader(open("experiment_reports/meta_learning_selector_20260704/metafeatures.csv"))}
allcols=[c for c in next(iter(rd.values())).keys() if c!="dataset"]
knownK_cols=allcols                                    # 含 k / cells_per_cluster
deploy_cols=[c for c in allcols if c not in ("k","cells_per_cluster")]  # 去掉真K
def feat(cols): return np.array([[float(rd[d][c]) for c in cols] for d in uni])

def evaluate(selfn, name):
    sel=np.zeros(len(uni))
    for i in range(len(uni)):
        tr=[j for j in range(len(uni)) if j!=i]
        pick=selfn(i, tr)
        sel[i]=A[i, pick]
    gap=VBS-SBS
    capt=(sel.mean()-SBS)/(gap+1e-9)
    # paired bootstrap CI of (sel - always-SBS) per-dataset
    diff=sel - A[:,sbs_j]
    rng=np.random.default_rng(0); bs=[]
    for _ in range(2000):
        idx=rng.integers(0,len(uni),len(uni)); bs.append(diff[idx].mean())
    lo,hi=np.percentile(bs,[2.5,97.5])
    reg=VBS - sel.mean(); wc=float((A.max(1)-sel).max())
    win=float((sel>A[:,sbs_j]+1e-9).mean())
    print(f"  {name:26s} util={sel.mean():.4f} regret={reg:.4f} wc_regret={wc:.4f} "
          f"win={win:.2f} 捕获gap={capt*100:5.0f}% CI[{lo:+.4f},{hi:+.4f}]{' 不跨0' if lo>0 else ''}")

# ---- baselines ----
def always_sbs(i,tr): return sbs_j
def random_sel(i,tr): return np.random.default_rng(i).integers(0,len(M))
def oracle(i,tr): return int(np.argmax(A[i]))
print("\n[baselines]")
evaluate(always_sbs,"always-SBS")
evaluate(random_sel,"random")
evaluate(oracle,"ORACLE(VBS上限)")

for cols,tag in [(deploy_cols,"deployment(去真K)"),(knownK_cols,"known-K")]:
    X=feat(cols)
    print(f"\n[meta-learners: {tag}, {len(cols)}特征]")
    def z(i,tr):
        mu=X[tr].mean(0); sd=X[tr].std(0)+1e-8; return (X-mu)/sd
    # 1-NN LODO
    def knn(i,tr):
        Z=z(i,tr); dists=[(np.sum((Z[j]-Z[i])**2),j) for j in tr]; nn=min(dists)[1]
        return int(np.argmax(A[nn]))
    # ridge performance regression: 每方法 ARI~feat, 预测test各方法ARI取argmax
    def ridge(i,tr):
        Z=z(i,tr); best=None; bestpred=-9
        for mj in range(len(M)):
            r=Ridge(alpha=1.0).fit(Z[tr], A[tr,mj]); p=r.predict(Z[i:i+1])[0]
            if p>bestpred: bestpred=p; best=mj
        return best
    # decision tree classifier on best-method label
    def dtree(i,tr):
        Z=z(i,tr); y=np.array([int(np.argmax(A[j])) for j in tr])
        if len(set(y))<2: return sbs_j
        t=DecisionTreeClassifier(max_depth=3,random_state=0).fit(Z[tr],y)
        return int(t.predict(Z[i:i+1])[0])
    evaluate(knn,"1-NN")
    evaluate(ridge,"ridge-perf-regression")
    evaluate(dtree,"decision-tree(d3)")
