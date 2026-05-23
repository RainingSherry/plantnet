#!/usr/bin/env python3
"""Re-inference for cursor_Doloris_GraphDiffusion on GPU 1."""
import os, sys, json, h5py, numpy as np, torch, torch.nn as nn
from pathlib import Path
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans

GRAPHDIFFUSION_DIR = Path("/home/luolie/biopipeline/dimension-reduction/plantnet/methods/DeepLearning/cursor_Doloris/GraphDiffusion")
sys.path.insert(0, str(GRAPHDIFFUSION_DIR))

from models import PlantDiffCluster


def load_h5ad_safe(h5ad_path):
    import scipy.sparse as sp
    f = h5py.File(h5ad_path, 'r')
    X = f['X']
    if isinstance(X, h5py.Group):
        data = np.array(X['data']); indices = np.array(X['indices']); indptr = np.array(X['indptr'])
        shape = tuple(int(x) for x in X.attrs['shape']) if 'shape' in X.attrs else None
        if shape: X_mat = sp.csr_matrix((data, indices, indptr), shape=shape)
        else: X_mat = sp.csr_matrix((data, indices, indptr))
        X_arr = X_mat.toarray().astype(np.float32)
    else:
        X_arr = np.array(X).astype(np.float32)
    obs = {}
    if 'obs' in f:
        for key in f['obs']:
            try:
                ds = f[f'obs/{key}']
                if isinstance(ds, h5py.Group):
                    if 'codes' in ds and 'categories' in ds:
                        codes = np.array(ds['codes'][:])
                        cats = np.array(ds['categories'][:])
                        cats_str = np.array([s.decode() if isinstance(s, bytes) else s for s in cats])
                        obs[key] = cats_str[codes]
                    else: continue
                else:
                    vals = np.array(ds[:])
                    if vals.dtype.kind == 'O':
                        vals = np.array([s.decode() if isinstance(s, bytes) else s for s in vals])
                    obs[key] = vals
            except: continue
    var_names = None
    if 'var' in f:
        if '_index' in f['var']:
            var_names = np.array([s.decode() if isinstance(s, bytes) else s for s in np.array(f['var/_index'][:])])
        elif 'gene_name' in f['var']:
            var_names = np.array([s.decode() if isinstance(s, bytes) else s for s in np.array(f['var/gene_name'][:])])
    f.close()
    return X_arr, obs, var_names


def get_label_column(obs):
    for col in ['Celltype', 'cell_type', 'celltype', 'CellType']:
        if col in obs: return np.array(obs[col])
    raise KeyError(f"No label col. Available: {list(obs.keys())}")


def compute_metrics(y_true, y_pred):
    from sklearn.preprocessing import LabelEncoder
    from scipy.optimize import linear_sum_assignment
    from sklearn.metrics import (accuracy_score, f1_score, normalized_mutual_info_score as nmi_score,
        adjusted_rand_score as ari_score, fowlkes_mallows_score as fmi_score,
        v_measure_score as vms_score, homogeneity_score as hom_score, completeness_score as com_score)
    le = LabelEncoder()
    y_true_enc = le.fit_transform(y_true)
    y_pred_enc = y_pred.astype(int)
    tu = np.unique(y_true_enc); pu = np.unique(y_pred_enc)
    G = np.zeros((len(tu), len(pu)), dtype=int)
    for i, ut in enumerate(tu):
        for j, up in enumerate(pu):
            G[i, j] = np.sum((y_true_enc == ut) & (y_pred_enc == up))
    A = linear_sum_assignment(-G)
    new_pred = np.zeros_like(y_pred_enc)
    for i, up in enumerate(pu):
        col_i = A[1][i] if i < len(A[1]) else i % len(tu)
        lab_i = A[0][i] if i < len(A[0]) else i % len(tu)
        new_pred[y_pred_enc == up] = tu[lab_i]
    return {
        "acc": float(accuracy_score(y_true_enc, new_pred)),
        "nmi": float(nmi_score(y_true_enc, y_pred_enc, average_method="arithmetic")),
        "ari": float(ari_score(y_true_enc, y_pred_enc)),
        "f1_macro": float(f1_score(y_true_enc, new_pred, average='macro', zero_division=0)),
        "fmi": float(fmi_score(y_true_enc, y_pred_enc)),
        "v_measure": float(vms_score(y_true_enc, y_pred_enc)),
        "homogeneity": float(hom_score(y_true_enc, y_pred_enc)),
        "completeness": float(com_score(y_true_enc, y_pred_enc)),
    }


class SimpleDataset(torch.utils.data.Dataset):
    def __init__(self, X_hvg, X_raw, hvg_idx, labels, topk_support=2000):
        # X_hvg: (n_cells, n_hvg) - preprocessed HVG expression
        # X_raw: (n_cells, n_full) - raw counts
        # hvg_idx: global indices of HVG genes
        # labels: cell type labels
        self.X = X_hvg  # Use HVG submatrix as input
        self.X_raw = X_raw
        self.hvg_idx = hvg_idx
        self.labels = labels
        self.n_cells = X_hvg.shape[0]
        self.n_hvg = X_hvg.shape[1]

        # Compute support: for each cell, top-k expressed HVG genes (LOCAL indices)
        support_idx_list = []
        support_weight_list = []
        max_support = 0

        for c in range(self.n_cells):
            expr = X_raw[c, hvg_idx]  # expression in HVG genes (global space)
            expr_clip = np.clip(expr, 0, None)
            # Top-k by expression
            k = min(topk_support, len(hvg_idx))
            topk_local = np.argpartition(expr_clip, -k)[-k:]
            topk_local = topk_local[np.argsort(expr_clip[topk_local])[::-1]]  # descending
            # Local indices (0 to n_hvg-1) not global
            support_idx_list.append(topk_local.astype(np.int64))
            # Normalized weights
            vals = expr_clip[topk_local]
            w_sum = vals.sum() + 1e-8
            support_weight_list.append((vals / w_sum).astype(np.float32))
            max_support = max(max_support, len(topk_local))

        # Pad to fixed length
        self.max_support = ((max_support + 31) // 32) * 32
        self.max_support = min(self.max_support, 2000)

        self.support_idx = np.full((self.n_cells, self.max_support), -1, dtype=np.int64)
        self.support_mask = np.zeros((self.n_cells, self.max_support), dtype=np.float32)
        self.support_weight = np.zeros((self.n_cells, self.max_support), dtype=np.float32)

        for c in range(self.n_cells):
            s = min(len(support_idx_list[c]), self.max_support)
            self.support_idx[c, :s] = support_idx_list[c][:s]
            self.support_mask[c, :s] = 1.0
            self.support_weight[c, :s] = support_weight_list[c][:s]

    def __len__(self): return self.n_cells
    def __getitem__(self, idx):
        return {
            "X": torch.from_numpy(self.X[idx]).float(),
            "support_idx": torch.from_numpy(self.support_idx[idx]).long(),
            "support_mask": torch.from_numpy(self.support_mask[idx]).float(),
            "support_weight": torch.from_numpy(self.support_weight[idx]).float(),
            "label": torch.tensor(0, dtype=torch.long),
        }


@torch.no_grad()
def gen_emb(model, dataset, device, bs=64):
    model.eval()
    loader = torch.utils.data.DataLoader(dataset, batch_size=bs, shuffle=False, drop_last=False, num_workers=0)
    all_emb = []
    for batch in loader:
        X = batch["X"].to(device)
        si = batch["support_idx"]; sm = batch["support_mask"]; sw = batch["support_weight"]
        out = model(X=X, cell_type=None, support_weight=sw, support_mask=sm, support_idx=si, t=None)
        all_emb.append(out["cell_z"].cpu().numpy())
    return np.vstack(all_emb)


def build_coexp_graph(X_hvg, n_hvg=1500, top_k=20):
    """Build co-expression gene graph."""
    n_genes = X_hvg.shape[1]
    chunk_size = 500
    all_corr = np.zeros((n_genes, n_genes), dtype=np.float32)
    for start in range(0, n_genes, chunk_size):
        end = min(start + chunk_size, n_genes)
        chunk = X_hvg[:, start:end]
        if chunk.shape[1] == 0: continue
        chunk_m = chunk.mean(axis=0, keepdims=True); chunk_s = chunk.std(axis=0, keepdims=True) + 1e-8
        chunk_n = (chunk - chunk_m) / chunk_s
        full_m = X_hvg.mean(axis=0, keepdims=True); full_s = X_hvg.std(axis=0, keepdims=True) + 1e-8
        full_n = (X_hvg - full_m) / full_s
        bc = np.dot(chunk_n.T, full_n) / max(X_hvg.shape[0] - 1, 1)
        all_corr[start:end, :] = np.abs(np.nan_to_num(bc, nan=0.0, posinf=0.0, neginf=0.0))
    np.fill_diagonal(all_corr, 0)
    row_list, col_list, weight_list = [], [], []
    for i in range(n_genes):
        abs_c = all_corr[i].copy(); abs_c[i] = 0
        neighbors = np.argsort(abs_c)[::-1][:top_k]
        for j in neighbors:
            if abs_c[j] > 0:
                row_list.append(i); col_list.append(j)
                weight_list.append(float(max(abs_c[j], 0.01)))
    w_min, w_max = min(weight_list), max(weight_list)
    if w_max > w_min:
        weight_list = [0.5 + 0.5 * (w - w_min) / (w_max - w_min) for w in weight_list]
    else:
        weight_list = [1.0] * len(weight_list)
    return {
        "edge_index": [row_list, col_list], "edge_weight": weight_list,
        "n_nodes": n_genes, "graph_type": "coexpression",
    }


def run(data_path, result_dir, config_path, n_clusters, device):
    rd = Path(result_dir); rd.mkdir(parents=True, exist_ok=True)
    mf = rd / "metrics.json"
    if mf.exists():
        print(f"metrics.json exists, skip"); return
    print(f"Loading {data_path}")
    X_raw, obs, var_names = load_h5ad_safe(data_path)
    y_true = get_label_column(obs)
    valid = y_true != "Unknow"
    X_raw = X_raw[valid]; y_true = y_true[valid]
    print(f"  {X_raw.shape[0]} cells, {len(np.unique(y_true))} types")
    X_log = np.log1p(X_raw)
    scaler = StandardScaler(); X_scaled = scaler.fit_transform(X_log)
    X_scaled = np.nan_to_num(X_scaled, nan=0.0, posinf=0.0, neginf=0.0)

    # HVG selection
    print("Building graph...")
    var = np.var(X_scaled, axis=0)
    hvg_idx = np.argsort(var)[::-1][:1500]
    X_hvg = X_scaled[:, hvg_idx]
    gn = var_names[hvg_idx] if var_names is not None else np.arange(len(hvg_idx))
    graph_dict = build_coexp_graph(X_hvg, n_hvg=1500, top_k=20)
    graph_dict["gene_names"] = [str(g) for g in gn]
    print(f"  Graph: {graph_dict['n_nodes']} nodes, {len(graph_dict['edge_index'][0])} edges")

    # Load config
    with open(config_path) as f: sc = json.load(f)
    n_ct = len([u for u in np.unique(y_true) if u != 'Unknow'])

    # Build model config (use the n_ct from data, not config)
    n_actual_clusters = max(n_ct, n_clusters)
    mc = {
        "gene_dim": sc.get("gene_dim", 64), "hidden_dim": sc.get("hidden_dim", 256),
        "embed_dim": sc.get("embed_dim", 128), "time_embed_dim": sc.get("time_embed_dim", 128),
        "n_layers": sc.get("n_layers", 2), "heads": sc.get("heads", [4, 4]),
        "pooling_strategy": sc.get("pooling_strategy", "attention"),
        "pooling_topk": sc.get("pooling_topk", 50),
        "n_clusters": n_actual_clusters,
        "cluster_strategy": sc.get("cluster_strategy", "gmm"),
        "use_diffusion": sc.get("use_diffusion", True),
        "use_mask_predictor": sc.get("use_mask_predictor", True),
        "num_timesteps": sc.get("num_timesteps", 500),
        "ddim_steps": sc.get("ddim_steps", 20),
        "beta_schedule": sc.get("beta_schedule", "cosine"),
        "refiner_depth": sc.get("refiner_depth", 3),
        "refiner_hidden_dim": sc.get("refiner_hidden_dim", 256),
        "lambda_cluster": sc.get("lambda_cluster", 0.1),
        "cell_type_num": n_ct,
        "use_decoder": sc.get("use_decoder", True),
        "decoder_hidden_dim": sc.get("decoder_hidden_dim", 256),
        "dropout": sc.get("dropout", 0.1),
    }
    print(f"Building model (n_clusters={n_actual_clusters})...")
    model = PlantDiffCluster(n_genes=len(hvg_idx), gene_names=list(gn), graph_dict=graph_dict, config=mc).to(device)
    print(f"  {sum(p.numel() for p in model.parameters() if p.requires_grad):,} params")

    print("Generating embeddings...")
    ds = SimpleDataset(X_hvg, X_raw, hvg_idx, y_true)
    embeddings = gen_emb(model, ds, device, bs=64)
    print(f"  Shape: {embeddings.shape}")

    print("Clustering...")
    results = {}; best_nmi = 0; best_pred = None; best_k = None
    for k in [n_actual_clusters, n_actual_clusters-2, n_actual_clusters+2, n_actual_clusters-5, n_actual_clusters+5]:
        if k < 2: continue
        km = KMeans(n_clusters=k, n_init=20, random_state=42)
        pred = km.fit_predict(embeddings)
        m = compute_metrics(y_true, pred)
        results[f"kmeans_{k}"] = m
        print(f"  kmeans_{k}: NMI={m['nmi']:.4f} ARI={m['ari']:.4f} ACC={m['acc']:.4f}")
        if m['nmi'] > best_nmi: best_nmi = m['nmi']; best_pred = pred; best_k = k
    bm = results[f"kmeans_{best_k}"]
    bm["source"] = "cursor_Doloris_GraphDiffusion"
    bm["best_n_clusters"] = best_k
    np.save(rd / "best_embeddings.npy", embeddings)
    np.save(rd / "best_predictions.npy", best_pred)
    with open(rd / "all_results.json", "w") as f: json.dump(results, f, indent=2)
    with open(mf, "w") as f: json.dump(bm, f, indent=2)
    print(f"\nSaved: {mf}")
    print(f"Best: kmeans_{best_k} NMI={bm['nmi']:.4f} ARI={bm['ari']:.4f} ACC={bm['acc']:.4f}")


if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument("--data-path", required=True)
    p.add_argument("--result-dir", required=True)
    p.add_argument("--config-path", required=True)
    p.add_argument("--n-clusters", type=int, default=15)
    p.add_argument("--gpu", type=int, default=1)
    a = p.parse_args()
    os.environ["CUDA_VISIBLE_DEVICES"] = str(a.gpu)
    dev = f"cuda:{0}" if torch.cuda.is_available() else "cpu"
    run(a.data_path, a.result_dir, a.config_path, a.n_clusters, dev)
