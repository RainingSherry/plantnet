#!/usr/bin/env python3
"""Analyze benchmark results from PlantSPADE_LGCL_protocol - Detailed version."""

import json
import os
from pathlib import Path
from collections import defaultdict

BASE_DIR = Path("/home/luolie/biopipeline/dimension-reduction/plantnet/results/PlantSPADE_LGCL_protocol")

METHODS = [
    "traditional_pca", "phytocluster", "scvi", "scmae",
    "plantspade_lgcl_baseline", "plantspade_lgcl_support_attention"
]

DATASETS = [
    "SRP182008", "SRP235541", "SRP171040", "SRP309176",
    "SRP145013", "CRA002977_1", "SRP224648", "CRA007122"
]

# Expected clusters for each dataset
EXPECTED_CLUSTERS = {
    "SRP182008": 15,
    "SRP235541": 18,
    "SRP171040": 12,
    "SRP309176": 13,
    "SRP145013": 9,
    "CRA002977_1": 7,
    "SRP224648": 4,
    "CRA007122": 7,
}

SEEDS = [1, 2]

def check_log_status(dataset, method, seed):
    """Check log file for success/failure status."""
    log_path = BASE_DIR / dataset / method / f"seed_{seed}" / "logs" / "run_suite.log"
    if not log_path.exists():
        return "NO_LOG", None
    
    try:
        with open(log_path) as f:
            content = f.read()
        
        if "Traceback (most recent call last)" in content:
            lines = content.strip().split('\n')
            error_lines = [l for l in lines if 'Error' in l or 'Exception' in l or 'ValueError' in l or 'Traceback' in l]
            return "FAILED", '\n'.join(error_lines[-5:]) if error_lines else "Unknown error"
        
        if "CUDA out of memory" in content:
            return "OOM", "CUDA Out of Memory"
        
        if "done] ok" in content:
            return "OK", None
        
        # Check for subprocess exit codes
        if "rc=1" in content:
            return "FAILED", "Exit code 1"
        if "rc=-11" in content:
            return "FAILED", "Segmentation fault (rc=-11)"
        
        return "UNKNOWN", "Cannot determine status"
        
    except Exception as e:
        return "ERROR", str(e)

def read_metrics(dataset, method, seed):
    """Read metrics.json if exists."""
    path = BASE_DIR / dataset / method / f"seed_{seed}" / "metrics.json"
    if not path.exists():
        return None
    try:
        with open(path) as f:
            return json.load(f)
    except Exception as e:
        return {"error": str(e)}

def get_primary_metrics(metrics):
    """Extract primary metrics from metrics.json.
    
    Different methods store metrics differently:
    - scvi/scmae: flat structure with acc, nmi, ari, f1_macro
    - plantspade_lgcl_*: nested structure with kmeans_known_k, leiden_fixed, etc.
    """
    if not metrics or "error" in metrics:
        return None
    
    # If flat structure (scvi, scmae style)
    if "nmi" in metrics and "ari" in metrics:
        return {
            "acc": metrics.get("acc"),
            "nmi": metrics.get("nmi"),
            "ari": metrics.get("ari"),
            "f1_macro": metrics.get("f1_macro"),
        }
    
    # If nested structure (plantspade style), prefer kmeans_known_k
    if "kmeans_known_k" in metrics:
        m = metrics["kmeans_known_k"]
        return {
            "acc": m.get("acc"),
            "nmi": m.get("nmi"),
            "ari": m.get("ari"),
            "f1_macro": m.get("f1_macro"),
        }
    
    # Try other keys
    for key in ["leiden_fixed", "louvain_fixed"]:
        if key in metrics:
            m = metrics[key]
            return {
                "acc": m.get("acc"),
                "nmi": m.get("nmi"),
                "ari": m.get("ari"),
                "f1_macro": m.get("f1_macro"),
            }
    
    return None

def main():
    print("=" * 100)
    print("PlantSPADE_LGCL Benchmark Results Summary")
    print("Run command: --datasets all --seeds 1,2 --gpus 1,2,3,4,5,6 --methods all")
    print("=" * 100)
    
    # Build status matrix
    results = defaultdict(lambda: defaultdict(dict))
    
    for dataset in DATASETS:
        for method in METHODS:
            for seed in SEEDS:
                status, error = check_log_status(dataset, method, seed)
                metrics = read_metrics(dataset, method, seed)
                results[dataset][method][seed] = {
                    "status": status,
                    "error": error,
                    "metrics": metrics,
                    "primary_metrics": get_primary_metrics(metrics)
                }
    
    # =========================================================================
    # PART 1: STATUS OVERVIEW
    # =========================================================================
    print("\n" + "=" * 100)
    print("PART 1: STATUS OVERVIEW")
    print("=" * 100)
    
    # Table header
    method_short = {
        "traditional_pca": "pca",
        "phytocluster": "phyt",
        "scvi": "scvi",
        "scmae": "scmae",
        "plantspade_lgcl_baseline": "lgcl_base",
        "plantspade_lgcl_support_attention": "lgcl_attn",
    }
    
    print(f"\n{'Dataset':<15} | {'Species':<20} | K | " + " | ".join([f"{method_short[m]:^10}" for m in METHODS]))
    print("-" * 100)
    
    for dataset in DATASETS:
        species = {
            "SRP182008": "A. thaliana Root",
            "SRP235541": "A. thaliana Root",
            "SRP171040": "A. thaliana Root",
            "SRP309176": "O. sativa Root",
            "SRP145013": "Z. mays Root",
            "CRA002977_1": "A. thaliana Leaf",
            "SRP224648": "Z. mays Leaf",
            "CRA007122": "G. max Nodule",
        }[dataset]
        
        statuses = []
        for method in METHODS:
            s_list = [results[dataset][method][s]["status"][0] for s in SEEDS]
            if s_list == ["O", "O"]:
                statuses.append("\033[92m✓\033[0m")
            elif s_list == ["F", "F"]:
                statuses.append("\033[91m✗\033[0m")
            elif "O" in s_list:
                statuses.append("\033[93m½\033[0m")
            else:
                statuses.append("\033[90m?\033[0m")
        
        print(f"{dataset:<15} | {species:<20} | {EXPECTED_CLUSTERS[dataset]:<2} | " + " | ".join(statuses))
    
    # =========================================================================
    # PART 2: SUCCESSFUL RUNS METRICS
    # =========================================================================
    print("\n\n" + "=" * 100)
    print("PART 2: SUCCESSFUL RUNS - PRIMARY METRICS (K-Means, k=known_k)")
    print("=" * 100)
    print("(Metrics: ARI = Adjusted Rand Index, NMI = Normalized Mutual Information, F1 = F1-macro)")
    
    for dataset in DATASETS:
        has_success = any(
            results[dataset][m][s]["status"] == "OK"
            for m in METHODS for s in SEEDS
        )
        if not has_success:
            continue
            
        print(f"\n### {dataset} (Expected K={EXPECTED_CLUSTERS[dataset]})")
        print("-" * 95)
        print(f"{'Method':<30} | {'Seed':>4} | {'ARI':>8} | {'NMI':>8} | {'F1':>8} | {'ACC':>8}")
        print("-" * 95)
        
        for method in METHODS:
            for seed in SEEDS:
                r = results[dataset][method][seed]
                pm = r["primary_metrics"]
                
                if r["status"] == "OK" and pm:
                    ari = f"{pm.get('ari', 0):.4f}" if pm.get('ari') else "N/A"
                    nmi = f"{pm.get('nmi', 0):.4f}" if pm.get('nmi') else "N/A"
                    f1 = f"{pm.get('f1_macro', 0):.4f}" if pm.get('f1_macro') else "N/A"
                    acc = f"{pm.get('acc', 0):.4f}" if pm.get('acc') else "N/A"
                    print(f"{method:<30} | {seed:>4} | {ari:>8} | {nmi:>8} | {f1:>8} | {acc:>8}")
                else:
                    status_str = r["status"] if r["status"] != "NO_LOG" else "not_run"
                    print(f"{method:<30} | {seed:>4} | \033[90m{status_str}\033[0m")
        print()
    
    # =========================================================================
    # PART 3: MEAN ± STD ACROSS SEEDS
    # =========================================================================
    print("\n" + "=" * 100)
    print("PART 3: AGGREGATED METRICS (Mean ± Std across seeds)")
    print("=" * 100)
    
    for dataset in DATASETS:
        has_success = any(
            results[dataset][m][s]["status"] == "OK"
            for m in METHODS for s in SEEDS
        )
        if not has_success:
            continue
            
        print(f"\n### {dataset}")
        print("-" * 75)
        
        for method in METHODS:
            metrics_list = []
            for seed in SEEDS:
                r = results[dataset][method][seed]
                if r["status"] == "OK" and r["primary_metrics"]:
                    metrics_list.append(r["primary_metrics"])
            
            if not metrics_list:
                continue
            
            def fmt(values):
                if not values:
                    return "N/A"
                if len(values) == 1:
                    return f"{values[0]:.4f}"
                mean = sum(values) / len(values)
                if len(values) > 1:
                    std = (sum((v - mean)**2 for v in values) / len(values)) ** 0.5
                    return f"{mean:.4f}±{std:.4f}"
                return f"{mean:.4f}"
            
            aris = [m.get('ari') for m in metrics_list if m.get('ari') is not None]
            nmis = [m.get('nmi') for m in metrics_list if m.get('nmi') is not None]
            f1s = [m.get('f1_macro') for m in metrics_list if m.get('f1_macro') is not None]
            
            print(f"  {method:<35} | ARI: {fmt(aris):<15} | NMI: {fmt(nmis):<15} | F1: {fmt(f1s)}")
    
    # =========================================================================
    # PART 4: FAILURE ANALYSIS
    # =========================================================================
    print("\n\n" + "=" * 100)
    print("PART 4: FAILURE ANALYSIS")
    print("=" * 100)
    
    # Categorize failures by error type
    error_categories = defaultdict(list)
    
    for dataset in DATASETS:
        for method in METHODS:
            for seed in SEEDS:
                r = results[dataset][method][seed]
                if r["status"] == "FAILED":
                    error = r["error"] or "Unknown"
                    # Categorize
                    if "_index" in error and "reserved name" in error:
                        error_categories["_index reserved name (h5py/anndata)"].append((dataset, method, seed))
                    elif "CalledProcessError" in error:
                        error_categories["Subprocess failed (external script error)"].append((dataset, method, seed))
                    elif "Segmentation" in error or "rc=-11" in error:
                        error_categories["Segmentation fault (rc=-11)"].append((dataset, method, seed))
                    elif "CUDA" in error or "OOM" in error:
                        error_categories["CUDA out of memory"].append((dataset, method, seed))
                    else:
                        error_categories[f"Other: {error[:50]}"].append((dataset, method, seed))
    
    for error_type, failures in sorted(error_categories.items(), key=lambda x: -len(x[1])):
        print(f"\n[{len(failures)} failures] {error_type}")
        for dataset, method, seed in failures[:5]:  # Show first 5
            print(f"  - {dataset} / {method} / seed_{seed}")
        if len(failures) > 5:
            print(f"  ... and {len(failures) - 5} more")
    
    # =========================================================================
    # PART 5: OVERALL SUMMARY
    # =========================================================================
    print("\n\n" + "=" * 100)
    print("PART 5: OVERALL SUMMARY")
    print("=" * 100)
    
    total_runs = len(DATASETS) * len(METHODS) * len(SEEDS)
    successful = sum(1 for d in DATASETS for m in METHODS for s in SEEDS 
                     if results[d][m][s]["status"] == "OK")
    failed = sum(1 for d in DATASETS for m in METHODS for s in SEEDS 
                 if results[d][m][s]["status"] == "FAILED")
    
    print(f"\nTotal runs: {total_runs} ({len(DATASETS)} datasets × {len(METHODS)} methods × {len(SEEDS)} seeds)")
    print(f"  ✓ Successful: {successful} ({100*successful/total_runs:.1f}%)")
    print(f"  ✗ Failed:      {failed} ({100*failed/total_runs:.1f}%)")
    
    print("\n### Success rate by method:")
    for method in METHODS:
        ok = sum(1 for d in DATASETS for s in SEEDS if results[d][method][s]["status"] == "OK")
        total_method = len(DATASETS) * len(SEEDS)
        pct = 100 * ok / total_method
        bar = "█" * int(pct / 5) + "░" * (20 - int(pct / 5))
        print(f"  {method:<35} {bar} {ok:>2}/{total_method} ({pct:>3.0f}%)")
    
    print("\n### Success rate by dataset:")
    for dataset in DATASETS:
        ok = sum(1 for m in METHODS for s in SEEDS if results[dataset][m][s]["status"] == "OK")
        total_ds = len(METHODS) * len(SEEDS)
        pct = 100 * ok / total_ds
        bar = "█" * int(pct / 5) + "░" * (20 - int(pct / 5))
        print(f"  {dataset:<15} {bar} {ok:>2}/{total_ds} ({pct:>3.0f}%)")

if __name__ == "__main__":
    main()
