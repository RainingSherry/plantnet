#!/usr/bin/env python3
"""Analyze benchmark results from PlantSPADE_LGCL_protocol - Fixed version."""

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

DATASET_INFO = {
    "SRP182008": ("A. thaliana", "Root tip"),
    "SRP235541": ("A. thaliana", "Root tip"),
    "SRP171040": ("A. thaliana", "Root tip"),
    "SRP309176": ("O. sativa", "Root tip"),
    "SRP145013": ("Z. mays", "Root tip"),
    "CRA002977_1": ("A. thaliana", "Leaf"),
    "SRP224648": ("Z. mays", "Leaf"),
    "CRA007122": ("G. max", "Nodule"),
}

SEEDS = [1, 2]

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

def read_log_error(dataset, method, seed):
    """Read log file for error information."""
    log_path = BASE_DIR / dataset / method / f"seed_{seed}" / "logs" / "run_suite.log"
    if not log_path.exists():
        return "NO_LOG"
    
    try:
        with open(log_path) as f:
            content = f.read()
        
        if "Traceback (most recent call last)" in content:
            lines = content.strip().split('\n')
            error_lines = [l for l in lines if 'Error' in l or 'Exception' in l or 'ValueError' in l or 'Traceback' in l or 'CalledProcessError' in l]
            return '\n'.join(error_lines[-5:]) if error_lines else "Unknown error"
        
        if "_index" in content and "reserved name" in content:
            return "ValueError: '_index' is a reserved name for dataframe columns"
        
        if "CUDA out of memory" in content:
            return "CUDA Out of Memory"
        
        if "rc=-11" in content:
            return "Segmentation fault (rc=-11)"
        
        return None  # No error found
        
    except Exception as e:
        return str(e)

def get_status(dataset, method, seed):
    """Determine run status based on metrics.json and log."""
    metrics = read_metrics(dataset, method, seed)
    if metrics and "error" not in metrics:
        return "OK", metrics, None
    
    error = read_log_error(dataset, method, seed)
    if error:
        if "_index" in str(error):
            return "FAILED_H5AD", metrics, error
        return "FAILED", metrics, error
    
    return "NOT_STARTED", metrics, None

def get_primary_metrics(metrics):
    """Extract primary metrics from metrics.json."""
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
    print("=" * 110)
    print("PlantSPADE_LGCL Benchmark Results Summary")
    print("Command: --datasets all --seeds 1,2 --gpus 1,2,3,4,5,6 --methods all --keep_going")
    print("=" * 110)
    
    # Build results matrix
    results = defaultdict(lambda: defaultdict(dict))
    
    for dataset in DATASETS:
        for method in METHODS:
            for seed in SEEDS:
                status, metrics, error = get_status(dataset, method, seed)
                pm = get_primary_metrics(metrics)
                results[dataset][method][seed] = {
                    "status": status,
                    "metrics": metrics,
                    "primary_metrics": pm,
                    "error": error
                }
    
    # =========================================================================
    # PART 1: STATUS OVERVIEW TABLE
    # =========================================================================
    print("\n" + "=" * 110)
    print("PART 1: RUN STATUS OVERVIEW")
    print("=" * 110)
    print("Legend: ✓ = Success, ✗ = Failed (H5AD error), ? = Failed (other), ○ = Not started")
    
    print(f"\n{'Dataset':<15} | {'Species':<12} | {'Tissue':<8} | K | " + 
          " | ".join([f"{m:^10}" for m in METHODS]))
    print("-" * 110)
    
    for dataset in DATASETS:
        species, tissue = DATASET_INFO[dataset]
        statuses = []
        for method in METHODS:
            s_list = [results[dataset][method][s]["status"] for s in SEEDS]
            
            # Determine combined status
            if all(s == "OK" for s in s_list):
                statuses.append("\033[92m✓\033[0m")
            elif all(s in ["FAILED_H5AD", "FAILED"] for s in s_list):
                statuses.append("\033[91m✗\033[0m")
            elif all(s == "NOT_STARTED" for s in s_list):
                statuses.append("\033[90m○\033[0m")
            elif any(s == "OK" for s in s_list):
                statuses.append("\033[93m½\033[0m")
            else:
                statuses.append("\033[93m?\033[0m")
        
        print(f"{dataset:<15} | {species:<12} | {tissue:<8} | {EXPECTED_CLUSTERS[dataset]:<2} | " + 
              " | ".join(statuses))
    
    # =========================================================================
    # PART 2: SUCCESSFUL RUNS DETAILS
    # =========================================================================
    print("\n\n" + "=" * 110)
    print("PART 2: SUCCESSFUL RUNS - DETAILED METRICS")
    print("=" * 110)
    print("Metrics: ARI (Adjusted Rand Index), NMI (Normalized Mutual Information), F1 (F1-macro), ACC (Accuracy)")
    print("Best clustering method shown: kmeans with known K")
    
    for dataset in DATASETS:
        has_success = any(
            results[dataset][m][s]["status"] == "OK"
            for m in METHODS for s in SEEDS
        )
        if not has_success:
            continue
            
        print(f"\n### {dataset} ({DATASET_INFO[dataset][0]} - {DATASET_INFO[dataset][1]}, Expected K={EXPECTED_CLUSTERS[dataset]})")
        print("-" * 105)
        print(f"{'Method':<35} | {'Seed':>4} | {'ARI':>8} | {'NMI':>8} | {'F1':>8} | {'ACC':>8} | {'Status':<10}")
        print("-" * 105)
        
        for method in METHODS:
            for seed in SEEDS:
                r = results[dataset][method][seed]
                pm = r["primary_metrics"]
                
                if r["status"] == "OK" and pm:
                    ari = f"{pm.get('ari', 0):.4f}" if pm.get('ari') is not None else "N/A"
                    nmi = f"{pm.get('nmi', 0):.4f}" if pm.get('nmi') is not None else "N/A"
                    f1 = f"{pm.get('f1_macro', 0):.4f}" if pm.get('f1_macro') is not None else "N/A"
                    acc = f"{pm.get('acc', 0):.4f}" if pm.get('acc') is not None else "N/A"
                    print(f"{method:<35} | {seed:>4} | {ari:>8} | {nmi:>8} | {f1:>8} | {acc:>8} | \033[92mOK\033[0m")
                elif r["status"] == "OK":
                    print(f"{method:<35} | {seed:>4} | \033[90m(no comparable metrics)\033[0m")
                else:
                    status_str = r["status"].replace("_", " ")
                    print(f"{method:<35} | {seed:>4} | \033[90m{status_str}\033[0m")
    
    # =========================================================================
    # PART 3: AGGREGATED COMPARISON (Mean ± Std)
    # =========================================================================
    print("\n\n" + "=" * 110)
    print("PART 3: AGGREGATED METRICS BY METHOD (Mean ± Std across seeds)")
    print("=" * 110)
    
    for dataset in DATASETS:
        has_success = any(
            results[dataset][m][s]["status"] == "OK"
            for m in METHODS for s in SEEDS
        )
        if not has_success:
            continue
            
        print(f"\n### {dataset}")
        
        def fmt_stats(values):
            if not values:
                return "N/A"
            if len(values) == 1:
                return f"{values[0]:.4f}"
            mean = sum(values) / len(values)
            std = (sum((v - mean)**2 for v in values) / len(values)) ** 0.5
            return f"{mean:.4f}±{std:.4f}"
        
        rows = []
        for method in METHODS:
            metrics_list = []
            for seed in SEEDS:
                r = results[dataset][method][seed]
                if r["status"] == "OK" and r["primary_metrics"]:
                    metrics_list.append(r["primary_metrics"])
            
            if not metrics_list:
                continue
            
            aris = [m.get('ari') for m in metrics_list if m.get('ari') is not None]
            nmis = [m.get('nmi') for m in metrics_list if m.get('nmi') is not None]
            f1s = [m.get('f1_macro') for m in metrics_list if m.get('f1_macro') is not None]
            
            rows.append((method, fmt_stats(aris), fmt_stats(nmis), fmt_stats(f1s)))
        
        if rows:
            print(f"{'Method':<35} | {'ARI':^15} | {'NMI':^15} | {'F1-macro':^15}")
            print("-" * 85)
            for method, ari, nmi, f1 in rows:
                print(f"{method:<35} | {ari:>15} | {nmi:>15} | {f1:>15}")
    
    # =========================================================================
    # PART 4: CROSS-DATASET METHOD RANKING
    # =========================================================================
    print("\n\n" + "=" * 110)
    print("PART 4: METHOD RANKING (Average ARI across successful datasets)")
    print("=" * 110)
    
    method_aris = defaultdict(list)
    for dataset in DATASETS:
        for method in METHODS:
            for seed in SEEDS:
                r = results[dataset][method][seed]
                if r["status"] == "OK" and r["primary_metrics"] and r["primary_metrics"].get('ari') is not None:
                    method_aris[method].append(r["primary_metrics"]['ari'])
    
    method_avg = []
    for method, aris in method_aris.items():
        if aris:
            avg = sum(aris) / len(aris)
            std = (sum((v - avg)**2 for v in aris) / len(aris)) ** 0.5
            method_avg.append((method, avg, std, len(aris)))
    
    method_avg.sort(key=lambda x: -x[1])
    
    print(f"{'Rank':<6} | {'Method':<40} | {'Mean ARI':>12} | {'Std':>8} | {'N runs':>8}")
    print("-" * 80)
    for i, (method, avg, std, n) in enumerate(method_avg, 1):
        print(f"{i:<6} | {method:<40} | {avg:>12.4f} | {std:>8.4f} | {n:>8}")
    
    # =========================================================================
    # PART 5: FAILURE ANALYSIS
    # =========================================================================
    print("\n\n" + "=" * 110)
    print("PART 5: FAILURE ANALYSIS")
    print("=" * 110)
    
    # Categorize failures
    error_categories = defaultdict(list)
    
    for dataset in DATASETS:
        for method in METHODS:
            for seed in SEEDS:
                r = results[dataset][method][seed]
                if r["status"] in ["FAILED", "FAILED_H5AD"]:
                    error = r["error"] or "Unknown"
                    
                    if "_index" in str(error):
                        error_categories["H5AD '_index' reserved name error"].append((dataset, method, seed))
                    elif "CUDA" in str(error) or "OOM" in str(error):
                        error_categories["CUDA Out of Memory"].append((dataset, method, seed))
                    elif "Segmentation" in str(error):
                        error_categories["Segmentation fault"].append((dataset, method, seed))
                    else:
                        error_categories[f"Other: {str(error)[:60]}"].append((dataset, method, seed))
    
    for error_type, failures in sorted(error_categories.items(), key=lambda x: -len(x[1])):
        print(f"\n[{len(failures)} failures] {error_type}")
        shown = 0
        for dataset, method, seed in failures:
            print(f"  - {dataset} / {method} / seed_{seed}")
            shown += 1
            if shown >= 3:
                break
        if len(failures) > 3:
            print(f"  ... and {len(failures) - 3} more")
    
    # =========================================================================
    # PART 6: OVERALL SUMMARY
    # =========================================================================
    print("\n\n" + "=" * 110)
    print("PART 6: OVERALL SUMMARY")
    print("=" * 110)
    
    total_runs = len(DATASETS) * len(METHODS) * len(SEEDS)
    successful = sum(1 for d in DATASETS for m in METHODS for s in SEEDS 
                     if results[d][m][s]["status"] == "OK")
    failed = sum(1 for d in DATASETS for m in METHODS for s in SEEDS 
                 if results[d][m][s]["status"] in ["FAILED", "FAILED_H5AD"])
    not_started = sum(1 for d in DATASETS for m in METHODS for s in SEEDS 
                      if results[d][m][s]["status"] == "NOT_STARTED")
    
    print(f"\nTotal runs: {total_runs} (8 datasets × 6 methods × 2 seeds)")
    print(f"  ✓ Successful:  {successful:>3} ({100*successful/total_runs:>5.1f}%)")
    print(f"  ✗ Failed:       {failed:>3} ({100*failed/total_runs:>5.1f}%)")
    print(f"  ○ Not started: {not_started:>3} ({100*not_started/total_runs:>5.1f}%)")
    
    print("\n### Success rate by method:")
    for method in METHODS:
        ok = sum(1 for d in DATASETS for s in SEEDS if results[d][method][s]["status"] == "OK")
        total_method = len(DATASETS) * len(SEEDS)
        pct = 100 * ok / total_method
        bar = "█" * max(0, int(pct / 5)) + "░" * max(0, 20 - int(pct / 5))
        print(f"  {method:<35} {bar} {ok:>2}/{total_method} ({pct:>3.0f}%)")
    
    print("\n### Success rate by dataset:")
    for dataset in DATASETS:
        ok = sum(1 for m in METHODS for s in SEEDS if results[dataset][m][s]["status"] == "OK")
        total_ds = len(METHODS) * len(SEEDS)
        pct = 100 * ok / total_ds
        bar = "█" * max(0, int(pct / 5)) + "░" * max(0, 20 - int(pct / 5))
        print(f"  {dataset:<15} {bar} {ok:>2}/{total_ds} ({pct:>3.0f}%)")

if __name__ == "__main__":
    main()
