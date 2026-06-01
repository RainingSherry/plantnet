#!/usr/bin/env python3
"""Analyze benchmark results from PlantSPADE_LGCL_protocol."""

import json
import os
from pathlib import Path
from collections import defaultdict
import glob

BASE_DIR = Path("/home/luolie/biopipeline/dimension-reduction/plantnet/results/PlantSPADE_LGCL_protocol")

METHODS = [
    "traditional_pca", "phytocluster", "scvi", "scmae",
    "plantspade_lgcl_baseline", "plantspade_lgcl_support_attention"
]

DATASETS = [
    "SRP182008", "SRP235541", "SRP171040", "SRP309176",
    "SRP145013", "CRA002977_1", "SRP224648", "CRA007122"
]

SEEDS = [1, 2]  # The run used seeds 1,2

def check_metrics_exists(dataset, method, seed):
    """Check if metrics.json exists."""
    path = BASE_DIR / dataset / method / f"seed_{seed}" / "metrics.json"
    return path.exists(), path

def read_metrics(dataset, method, seed):
    """Read metrics.json if exists."""
    exists, path = check_metrics_exists(dataset, method, seed)
    if not exists:
        return None
    try:
        with open(path) as f:
            return json.load(f)
    except Exception as e:
        return {"error": str(e)}

def check_log_success(dataset, method, seed):
    """Check log file for success/failure status."""
    log_path = BASE_DIR / dataset / method / f"seed_{seed}" / "logs" / "run_suite.log"
    if not log_path.exists():
        return "NO_LOG", None
    
    try:
        with open(log_path) as f:
            content = f.read()
        
        # Check for common failure indicators
        if "Traceback (most recent call last)" in content:
            # Get last few lines of error
            lines = content.strip().split('\n')
            error_lines = [l for l in lines if 'Error' in l or 'Exception' in l or 'Traceback' in l]
            if error_lines:
                return "FAILED", '\n'.join(error_lines[-3:])
            return "FAILED", "Unknown error (check log)"
        
        if "CUDA out of memory" in content:
            return "OOM", "CUDA Out of Memory"
        
        if "done] ok" in content or "completed successfully" in content.lower():
            return "OK", None
        
        # Check for specific error patterns
        if "rc=1" in content or "rc=-11" in content:
            return "FAILED", "Non-zero exit code"
        
        # If we have metrics.json and no obvious error, consider it a success
        exists, _ = check_metrics_exists(dataset, method, seed)
        if exists:
            return "OK", None
        
        return "UNKNOWN", "Cannot determine status"
        
    except Exception as e:
        return "ERROR", str(e)

def main():
    print("=" * 80)
    print("PlantSPADE_LGCL Benchmark Results Summary")
    print("=" * 80)
    
    # Build status matrix
    results = defaultdict(lambda: defaultdict(dict))
    
    for dataset in DATASETS:
        for method in METHODS:
            for seed in SEEDS:
                status, error = check_log_success(dataset, method, seed)
                metrics = read_metrics(dataset, method, seed)
                results[dataset][method][seed] = {
                    "status": status,
                    "error": error,
                    "metrics": metrics
                }
    
    # Print summary by dataset
    print("\n" + "=" * 80)
    print("RESULTS BY DATASET")
    print("=" * 80)
    
    for dataset in DATASETS:
        print(f"\n### {dataset}")
        print("-" * 60)
        
        dataset_ok = 0
        dataset_failed = 0
        dataset_no_log = 0
        
        for method in METHODS:
            method_statuses = []
            for seed in SEEDS:
                r = results[dataset][method][seed]
                method_statuses.append(r["status"])
            
            ok_count = method_statuses.count("OK")
            fail_count = method_statuses.count("FAILED") + method_statuses.count("OOM")
            no_log_count = method_statuses.count("NO_LOG")
            
            if ok_count == len(SEEDS):
                status_icon = "✓"
                dataset_ok += 1
            elif ok_count == 0 and no_log_count == len(SEEDS):
                status_icon = "○"
                dataset_no_log += 1
            else:
                status_icon = "✗"
                dataset_failed += 1
            
            status_str = " | ".join(method_statuses)
            print(f"  {status_icon} {method:40s} [{status_str}]")
        
        total = len(METHODS)
        print(f"\n  Summary: {dataset_ok}/{total} methods complete, {dataset_failed}/{total} failed, {dataset_no_log}/{total} not started")
    
    # Print metrics comparison for successful runs
    print("\n" + "=" * 80)
    print("METRICS COMPARISON (Mean ± Std across seeds)")
    print("=" * 80)
    
    metrics_to_show = ["ari", "nmi", "f1", "accuracy"]
    
    for dataset in DATASETS:
        print(f"\n### {dataset}")
        
        header = f"{'Method':<35} | ARI            | NMI            | F1             | Accuracy"
        print(header)
        print("-" * 90)
        
        for method in METHODS:
            metrics_list = []
            for seed in SEEDS:
                r = results[dataset][method][seed]
                if r["status"] == "OK" and r["metrics"] and "error" not in r["metrics"]:
                    metrics_list.append(r["metrics"])
            
            if not metrics_list:
                print(f"{method:<35} | N/A")
                continue
            
            # Calculate mean and std for each metric
            row_parts = [f"{method:<35}"]
            has_data = False
            
            for metric in metrics_to_show:
                values = [m.get(metric) for m in metrics_list if m.get(metric) is not None]
                if values:
                    mean_val = sum(values) / len(values)
                    std_val = (sum((v - mean_val)**2 for v in values) / len(values)) ** 0.5 if len(values) > 1 else 0
                    row_parts.append(f"{mean_val:.4f}±{std_val:.4f}")
                    has_data = True
                else:
                    row_parts.append("N/A")
            
            if has_data:
                print(" | ".join(row_parts))
            else:
                print(f"{method:<35} | No comparable metrics")
    
    # Detailed failure analysis
    print("\n" + "=" * 80)
    print("FAILURE ANALYSIS")
    print("=" * 80)
    
    failures = []
    for dataset in DATASETS:
        for method in METHODS:
            for seed in SEEDS:
                r = results[dataset][method][seed]
                if r["status"] in ["FAILED", "OOM", "ERROR"]:
                    failures.append((dataset, method, seed, r["status"], r["error"]))
    
    if failures:
        for dataset, method, seed, status, error in failures:
            print(f"\n• {dataset} / {method} / seed_{seed}")
            print(f"  Status: {status}")
            if error:
                print(f"  Error: {error[:200]}")
    else:
        print("\nNo failures detected!")
    
    # Summary statistics
    print("\n" + "=" * 80)
    print("OVERALL SUMMARY")
    print("=" * 80)
    
    total_runs = len(DATASETS) * len(METHODS) * len(SEEDS)
    successful = sum(1 for d in DATASETS for m in METHODS for s in SEEDS 
                     if results[d][m][s]["status"] == "OK")
    failed = sum(1 for d in DATASETS for m in METHODS for s in SEEDS 
                 if results[d][m][s]["status"] in ["FAILED", "OOM", "ERROR"])
    no_log = sum(1 for d in DATASETS for m in METHODS for s in SEEDS 
                  if results[d][m][s]["status"] == "NO_LOG")
    
    print(f"\nTotal runs attempted: {total_runs}")
    print(f"  ✓ Successful: {successful} ({100*successful/total_runs:.1f}%)")
    print(f"  ✗ Failed:      {failed} ({100*failed/total_runs:.1f}%)")
    print(f"  ○ Not started:  {no_log} ({100*no_log/total_runs:.1f}%)")
    
    # Method-wise summary
    print("\n### Success rate by method:")
    for method in METHODS:
        ok = sum(1 for d in DATASETS for s in SEEDS if results[d][method][s]["status"] == "OK")
        total_method = len(DATASETS) * len(SEEDS)
        pct = 100 * ok / total_method
        bar = "█" * int(pct / 5) + "░" * (20 - int(pct / 5))
        print(f"  {method:<35} {bar} {ok}/{total_method} ({pct:.0f}%)")

if __name__ == "__main__":
    main()
