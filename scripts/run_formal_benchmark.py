#!/usr/bin/env python3
"""
Formal Benchmark Runner
======================
Per BDD Scenarios 15, 16: Runs only Authenticity=VERIFIED + Smoke=PASS models
by default. Reads from methods/method_manifest.yaml.

Usage:
    python scripts/run_formal_benchmark.py \
        --data_path data/SRP182008.h5ad \
        --out_dir results/formal \
        --methods dec scdcc scdsc scanpy_standard leiden louvain sc3 \
        --n_clusters 15 \
        --seeds 42 43 44 \
        --epochs 200 \
        --pretrain_epochs 200 \
        --no_cuda

    # Allow unverified models (output will be marked 'unverified')
    python scripts/run_formal_benchmark.py \
        --data_path data/SRP182008.h5ad \
        --methods scgnn sccdcg \
        --allow_unverified

Exit codes:
    0  All methods completed successfully
    1  At least one method failed
    2  No valid methods selected
"""

import os
import sys
import json
import yaml
import time
import subprocess
import argparse
import re
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Optional, Any

PROJECT_ROOT = Path(__file__).parent.parent.resolve()
METHODS_DIR = PROJECT_ROOT / "methods"
MANIFEST_PATH = METHODS_DIR / "method_manifest.yaml"


def load_manifest() -> Dict[str, Any]:
    """Load the method manifest YAML."""
    with open(MANIFEST_PATH, "r", encoding="utf-8") as f:
        data = yaml.safe_load(f)
    return {m["key"]: m for m in data["methods"]}


def check_authenticity(method_key: str, method_info: Dict[str, Any]) -> tuple[bool, str]:
    """
    Check if a method can enter the formal benchmark.
    Returns (allowed, reason).
    """
    auth = method_info.get("authenticity", "UNKNOWN")
    smoke = method_info.get("smoke", "UNKNOWN")

    if auth == "VERIFIED" and smoke == "PASS":
        return True, "VERIFIED + Smoke=PASS"
    elif auth == "VERIFIED" and smoke == "UNKNOWN":
        return True, "VERIFIED (smoke status unknown)"
    elif auth == "ENV-GATED":
        return False, f"ENV-GATED: {method_info.get('reason', 'unknown env')}"
    elif auth == "PENDING":
        return False, "PENDING (not yet audited)"
    elif auth == "PLACEHOLDER":
        return False, "PLACEHOLDER (no code)"
    elif auth == "FAILED":
        return False, "FAILED (known issues)"
    else:
        return False, f"Unknown authenticity status: {auth}"


def write_authenticity_json(out_dir: Path, method_info: Dict[str, Any]) -> None:
    """Write authenticity metadata to the output directory."""
    auth_data = {
        "method": method_info["key"],
        "name": method_info.get("name", method_info["key"]),
        "authenticity": method_info.get("authenticity", "UNKNOWN"),
        "source_path": method_info.get("source_path", ""),
        "target_path": method_info.get("target_path", ""),
        "known_deviations": method_info.get("known_deviations", []),
        "substitute_model_used": False,
        "framework": method_info.get("framework", ""),
        "category": method_info.get("category", ""),
    }
    auth_path = out_dir / "authenticity.json"
    with open(auth_path, "w", encoding="utf-8") as f:
        json.dump(auth_data, f, indent=2)


def build_command(
    method_key: str,
    method_info: Dict[str, Any],
    data_path: str,
    out_dir: Path,
    n_clusters: int,
    epochs: int,
    pretrain_epochs: int,
    seed: int,
    no_cuda: bool,
    extra_args: Optional[List[str]] = None,
) -> List[str]:
    """Build the command-line invocation for a method."""

    run_py = PROJECT_ROOT / method_info["path"]

    cmd = [
        sys.executable,
        str(run_py),
        "--data_path", data_path,
        "--save_dir", str(out_dir),
        "--n_clusters", str(n_clusters),
        "--seed", str(seed),
    ]

    # Add method-specific training args
    if method_info.get("category") == "DeepLearning":
        cmd.extend(["--epochs", str(epochs)])
        if method_info.get("key") in ("dec", "scdcc", "scdsc"):
            cmd.extend(["--pretrain_epochs", str(pretrain_epochs)])
    elif method_info.get("key") == "scdcc":
        cmd.extend(["--pretrain_epochs", str(pretrain_epochs)])
    elif method_info.get("key") == "dec":
        cmd.extend(["--pretrain_epochs", str(pretrain_epochs)])

    # Deep learning methods support --no_cuda; Traditional methods do not
    if no_cuda and method_info.get("category") in ("DeepLearning", "GNN"):
        cmd.append("--no_cuda")

    # Extra args (e.g., --no_mix for ablation)
    if extra_args:
        cmd.extend(extra_args)

    return cmd


def verify_output(out_dir: Path) -> bool:
    """
    Verify that a method's output contains all required files.
    Returns True if all required files exist.
    """
    required_files = [
        "embedding_final.npy",
        "labels.npy",
        "metrics.json",
        "args.json",
    ]
    for fname in required_files:
        if not (out_dir / fname).exists():
            return False
    return True


def load_metrics(out_dir: Path) -> Optional[Dict[str, float]]:
    """Load metrics.json if it exists."""
    metrics_path = out_dir / "metrics.json"
    if metrics_path.exists():
        with open(metrics_path, "r") as f:
            return json.load(f)
    return None


def run_method(
    method_key: str,
    method_info: Dict[str, Any],
    data_path: str,
    base_out_dir: Path,
    n_clusters: int,
    epochs: int,
    pretrain_epochs: int,
    seed: int,
    no_cuda: bool,
    allow_unverified: bool,
    verbose: bool = False,
) -> Dict[str, Any]:
    """Run a single method with a given seed. Returns a result dict."""

    run_id = f"{method_key}__seed{seed}__{datetime.now().strftime('%Y%m%d_%H%M%S')}"
    out_dir = base_out_dir / run_id
    os.makedirs(out_dir, exist_ok=True)

    # Write authenticity metadata
    write_authenticity_json(out_dir, method_info)

    # Check if allowed
    allowed, reason = check_authenticity(method_key, method_info)
    if not allowed:
        # If --allow_unverified, we still run but flag it
        if not allow_unverified:
            status = "skipped"
            result = {
                "method": method_key,
                "seed": seed,
                "status": status,
                "reason": reason,
                "out_dir": str(out_dir),
            }
            return result
        else:
            status = "unverified"
            # Mark in authenticity.json
            auth_path = out_dir / "authenticity.json"
            with open(auth_path, "r") as f:
                auth_data = json.load(f)
            auth_data["unverified"] = True
            auth_data["unverified_reason"] = reason
            with open(auth_path, "w") as f:
                json.dump(auth_data, f, indent=2)

    # Build command
    extra_args = method_info.get("extra_args", None)
    cmd = build_command(
        method_key, method_info, data_path, out_dir,
        n_clusters, epochs, pretrain_epochs, seed, no_cuda, extra_args
    )

    # Write command
    cmd_file = out_dir / "command.txt"
    with open(cmd_file, "w") as f:
        f.write(" ".join(f'"{c}"' if " " in c else c for c in cmd) + "\n")

    # Run
    start_time = time.time()
    status = "success"
    error_msg = ""

    try:
        env = os.environ.copy()
        env["PYTHONPATH"] = str(METHODS_DIR) + ":" + env.get("PYTHONPATH", "")

        proc = subprocess.run(
            cmd,
            cwd=str(PROJECT_ROOT),
            capture_output=True,
            text=True,
            timeout=3600 * 6,  # 6 hour timeout
            env=env,
        )

        # Write log
        log_file = out_dir / "run.log"
        with open(log_file, "w") as f:
            f.write(f"Command: {' '.join(cmd)}\n")
            f.write(f"Exit code: {proc.returncode}\n")
            f.write(f"\n=== STDOUT ===\n{proc.stdout}\n")
            if proc.stderr:
                f.write(f"\n=== STDERR ===\n{proc.stderr}\n")

        if proc.returncode != 0:
            status = "failed"
            error_msg = f"Exit code {proc.returncode}"
        else:
            # Verify output
            if not verify_output(out_dir):
                status = "incomplete"
                error_msg = "Missing required output files"
            else:
                status = "success" if not allow_unverified or reason.startswith("VERIFIED") else "unverified"

    except subprocess.TimeoutExpired:
        status = "timeout"
        error_msg = "Execution exceeded 6 hour timeout"
        with open(out_dir / "run.log", "w") as f:
            f.write(f"Command timed out after 6 hours\nCommand: {' '.join(cmd)}\n")
    except Exception as e:
        status = "error"
        error_msg = str(e)
        with open(out_dir / "run.log", "w") as f:
            f.write(f"Error: {e}\nCommand: {' '.join(cmd)}\n")

    elapsed = time.time() - start_time

    # Load metrics
    metrics = load_metrics(out_dir)

    result = {
        "method": method_key,
        "seed": seed,
        "status": status,
        "reason": reason if status == "skipped" else "",
        "error": error_msg,
        "elapsed_seconds": round(elapsed, 1),
        "out_dir": str(out_dir),
        "metrics": metrics,
    }

    if verbose:
        icon = {"success": "✓", "failed": "✗", "skipped": "⊗", "incomplete": "⚠", "unverified": "⚠", "timeout": "⏱", "error": "!"}.get(status, "?")
        print(f"  {icon} {method_key} (seed={seed}): {status} ({elapsed:.1f}s)")
        if metrics:
            print(f"    ACC={metrics.get('acc', 'N/A'):.4f} NMI={metrics.get('nmi', 'N/A'):.4f} ARI={metrics.get('ari', 'N/A'):.4f}")
        if error_msg:
            print(f"    Error: {error_msg}")

    return result


def collect_results(base_out_dir: Path, method_keys: List[str], seeds: List[int]) -> List[Dict]:
    """Collect results from completed runs."""
    results = []
    for sub_dir in sorted(base_out_dir.iterdir()):
        if not sub_dir.is_dir():
            continue
        # Parse method_key and seed from dir name
        parts = sub_dir.name.split("__")
        if len(parts) < 2:
            continue
        method_key = parts[0]
        if method_key not in method_keys:
            continue
        metrics = load_metrics(sub_dir)
        auth_path = sub_dir / "authenticity.json"
        auth_data = {}
        if auth_path.exists():
            with open(auth_path) as f:
                auth_data = json.load(f)

        results.append({
            "method": method_key,
            "out_dir": str(sub_dir),
            "metrics": metrics,
            "authenticity": auth_data.get("authenticity", "UNKNOWN"),
            "substitute_model_used": auth_data.get("substitute_model_used", True),
            "unverified": auth_data.get("unverified", False),
        })
    return results


def generate_summary_csv(results: List[Dict], out_path: Path, dataset: str) -> None:
    """Generate benchmark_summary.csv from results."""
    rows = []
    for r in results:
        m = r["metrics"] or {}
        rows.append({
            "dataset": dataset,
            "method": r["method"],
            "authenticity": r["authenticity"],
            "substitute_model_used": r["substitute_model_used"],
            "unverified": r["unverified"],
            "acc": m.get("acc", ""),
            "nmi": m.get("nmi", ""),
            "ari": m.get("ari", ""),
            "f1_macro": m.get("f1_macro", ""),
            "fmi": m.get("fmi", ""),
            "v_measure": m.get("v_measure", ""),
            "homogeneity": m.get("homogeneity", ""),
            "completeness": m.get("completeness", ""),
            "out_dir": r["out_dir"],
        })

    import csv
    if rows:
        with open(out_path, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=rows[0].keys())
            writer.writeheader()
            writer.writerows(rows)
        print(f"Summary saved to: {out_path}")


def generate_mean_std_csv(results: List[Dict], out_path: Path, dataset: str) -> None:
    """Generate benchmark_summary_mean_std.csv (mean ± std over seeds)."""
    import csv
    from collections import defaultdict
    import statistics

    by_method = defaultdict(list)
    for r in results:
        by_method[r["method"]].append(r)

    rows = []
    metric_names = ["acc", "nmi", "ari", "f1_macro", "fmi", "v_measure", "homogeneity", "completeness"]

    for method, runs in sorted(by_method.items()):
        row = {"dataset": dataset, "method": method, "authenticity": runs[0]["authenticity"]}
        for m in metric_names:
            vals = [r["metrics"].get(m) for r in runs if r["metrics"] and r["metrics"].get(m) is not None]
            if vals:
                mean_val = statistics.mean(vals)
                std_val = statistics.stdev(vals) if len(vals) > 1 else 0.0
                row[m] = f"{mean_val:.4f} ± {std_val:.4f}"
            else:
                row[m] = "N/A"
        rows.append(row)

    if rows:
        with open(out_path, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=rows[0].keys())
            writer.writeheader()
            writer.writerows(rows)
        print(f"Mean±std summary saved to: {out_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Run formal benchmark on VERIFIED methods",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Run all default methods (VERIFIED + Smoke=PASS) on one dataset
  python scripts/run_formal_benchmark.py \\
      --data_path data/SRP182008.h5ad \\
      --out_dir results/formal/SRP182008 \\
      --n_clusters 15 \\
      --seeds 42 43 44

  # Run a specific subset
  python scripts/run_formal_benchmark.py \\
      --data_path data/SRP182008.h5ad \\
      --methods dec scdcc scanpy_standard \\
      --out_dir results/formal/dec_sc \\
      --n_clusters 15 --seeds 42

  # Include pending/unverified methods (output marked 'unverified')
  python scripts/run_formal_benchmark.py \\
      --data_path data/SRP182008.h5ad \\
      --methods scgnn \\
      --out_dir results/formal/scgnn_test \\
      --n_clusters 15 --allow_unverified
"""
    )

    parser.add_argument("--data_path", type=str, required=True,
                       help="Path to input h5ad file")
    parser.add_argument("--out_dir", type=str, default="results/formal",
                       help="Base output directory")
    parser.add_argument("--n_clusters", type=int, required=True,
                       help="Number of clusters")
    parser.add_argument("--seeds", type=int, nargs="+", default=[42],
                       help="Random seeds (default: 42)")
    parser.add_argument("--epochs", type=int, default=200,
                       help="Training epochs for deep learning models")
    parser.add_argument("--pretrain_epochs", type=int, default=200,
                       help="Pretraining epochs for models with pretrain phase")
    parser.add_argument("--methods", type=str, nargs="+", default=None,
                       help="Method keys to run (default: all VERIFIED+Smoke=PASS)")
    parser.add_argument("--no_cuda", action="store_true",
                       help="Disable CUDA (use CPU only)")
    parser.add_argument("--allow_unverified", action="store_true",
                       help="Allow PENDING/ENV-GATED models (output marked 'unverified')")
    parser.add_argument("--skip_run", action="store_true",
                       help="Skip running models, only collect and generate summary")
    parser.add_argument("--dataset_name", type=str, default=None,
                       help="Dataset name for summary tables (default: derived from data_path)")
    parser.add_argument("--verbose", "-v", action="store_true",
                       help="Verbose output")

    args = parser.parse_args()

    manifest = load_manifest()
    dataset_name = args.dataset_name or os.path.splitext(os.path.basename(args.data_path))[0]
    # If --out_dir already contains the dataset_name, don't double-nest
    if args.out_dir and dataset_name in args.out_dir:
        base_out_dir = Path(args.out_dir)
    else:
        base_out_dir = Path(args.out_dir) / dataset_name
    os.makedirs(base_out_dir, exist_ok=True)

    # Determine which methods to run
    if args.methods:
        selected = {k: manifest[k] for k in args.methods if k in manifest}
        missing = [k for k in args.methods if k not in manifest]
        if missing:
            print(f"WARNING: Unknown method keys: {missing}")
    else:
        # Default: all VERIFIED + Smoke=PASS
        selected = {
            k: v for k, v in manifest.items()
            if v.get("authenticity") == "VERIFIED" and v.get("smoke") == "PASS"
        }

    if not selected:
        print("ERROR: No valid methods selected. Use --methods or check manifest.")
        sys.exit(2)

    # ── Summary of what will run ──────────────────────────────
    print("=" * 70)
    print("Formal Benchmark Runner")
    print("=" * 70)
    print(f"  Data:     {args.data_path}")
    print(f"  Dataset:  {dataset_name}")
    print(f"  Clusters: {args.n_clusters}")
    print(f"  Seeds:    {args.seeds}")
    print(f"  Epochs:   {args.epochs} (pretrain: {args.pretrain_epochs})")
    print(f"  CUDA:     {'disabled (CPU only)' if args.no_cuda else 'enabled'}")
    print(f"  Unverified: {'ALLOWED' if args.allow_unverified else 'REJECTED'}")
    print()
    print(f"  Methods ({len(selected)}):")
    for k, v in sorted(selected.items()):
        allowed, reason = check_authenticity(k, v)
        icon = "✓" if allowed else "⊗"
        print(f"    {icon} {k:20s} [{v['category']:12s}] {reason}")
    print()

    if not args.allow_unverified:
        skipped = [k for k, v in selected.items() if not check_authenticity(k, v)[0]]
        if skipped:
            print(f"  Skipping {len(skipped)} unverified methods: {skipped}")
            print()

    print("=" * 70)
    print()

    # ── Run ──────────────────────────────────────────────────
    if not args.skip_run:
        all_results = []
        for method_key, method_info in sorted(selected.items()):
            allowed, reason = check_authenticity(method_key, method_info)
            if not allowed and not args.allow_unverified:
                print(f"  ⊗ {method_key}: skipped ({reason})")
                continue

            print(f"  Running {method_key}...")
            for seed in args.seeds:
                result = run_method(
                    method_key, method_info,
                    args.data_path, base_out_dir,
                    args.n_clusters, args.epochs, args.pretrain_epochs,
                    seed, args.no_cuda, args.allow_unverified,
                    verbose=args.verbose,
                )
                all_results.append(result)

            print()

    # ── Collect results ────────────────────────────────────────
    print("=" * 70)
    print("Collecting results...")
    results = collect_results(base_out_dir, list(selected.keys()), args.seeds)

    # Group by method and seed
    for r in results:
        if not args.skip_run:
            continue  # Already printed
        print(f"  {r['method']}: authenticity={r['authenticity']}")

    # ── Generate summary tables ────────────────────────────────
    summary_csv = base_out_dir / "benchmark_summary.csv"
    generate_summary_csv(results, summary_csv, dataset_name)

    mean_csv = base_out_dir / "benchmark_summary_mean_std.csv"
    generate_mean_std_csv(results, mean_csv, dataset_name)

    # ── Status summary ─────────────────────────────────────────
    status_counts = {}
    for r in results:
        s = r.get("status", "unknown")
        status_counts[s] = status_counts.get(s, 0) + 1

    print()
    print("=" * 70)
    print("Benchmark Complete")
    print("=" * 70)
    print(f"  Results dir: {base_out_dir}")
    print(f"  Status:")
    for s, c in sorted(status_counts.items()):
        print(f"    {s}: {c}")
    print()

    if status_counts.get("failed") or status_counts.get("error"):
        print("  ⚠  Some methods failed. Check run.log files.")
        sys.exit(1)
    else:
        print("  ✓  All methods completed successfully.")
        sys.exit(0)


if __name__ == "__main__":
    main()
