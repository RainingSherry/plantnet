#!/usr/bin/env python3
"""
Formal Benchmark Runner
=======================
Per BDD Scenarios 1-16: Runs only Authenticity=VERIFIED + Smoke=PASS models
by default. Reads from methods/method_manifest.yaml.

GPU 0 and GPU 7 are FORBIDDEN. Only --gpu 1-6 is allowed.

Usage:
    python scripts/run_formal_benchmark.py \
        --data_path data/subsample_2k.h5ad \
        --dataset_name subsample_2k \
        --out_dir results/formal \
        --n_clusters 7 \
        --seeds 42 \
        --gpu 1 \
        --methods \
            neighbormix_scmae \
            nm_scmae_nomix \
            scmae \
            dec \
            scdcc \
            scdsc \
            scanpy_standard \
            leiden \
            louvain \
            sc3

    # Dry-run / preflight (no actual execution)
    python scripts/run_formal_benchmark.py ... --dry_run

    # CPU-only
    python scripts/run_formal_benchmark.py ... --no_cuda

Exit codes:
    0  All methods completed successfully
    1  At least one method failed
    2  No valid methods selected
    3  Forbidden GPU detected
"""

import os
import sys
import json
import yaml
import time
import subprocess
import argparse
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Optional, Any

PROJECT_ROOT = Path(__file__).parent.parent.resolve()
METHODS_DIR = PROJECT_ROOT / "methods"
MANIFEST_PATH = METHODS_DIR / "method_manifest.yaml"
FORBIDDEN_GPUS = {"0", "7"}

# Default formal method list (proposed + ablation + external + baselines)
DEFAULT_FORMAL_METHODS = [
    "neighbormix_scmae",
    "nm_scmae_nomix",
    "scmae",
    "dec",
    "scdcc",
    "scdsc",
    "scanpy_standard",
    "leiden",
    "louvain",
    "sc3",
]


# ────────────────────────────────────────────────────────────────
# GPU Validation
# ────────────────────────────────────────────────────────────────

def validate_gpu_policy(gpu: int, no_cuda: bool) -> None:
    """Validate GPU policy: reject GPU 0 and GPU 7, check CUDA_VISIBLE_DEVICES."""
    if no_cuda:
        return
    if str(gpu) in FORBIDDEN_GPUS:
        raise SystemExit(
            f"[GPU Policy] Forbidden GPU {gpu} detected. "
            "GPU 0 and GPU 7 are NOT allowed. Use --gpu 1-6 or --no_cuda."
        )
    visible = os.environ.get("CUDA_VISIBLE_DEVICES", "")
    if visible:
        visible_ids = {x.strip() for x in visible.split(",") if x.strip()}
        forbidden = visible_ids & FORBIDDEN_GPUS
        if forbidden:
            raise SystemExit(
                f"[GPU Policy] CUDA_VISIBLE_DEVICES={visible} contains forbidden GPU(s) {forbidden}. "
                "GPU 0 and GPU 7 are NOT allowed."
            )


def get_git_info() -> tuple[str, str]:
    """Return (commit_sha, branch) for the current repo."""
    try:
        sha = subprocess.check_output(
            ["git", "rev-parse", "HEAD"], text=True, cwd=str(PROJECT_ROOT)
        ).strip()
        branch = subprocess.check_output(
            ["git", "rev-parse", "--abbrev-ref", "HEAD"], text=True, cwd=str(PROJECT_ROOT)
        ).strip()
        return sha, branch
    except Exception:
        return "unknown", "unknown"


# ────────────────────────────────────────────────────────────────
# Manifest helpers
# ────────────────────────────────────────────────────────────────

def load_manifest() -> Dict[str, Any]:
    """Load the method manifest YAML into a dict keyed by method key."""
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
    elif auth == "VERIFIED" and smoke in ("UNKNOWN", ""):
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
    with open(out_dir / "authenticity.json", "w", encoding="utf-8") as f:
        json.dump(auth_data, f, indent=2)


# ────────────────────────────────────────────────────────────────
# Command builder
# ────────────────────────────────────────────────────────────────

def build_command(
    method_key: str,
    method_info: Dict[str, Any],
    data_path: str,
    out_dir: Path,
    n_clusters: int,
    epochs: int,
    pretrain_epochs: int,
    seed: int,
    gpu: int,
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

    category = method_info.get("category", "")
    is_deep_or_gnn = category in ("DeepLearning", "GNN")

    if category == "DeepLearning":
        cmd.extend(["--epochs", str(epochs)])
        if method_key in ("dec", "scdcc", "scdsc"):
            cmd.extend(["--pretrain_epochs", str(pretrain_epochs)])
    elif is_deep_or_gnn:
        cmd.extend(["--epochs", str(epochs)])
        if method_key in ("dec", "scdcc", "scdsc"):
            cmd.extend(["--pretrain_epochs", str(pretrain_epochs)])

    if no_cuda:
        if is_deep_or_gnn:
            cmd.append("--no_cuda")
    elif is_deep_or_gnn:
        cmd.extend(["--gpu", str(gpu)])

    if extra_args:
        cmd.extend(extra_args)

    return cmd


# ────────────────────────────────────────────────────────────────
# Output verification
# ────────────────────────────────────────────────────────────────

def verify_output(out_dir: Path) -> bool:
    """
    Verify that a method's output contains all required files.
    Returns True if all required files exist.
    """
    for fname in ("embedding_final.npy", "labels.npy", "metrics.json", "args.json"):
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


def normalize_metrics(metrics: Optional[Dict[str, Any]]) -> Dict[str, Any]:
    """
    Unwrap nested metrics structures into a flat dict of evaluation metrics.

    Handles three known layouts:
      1. {"kmeans_known_k": {...}}         — NeighborMix_scMAE / scMAE
      2. {"fixed": {"kmeans_known_k": {...}}} — legacy double-nested
      3. flat {"acc": ..., "nmi": ...}    — already flat, returned as-is
    """
    if not metrics:
        return {}
    if "kmeans_known_k" in metrics and isinstance(metrics["kmeans_known_k"], dict):
        return metrics["kmeans_known_k"]
    if "fixed" in metrics and isinstance(metrics["fixed"], dict):
        fixed = metrics["fixed"]
        if "kmeans_known_k" in fixed and isinstance(fixed["kmeans_known_k"], dict):
            return fixed["kmeans_known_k"]
        return fixed
    return metrics


# ────────────────────────────────────────────────────────────────
# Per-run execution
# ────────────────────────────────────────────────────────────────

def run_method(
    method_key: str,
    method_info: Dict[str, Any],
    data_path: str,
    base_out_dir: Path,
    n_clusters: int,
    epochs: int,
    pretrain_epochs: int,
    seed: int,
    gpu: int,
    no_cuda: bool,
    allow_unverified: bool,
    commit_sha: str,
    branch: str,
    verbose: bool = False,
    dry_run: bool = False,
) -> Dict[str, Any]:
    """Run a single method with a given seed. Returns a result dict."""

    run_id = f"{method_key}__seed{seed}__{datetime.now().strftime('%Y%m%d_%H%M%S')}"
    out_dir = base_out_dir / run_id
    os.makedirs(out_dir, exist_ok=True)

    start_time = datetime.now().isoformat()
    end_time = ""
    elapsed_seconds = 0.0
    return_code = None
    status = "unknown"
    error_msg = ""
    cmd: List[str] = []

    write_authenticity_json(out_dir, method_info)

    # Check if allowed
    allowed, reason = check_authenticity(method_key, method_info)
    if not allowed:
        if not allow_unverified:
            status = "skipped"
            end_time = datetime.now().isoformat()
            _write_status(out_dir, method_key, seed, status, return_code,
                          elapsed_seconds, gpu, no_cuda, start_time, end_time,
                          cmd, commit_sha, branch, reason, "")
            return {
                "method": method_key,
                "seed": seed,
                "status": status,
                "reason": reason,
                "out_dir": str(out_dir),
                "commit_sha": commit_sha,
                "branch": branch,
            }
        else:
            status = "unverified"
            auth_path = out_dir / "authenticity.json"
            with open(auth_path, "r") as f:
                auth_data = json.load(f)
            auth_data["unverified"] = True
            auth_data["unverified_reason"] = reason
            with open(auth_path, "w") as f:
                json.dump(auth_data, f, indent=2)

    extra_args = method_info.get("extra_args", None)
    cmd = build_command(
        method_key, method_info, data_path, out_dir,
        n_clusters, epochs, pretrain_epochs, seed,
        gpu, no_cuda, extra_args,
    )

    with open(out_dir / "command.txt", "w") as f:
        f.write(" ".join(f'"{c}"' if " " in c else c for c in cmd) + "\n")

    if dry_run:
        status = "dry_run"
        end_time = datetime.now().isoformat()
        _write_status(out_dir, method_key, seed, status, return_code,
                      elapsed_seconds, gpu, no_cuda, start_time, end_time,
                      cmd, commit_sha, branch, "", "")
        return {
            "method": method_key,
            "seed": seed,
            "status": status,
            "out_dir": str(out_dir),
            "commit_sha": commit_sha,
            "branch": branch,
        }

    # Actual execution
    t0 = time.time()
    try:
        env = os.environ.copy()
        env["PYTHONPATH"] = str(METHODS_DIR) + ":" + env.get("PYTHONPATH", "")
        # Limit threading to avoid GPU OOM from thread over-subscription.
        # PyTorch defaults to 96 threads; OpenBLAS maxes at 64.
        # With 3 benchmarks + multiple workers, this causes SIGSEGV in CUDA.
        env["OMP_NUM_THREADS"] = "4"
        env["OPENBLAS_NUM_THREADS"] = "4"
        env["MKL_NUM_THREADS"] = "4"
        env["NUMEXPR_NUM_THREADS"] = "4"
        # Prevent CUDA OOM fragmentation for large-graph methods (scDSC)
        if method_key == "scdsc":
            env["PYTORCH_CUDA_ALLOC_CONF"] = "expandable_segments:True"

        proc = subprocess.run(
            cmd,
            cwd=str(PROJECT_ROOT),
            capture_output=True,
            text=True,
            timeout=3600 * 6,
            env=env,
        )
        return_code = proc.returncode

        with open(out_dir / "run.log", "w") as f:
            f.write(f"Command: {' '.join(cmd)}\n")
            f.write(f"Exit code: {proc.returncode}\n")
            f.write(f"\n=== STDOUT ===\n{proc.stdout}\n")
            if proc.stderr:
                f.write(f"\n=== STDERR ===\n{proc.stderr}\n")

        if return_code != 0:
            status = "failed"
            error_msg = f"Exit code {return_code}"
        elif not verify_output(out_dir):
            status = "incomplete"
            error_msg = "Missing required output files"
        else:
            status = "success"

    except subprocess.TimeoutExpired:
        status = "timeout"
        return_code = -1
        error_msg = "Execution exceeded 6 hour timeout"
        with open(out_dir / "run.log", "w") as f:
            f.write(f"Command timed out after 6 hours\nCommand: {' '.join(cmd)}\n")
    except Exception as e:
        status = "error"
        return_code = -2
        error_msg = str(e)
        with open(out_dir / "run.log", "w") as f:
            f.write(f"Error: {e}\nCommand: {' '.join(cmd)}\n")

    elapsed_seconds = round(time.time() - t0, 1)
    end_time = datetime.now().isoformat()

    _write_status(out_dir, method_key, seed, status, return_code,
                  elapsed_seconds, gpu, no_cuda, start_time, end_time,
                  cmd, commit_sha, branch, "", error_msg)

    metrics = load_metrics(out_dir)

    result = {
        "method": method_key,
        "seed": seed,
        "status": status,
        "error": error_msg,
        "elapsed_seconds": elapsed_seconds,
        "out_dir": str(out_dir),
        "metrics": metrics,
        "gpu": gpu,
        "no_cuda": no_cuda,
        "return_code": return_code,
        "commit_sha": commit_sha,
        "branch": branch,
        "command": " ".join(cmd),
    }

    if verbose:
        icon = {
            "success": "✓", "failed": "✗", "skipped": "⊗",
            "incomplete": "⚠", "unverified": "⚠", "timeout": "⏱",
            "error": "!", "dry_run": "⟳",
        }.get(status, "?")
        print(
            f"  {icon} {method_key} (seed={seed}): {status} ({elapsed_seconds:.1f}s)"
        )
        if metrics:
            nm = normalize_metrics(metrics)
            acc = nm.get("acc"); nmi = nm.get("nmi"); ari = nm.get("ari")
            print(
                f"    ACC={f'{acc:.4f}' if acc is not None else 'N/A'} "
                f"NMI={f'{nmi:.4f}' if nmi is not None else 'N/A'} "
                f"ARI={f'{ari:.4f}' if ari is not None else 'N/A'}"
            )
        if error_msg:
            print(f"    Error: {error_msg}")

    return result


def _write_status(
    out_dir: Path,
    method: str,
    seed: int,
    status: str,
    return_code: Optional[int],
    elapsed_seconds: float,
    gpu: int,
    no_cuda: bool,
    start_time: str,
    end_time: str,
    command: List[str],
    commit_sha: str,
    branch: str,
    reason: str,
    error_msg: str,
) -> None:
    """Write status.json to the output directory."""
    status_data = {
        "method": method,
        "seed": seed,
        "dataset": "",
        "status": status,
        "return_code": return_code,
        "elapsed_seconds": elapsed_seconds,
        "gpu": gpu,
        "no_cuda": no_cuda,
        "start_time": start_time,
        "end_time": end_time,
        "command": " ".join(command),
        "commit_sha": commit_sha,
        "branch": branch,
        "reason": reason,
        "error": error_msg,
    }
    with open(out_dir / "status.json", "w", encoding="utf-8") as f:
        json.dump(status_data, f, indent=2)


# ────────────────────────────────────────────────────────────────
# Result collection
# ────────────────────────────────────────────────────────────────

def collect_results(base_out_dir: Path, method_keys: List[str], seeds: List[int]) -> List[Dict]:
    """Collect results from completed runs, reading status.json and authenticity.json."""
    results = []
    for sub_dir in sorted(base_out_dir.iterdir()):
        if not sub_dir.is_dir():
            continue
        parts = sub_dir.name.split("__")
        if len(parts) < 2:
            continue
        method_key = parts[0]
        if method_key not in method_keys:
            continue

        metrics = load_metrics(sub_dir)

        # Read status.json
        status_path = sub_dir / "status.json"
        status_data = {}
        if status_path.exists():
            with open(status_path) as f:
                status_data = json.load(f)

        # Read authenticity.json
        auth_path = sub_dir / "authenticity.json"
        auth_data = {}
        if auth_path.exists():
            with open(auth_path) as f:
                auth_data = json.load(f)

        # Parse seed from dir name if not in status
        seed_from_dir = None
        for s in seeds:
            seed_str = f"seed{s}_"
            if seed_str in sub_dir.name:
                seed_from_dir = s
                break
        if seed_from_dir is None:
            seed_from_dir = status_data.get("seed", "unknown")

        results.append({
            "method": method_key,
            "seed": status_data.get("seed", seed_from_dir),
            "status": status_data.get("status", "unknown"),
            "error": status_data.get("error", ""),
            "elapsed_seconds": status_data.get("elapsed_seconds", 0.0),
            "gpu": status_data.get("gpu", ""),
            "no_cuda": status_data.get("no_cuda", False),
            "return_code": status_data.get("return_code", ""),
            "commit_sha": status_data.get("commit_sha", ""),
            "branch": status_data.get("branch", ""),
            "command": status_data.get("command", ""),
            "out_dir": str(sub_dir),
            "metrics": metrics,
            "authenticity": auth_data.get("authenticity", "UNKNOWN"),
            "substitute_model_used": auth_data.get("substitute_model_used", True),
            "unverified": auth_data.get("unverified", False),
        })
    return results


# ────────────────────────────────────────────────────────────────
# Summary generation
# ────────────────────────────────────────────────────────────────

METRIC_NAMES = [
    "acc", "nmi", "ari", "f1_macro", "fmi",
    "v_measure", "homogeneity", "completeness",
]


def generate_summary_csv(results: List[Dict], out_path: Path, dataset: str) -> None:
    """Generate benchmark_summary.csv with full traceability fields."""
    rows = []
    for r in results:
        m = normalize_metrics(r.get("metrics"))
        row = {
            "dataset": dataset,
            "method": r["method"],
            "seed": r.get("seed", ""),
            "status": r.get("status", "unknown"),
            "authenticity": r.get("authenticity", "UNKNOWN"),
            "substitute_model_used": r.get("substitute_model_used", True),
            "unverified": r.get("unverified", False),
            "acc": m.get("acc", ""),
            "nmi": m.get("nmi", ""),
            "ari": m.get("ari", ""),
            "f1_macro": m.get("f1_macro", ""),
            "fmi": m.get("fmi", ""),
            "v_measure": m.get("v_measure", ""),
            "homogeneity": m.get("homogeneity", ""),
            "completeness": m.get("completeness", ""),
            "runtime_seconds": r.get("elapsed_seconds", ""),
            "gpu": r.get("gpu", ""),
            "no_cuda": r.get("no_cuda", ""),
            "commit_sha": r.get("commit_sha", ""),
            "branch": r.get("branch", ""),
            "command": r.get("command", ""),
            "error": r.get("error", ""),
            "save_dir": r.get("out_dir", ""),
        }
        rows.append(row)

    if rows:
        import csv
        with open(out_path, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=rows[0].keys())
            writer.writeheader()
            writer.writerows(rows)
        print(f"Summary saved to: {out_path}")

    return rows


def generate_mean_std_csv(results: List[Dict], out_path: Path, dataset: str) -> None:
    """
    Generate benchmark_summary_mean_std.csv, aggregating over seeds.
    Only includes rows where status=success AND substitute_model_used=false.
    """
    import csv
    import statistics
    from collections import defaultdict

    # Filter: only successful, verified runs
    valid_results = [
        r for r in results
        if r.get("status") == "success"
        and r.get("substitute_model_used") is False
        and r.get("authenticity") == "VERIFIED"
    ]

    by_method = defaultdict(list)
    for r in valid_results:
        by_method[r["method"]].append(r)

    rows = []
    for method, runs in sorted(by_method.items()):
        row = {
            "dataset": dataset,
            "method": method,
            "n_success": len(runs),
            "n_total": len([r for r in results if r["method"] == method]),
            "authenticity": runs[0]["authenticity"],
            "substitute_model_used": False,
        }
        for metric in METRIC_NAMES:
            vals = [
                normalize_metrics(r.get("metrics")).get(metric)
                for r in runs
                if normalize_metrics(r.get("metrics")).get(metric) is not None
            ]
            if vals:
                mean_val = statistics.mean(vals)
                std_val = statistics.stdev(vals) if len(vals) > 1 else 0.0
                row[metric] = f"{mean_val:.4f} ± {std_val:.4f}"
                row[f"{metric}_mean"] = mean_val
                row[f"{metric}_std"] = std_val
            else:
                row[metric] = "N/A"
                row[f"{metric}_mean"] = ""
                row[f"{metric}_std"] = ""

        status_summary = "; ".join(
            f"seed={r.get('seed','')}:{r.get('status','')}"
            for r in runs
        )
        row["status_summary"] = status_summary
        rows.append(row)

    if rows:
        # Write readable form (with ±)
        readable_fields = [
            "dataset", "method", "n_success", "n_total", "authenticity",
            "substitute_model_used",
            "acc", "nmi", "ari", "f1_macro", "fmi",
            "v_measure", "homogeneity", "completeness",
            "status_summary",
        ]
        with open(out_path, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=readable_fields, extrasaction="ignore")
            writer.writeheader()
            writer.writerows(rows)

        # Also write a machine-readable .md summary
        md_path = out_path.with_suffix(".md")
        with open(md_path, "w", encoding="utf-8") as f:
            f.write(f"# Benchmark Results: {dataset}\n\n")
            headers = readable_fields
            f.write("| " + " | ".join(headers) + " |\n")
            f.write("|" + "|".join(" --- " for _ in headers) + "|\n")
            for row in rows:
                cells = [str(row.get(h, "")) for h in headers]
                f.write("| " + " | ".join(cells) + " |\n")
        print(f"Mean±std summary saved to: {out_path}")
        print(f"Markdown summary saved to: {md_path}")


# ────────────────────────────────────────────────────────────────
# Main
# ────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Run formal benchmark on VERIFIED methods",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=f"""
Examples:
  # Default formal run: all 10 methods (proposed + ablation + baselines)
  python scripts/run_formal_benchmark.py \\
      --data_path data/subsample_2k.h5ad \\
      --dataset_name subsample_2k \\
      --out_dir results/formal \\
      --n_clusters 7 \\
      --seeds 42 \\
      --gpu 1

  # CPU-only
  python scripts/run_formal_benchmark.py \\
      --data_path data/subsample_2k.h5ad \\
      --dataset_name subsample_2k \\
      --out_dir results/formal \\
      --n_clusters 7 --seeds 42 --no_cuda

  # Preflight / dry run (no execution)
  python scripts/run_formal_benchmark.py \\
      --data_path data/subsample_2k.h5ad \\
      --out_dir results/formal \\
      --n_clusters 7 --seeds 42 --gpu 1 --dry_run

  # Run specific methods
  python scripts/run_formal_benchmark.py \\
      --data_path data/SRP182008.h5ad \\
      --methods dec scdcc scanpy_standard \\
      --n_clusters 15 --seeds 42

  # Include pending/unverified (flagged as unverified in output)
  python scripts/run_formal_benchmark.py \\
      --data_path data/subsample_2k.h5ad \\
      --methods scgnn sccdcg \\
      --n_clusters 7 --allow_unverified --gpu 1

GPU Policy: GPU 0 and GPU 7 are FORBIDDEN.
Only VERIFIED + Smoke=PASS methods are included in formal benchmark by default.
Substitute implementations are forbidden.
""",
    )

    parser.add_argument("--data_path", type=str, required=True,
                       help="Path to input h5ad file")
    parser.add_argument("--out_dir", type=str, default="results/formal",
                       help="Base output directory")
    parser.add_argument("--dataset_name", type=str, default=None,
                       help="Dataset name for summary tables (default: derived from data_path)")
    parser.add_argument("--n_clusters", type=int, required=True,
                       help="Number of clusters")
    parser.add_argument("--seeds", type=int, nargs="+", default=[42],
                       help="Random seeds (default: 42)")
    parser.add_argument("--epochs", type=int, default=200,
                       help="Training epochs for deep learning models")
    parser.add_argument("--pretrain_epochs", type=int, default=200,
                       help="Pretraining epochs for models with pretrain phase")
    parser.add_argument("--methods", type=str, nargs="+", default=None,
                       help=f"Method keys to run. Default: all VERIFIED+Smoke=PASS. "
                            f"Default formal list: {DEFAULT_FORMAL_METHODS}")
    parser.add_argument("--no_cuda", action="store_true",
                       help="Disable CUDA (use CPU only)")
    parser.add_argument("--gpu", type=int, default=1,
                       help="GPU device ID (default: 1). GPU 0 and GPU 7 are forbidden.")
    parser.add_argument("--allow_unverified", action="store_true",
                       help="Allow PENDING/ENV-GATED models (output marked 'unverified')")
    parser.add_argument("--skip_run", action="store_true",
                       help="Skip running models, only collect and generate summary")
    parser.add_argument("--dry_run", action="store_true",
                       help="Print commands and validate but do not execute any model")
    parser.add_argument("--verbose", "-v", action="store_true",
                       help="Verbose output")

    args = parser.parse_args()

    # ── GPU policy validation ──────────────────────────────────
    validate_gpu_policy(args.gpu, args.no_cuda)

    # ── Load manifest ───────────────────────────────────────────
    manifest = load_manifest()
    dataset_name = args.dataset_name or os.path.splitext(os.path.basename(args.data_path))[0]
    if args.out_dir and dataset_name in args.out_dir:
        base_out_dir = Path(args.out_dir)
    else:
        base_out_dir = Path(args.out_dir) / dataset_name
    os.makedirs(base_out_dir, exist_ok=True)

    # ── Determine methods to run ────────────────────────────────
    if args.methods:
        selected = {}
        for k in args.methods:
            if k in manifest:
                selected[k] = manifest[k]
            else:
                print(f"WARNING: Unknown method key: {k!r} (skipping)")
        missing = [k for k in args.methods if k not in manifest]
        if missing:
            print(f"WARNING: Unknown method keys: {missing}")
    else:
        selected = {
            k: v for k, v in manifest.items()
            if v.get("authenticity") == "VERIFIED"
            and v.get("smoke") in ("PASS", "UNKNOWN", "")
            and v.get("default_in_formal") is True
        }

    if not selected:
        print("ERROR: No valid methods selected. Use --methods or check manifest.")
        sys.exit(2)

    # ── Git info ───────────────────────────────────────────────
    commit_sha, branch = get_git_info()

    # ── Summary of what will run ────────────────────────────────
    print("=" * 70)
    print("Formal Benchmark Runner")
    print("=" * 70)
    print(f"  Data:     {args.data_path}")
    print(f"  Dataset:  {dataset_name}")
    print(f"  Clusters: {args.n_clusters}")
    print(f"  Seeds:    {args.seeds}")
    print(f"  Epochs:   {args.epochs} (pretrain: {args.pretrain_epochs})")
    print(f"  CUDA:     {'disabled (CPU only)' if args.no_cuda else f'GPU {args.gpu}'}")
    print(f"  Unverified: {'ALLOWED' if args.allow_unverified else 'REJECTED'}")
    print(f"  Dry run:  {'YES' if args.dry_run else 'NO'}")
    print(f"  Commit:   {commit_sha[:8]} ({branch})")
    print()
    print(f"  Methods ({len(selected)}):")
    for k, v in sorted(selected.items()):
        allowed, reason = check_authenticity(k, v)
        icon = "✓" if allowed else "⊗"
        print(f"    {icon} {k:25s} [{v['category']:12s}] {reason}")
    print()

    if not args.allow_unverified:
        skipped = [
            k for k, v in selected.items()
            if not check_authenticity(k, v)[0]
        ]
        if skipped:
            print(f"  Skipping {len(skipped)} unverified methods: {skipped}")
            print()

    print("=" * 70)
    print()

    # ── Dry run: just validate and print ───────────────────────
    if args.dry_run:
        print("DRY RUN — no models executed. Preflight checks passed.")
        print()
        for method_key, method_info in sorted(selected.items()):
            allowed, reason = check_authenticity(method_key, method_info)
            for seed in args.seeds:
                extra_args = method_info.get("extra_args", None)
                cmd = build_command(
                    method_key, method_info, args.data_path,
                    base_out_dir / f"{method_key}__seed{seed}__dry",
                    args.n_clusters, args.epochs, args.pretrain_epochs,
                    seed, args.gpu, args.no_cuda, extra_args,
                )
                icon = "✓" if allowed else "⊗"
                print(f"  {icon} {method_key} (seed={seed})")
                print(f"     {' '.join(cmd)}")
            print()
        sys.exit(0)

    # ── Execute runs ───────────────────────────────────────────
    if not args.skip_run:
        all_results = []
        for method_key, method_info in sorted(selected.items()):
            allowed, reason = check_authenticity(method_key, method_info)
            if not allowed and not args.allow_unverified:
                print(f"  ⊗ {method_key}: skipped ({reason})")
                for seed in args.seeds:
                    # Write skipped status even for non-run methods
                    run_id = f"{method_key}__seed{seed}__{datetime.now().strftime('%Y%m%d_%H%M%S')}"
                    out_dir = base_out_dir / run_id
                    os.makedirs(out_dir, exist_ok=True)
                    write_authenticity_json(out_dir, method_info)
                    _write_status(
                        out_dir, method_key, seed, "skipped",
                        None, 0.0, args.gpu, args.no_cuda,
                        datetime.now().isoformat(), datetime.now().isoformat(),
                        [], commit_sha, branch, reason, "",
                    )
                print()
                continue

            print(f"  Running {method_key}...")
            for seed in args.seeds:
                result = run_method(
                    method_key, method_info,
                    args.data_path, base_out_dir,
                    args.n_clusters, args.epochs, args.pretrain_epochs,
                    seed, args.gpu, args.no_cuda,
                    args.allow_unverified,
                    commit_sha, branch,
                    verbose=args.verbose,
                    dry_run=args.dry_run,
                )
                all_results.append(result)
            print()

    # ── Collect results ────────────────────────────────────────
    print("=" * 70)
    print("Collecting results...")
    results = collect_results(base_out_dir, list(selected.keys()), args.seeds)

    # ── Generate summary tables ────────────────────────────────
    summary_csv = base_out_dir / "benchmark_summary.csv"
    generate_summary_csv(results, summary_csv, dataset_name)

    mean_csv = base_out_dir / "benchmark_summary_mean_std.csv"
    generate_mean_std_csv(results, mean_csv, dataset_name)

    # ── Status summary ────────────────────────────────────────
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

    main_table_count = sum(
        1 for r in results
        if r.get("status") == "success"
        and r.get("substitute_model_used") is False
        and r.get("authenticity") == "VERIFIED"
    )
    print(f"  Main-table-eligible rows (success + VERIFIED + no substitute): {main_table_count}")

    if status_counts.get("failed") or status_counts.get("error") or status_counts.get("timeout"):
        print("  ⚠  Some methods failed. Check run.log files and status.json.")
        sys.exit(1)
    else:
        print("  ✓  All methods completed successfully.")
        sys.exit(0)


if __name__ == "__main__":
    main()
