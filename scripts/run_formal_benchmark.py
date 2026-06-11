#!/usr/bin/env python3
"""
Formal Benchmark Runner
=======================
Per BDD Scenarios 1-16: Runs the curated DEFAULT_FORMAL_METHODS list
by default. Method metadata is read from methods/method_manifest.yaml.

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
import re
import shutil
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Optional, Any

PROJECT_ROOT = Path(__file__).parent.parent.resolve()
METHODS_DIR = PROJECT_ROOT / "methods"
MANIFEST_PATH = METHODS_DIR / "method_manifest.yaml"
SCRIPTS_DIR = PROJECT_ROOT / "scripts"
PROCESSED_DATA_DIR = PROJECT_ROOT / "data" / "processed"
FORBIDDEN_GPUS = {"0", "7"}
DEFAULT_FORMAL_METHODS = [
    "neighbormix_scmae",
    "nm_scmae_nomix",
    "scmae",
    "dec",
    "scdcc",
    "scdsc",
    "leiden",
    "louvain",
    "sc3",
]


# ────────────────────────────────────────────────────────────────
# Runtime Registry (Phase 1 / Scenario 1.2)
# ────────────────────────────────────────────────────────────────

def load_runtime_registry(registry_path: Optional[Path]) -> Dict[str, dict]:
    """Load runtime registry YAML. Returns empty dict if not found or not specified."""
    if not registry_path or not registry_path.exists():
        return {}
    try:
        with open(registry_path, "r", encoding="utf-8") as f:
            data = yaml.safe_load(f)
        return data.get("runtimes", {})
    except Exception as e:
        print(f"  WARNING: Could not load runtime registry {registry_path}: {e}", file=sys.stderr)
        return {}


def get_runtime_env(method_info: Dict[str, Any]) -> str:
    """
    Extract the logical runtime name from method_info.

    Supports both flat (runtime_env) and nested (runtime.env_name) formats.
    """
    # Flat format: runtime_env = "plantnet-tf1"
    if method_info.get("runtime_env"):
        return str(method_info["runtime_env"])
    # Nested format: runtime = {env_name: "plantnet-tf1", ...}
    runtime = method_info.get("runtime") or {}
    if isinstance(runtime, dict):
        return runtime.get("env_name", "")
    return ""


def resolve_python_executable(
    method_key: str,
    method_info: Dict[str, Any],
    runtime_registry: Dict[str, dict],
) -> str:
    """
    Per BDD Scenario 2.1: resolve the correct Python executable for a method.

    Priority:
      1. method_info["runtime"]["python"] (absolute path, if present and exists)
      2. runtime_registry[runtime_env]["python"] (from registry)
      3. sys.executable (fallback: current environment)

    The registry maps logical names (e.g. 'plantnet-tf1') to concrete Python paths.
    Logical name is read from runtime.env_name or runtime_env field.
    """
    runtime_env = get_runtime_env(method_info)

    # 1. Direct python path in manifest runtime field
    runtime = method_info.get("runtime") or {}
    if isinstance(runtime, dict):
        python_path = runtime.get("python", "")
        if python_path and Path(python_path).exists():
            return python_path

    # 2. Lookup via runtime_registry[runtime_env]
    if runtime_env and runtime_env in runtime_registry:
        entry = runtime_registry[runtime_env]
        python_path = entry.get("python", "")
        if python_path and Path(python_path).exists():
            return python_path
        print(f"  WARNING: runtime '{runtime_env}' python not found: {python_path}, using default", file=sys.stderr)

    # 3. Fallback
    return sys.executable


def write_environment_json(
    out_dir: Path,
    python_executable: str,
    runtime_env: str = "",
    framework: str = "",
) -> None:
    """
    Per BDD Scenario 2.3: write environment.json with version info.
    """
    env_data = {
        "runtime_backend": "conda",
        "runtime_env": runtime_env or "default",
        "python_executable": python_executable,
        "python_version": "",
        "cuda_visible_devices": os.environ.get("CUDA_VISIBLE_DEVICES", ""),
        "torch_version": "",
        "tensorflow_version": "",
        "scanpy_version": "",
        "anndata_version": "",
        "version_errors": {},
    }

    # Collect version info from the actual Python executable
    exe = python_executable or sys.executable
    version_cmds = [
        ("python_version", [exe, "-c", "import sys; print(sys.version.split()[0])"]),
    ]
    framework_lower = str(framework).lower()
    if any(token in framework_lower for token in ("torch", "pytorch")):
        version_cmds.append(("torch_version", [exe, "-c", "import torch; print(torch.__version__)"]))
    if any(token in framework_lower for token in ("tensorflow", "keras")):
        version_cmds.append(("tensorflow_version", [exe, "-c", "import tensorflow; print(tensorflow.__version__)"]))
    if "scanpy" in framework_lower:
        version_cmds.extend([
            ("scanpy_version", [exe, "-c", "import scanpy; print(scanpy.__version__)"]),
            ("anndata_version", [exe, "-c", "import anndata; print(anndata.__version__)"]),
        ])
    for key, cmd in version_cmds:
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=3)
            if result.returncode == 0:
                env_data[key] = result.stdout.strip()
            else:
                env_data["version_errors"][key] = result.stderr.strip()[-500:]
        except subprocess.TimeoutExpired:
            env_data["version_errors"][key] = "version probe timed out"
        except Exception as exc:
            env_data["version_errors"][key] = str(exc)

    with open(out_dir / "environment.json", "w", encoding="utf-8") as f:
        json.dump(env_data, f, indent=2)


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


def source_fingerprint(path: str) -> Dict[str, Any]:
    """Return stable source-file identity fields for conversion cache checks."""
    resolved = os.path.abspath(path)
    stat = os.stat(resolved)
    return {
        "input_path": resolved,
        "input_size": int(stat.st_size),
        "input_mtime_ns": int(stat.st_mtime_ns),
    }


def conversion_cache_matches(meta: Dict[str, Any], data_path: str, dataset_name: str) -> bool:
    """Check whether an existing converted h5ad was created from this exact source file."""
    expected = {
        "dataset_name": dataset_name,
        **source_fingerprint(data_path),
        "conversion_label_key": "auto",
        "conversion_matrix_key": "auto",
        "conversion_n_clusters": "auto",
    }
    return all(meta.get(key) == value for key, value in expected.items())


def auto_convert_h5(data_path: str, dataset_name: str) -> tuple[str, int]:
    """
    If data_path ends with .h5, convert to .h5ad via prepare_dataset.py and return
    the path to the converted file and inferred n_clusters.

    Returns (converted_path, n_clusters).

    If data_path ends with .h5ad, just loads and infers n_clusters (requires
    --n_clusters auto or a real value), returning (data_path, n_clusters).

    Exits on failure.
    """
    import scanpy as sc

    ext = os.path.splitext(data_path)[1].lower()

    if ext == ".h5":
        PROCESSED_DATA_DIR.mkdir(parents=True, exist_ok=True)
        output_path = PROCESSED_DATA_DIR / f"{dataset_name}.h5ad"
        meta_path = PROCESSED_DATA_DIR / f"{dataset_name}.meta.json"

        if output_path.exists() and meta_path.exists():
            with open(meta_path) as f:
                meta = json.load(f)
            if conversion_cache_matches(meta, data_path, dataset_name):
                print(f"  [auto-convert] Reusing cached: {output_path}")
                return str(output_path), meta["n_clusters"]
            print(f"  [auto-convert] Cache is stale for {dataset_name}; regenerating.")

        # Build prepare_dataset.py command
        prepare_script = SCRIPTS_DIR / "prepare_dataset.py"
        cmd = [
            sys.executable,
            str(prepare_script),
            "--input_path", data_path,
            "--dataset_name", dataset_name,
            "--output_dir", str(PROCESSED_DATA_DIR),
            "--force",  # always regenerate fresh conversion
        ]
        print(f"\n  [auto-convert] .h5 detected. Converting via prepare_dataset.py ...")
        print(f"    Input:  {data_path}")
        print(f"    Output: {output_path}")
        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=600,
                env=os.environ.copy(),
            )
            if result.returncode != 0:
                print(f"  [auto-convert] FAILED:\n{result.stderr}\n{result.stdout}")
                raise SystemExit(1)
            # Read back n_clusters from meta.json
            with open(meta_path) as f:
                meta = json.load(f)
            n_clusters = meta["n_clusters"]
            print(f"  [auto-convert] Done. shape=({meta['n_cells']}, {meta['n_genes']}), "
                  f"n_clusters={n_clusters}")
            return str(output_path), n_clusters
        except subprocess.TimeoutExpired:
            print("  [auto-convert] TIMEOUT (10 min)")
            raise SystemExit(1)
        except Exception as e:
            print(f"  [auto-convert] ERROR: {e}")
            raise SystemExit(1)

    elif ext == ".h5ad":
        # n_clusters is handled separately; just validate the file is readable
        try:
            ad = sc.read_h5ad(data_path)
            return data_path, None  # caller must provide n_clusters
        except Exception as e:
            print(f"  [data] ERROR reading {data_path}: {e}")
            raise SystemExit(1)

    else:
        print(f"  [data] Unsupported extension: {ext} (expected .h5 or .h5ad)")
        raise SystemExit(1)


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
    elif auth == "PENDING_AUDITED":
        return False, "PENDING_AUDITED (audited, smoke not yet run)"
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
    python_bin: str,
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
    """Build the command-line invocation for a method.

    Args:
        python_bin: Path to the Python executable to use.
                    Use resolved_python_executable(), not sys.executable.
    """

    run_py = PROJECT_ROOT / method_info["path"]

    cmd = [
        python_bin,
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

def required_artifacts(method_info: Dict[str, Any]) -> List[str]:
    """Return required output artifacts for a method."""
    required = method_info.get("required_artifacts")
    if required is None:
        return ["metrics.json"]
    if isinstance(required, str):
        return [required]
    return [str(item) for item in required]


def verify_output(out_dir: Path, method_info: Dict[str, Any]) -> tuple[bool, str]:
    """
    Verify that a method's output contains its required files.
    Returns (ok, reason).
    """
    missing = [fname for fname in required_artifacts(method_info) if not (out_dir / fname).exists()]
    if missing:
        return False, "Missing required output files: " + ", ".join(missing)
    return True, ""


def load_metrics(out_dir: Path) -> Optional[Dict[str, float]]:
    """Load metrics.json if it exists."""
    metrics_path = out_dir / "metrics.json"
    if metrics_path.exists():
        with open(metrics_path, "r") as f:
            return json.load(f)
    return None


def parse_seed_from_run_id(run_id: str) -> Optional[int]:
    """Parse seed from run IDs like method__seed42__20260611_120000."""
    match = re.search(r"(?:^|__)seed(-?\d+)(?:__|_|$)", run_id)
    if not match:
        return None
    return int(match.group(1))


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
    dataset_name: str = "",
    python_bin: str = "",
    runtime_env: str = "",
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
    write_environment_json(
        out_dir,
        python_bin or sys.executable,
        runtime_env,
        framework=method_info.get("framework", ""),
    )

    # Check if allowed
    allowed, reason = check_authenticity(method_key, method_info)
    if not allowed:
        if not allow_unverified:
            status = "skipped"
            end_time = datetime.now().isoformat()
            _write_status(out_dir, method_key, seed, status, return_code,
                          elapsed_seconds, gpu, no_cuda, start_time, end_time,
                          cmd, commit_sha, branch, reason, "",
                          dataset=dataset_name, n_clusters=n_clusters,
                          python_executable=python_bin or sys.executable,
                          runtime_env=runtime_env)
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
        python_bin or sys.executable,
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
                      cmd, commit_sha, branch, "", "",
                      dataset=dataset_name, n_clusters=n_clusters,
                      python_executable=python_bin or sys.executable,
                      runtime_env=runtime_env)
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
        env.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")
        env.setdefault("NUMBA_CACHE_DIR", "/tmp/numba-cache")
        os.makedirs(env["MPLCONFIGDIR"], exist_ok=True)
        os.makedirs(env["NUMBA_CACHE_DIR"], exist_ok=True)
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
        else:
            ok, reason = verify_output(out_dir, method_info)
            if ok:
                status = "success"
            else:
                status = "incomplete"
                error_msg = reason

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
                  cmd, commit_sha, branch, "", error_msg,
                  dataset=dataset_name, n_clusters=n_clusters,
                  python_executable=python_bin or sys.executable,
                  runtime_env=runtime_env)

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
        "python_executable": python_bin or sys.executable,
        "runtime_env": runtime_env,
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
    dataset: str = "",
    n_clusters: int = None,
    python_executable: str = "",
    runtime_env: str = "",
) -> None:
    """Write status.json to the output directory."""
    status_data = {
        "method": method,
        "seed": seed,
        "dataset": dataset,
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
        "runtime_env": runtime_env,
        "python_executable": python_executable,
    }
    if n_clusters is not None:
        status_data["n_clusters"] = n_clusters
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

        # Parse seed from dir name if not in status.
        seed_from_dir = parse_seed_from_run_id(sub_dir.name)
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
            # Runtime isolation fields (written by _write_status / run_method)
            "runtime_env": status_data.get("runtime_env", ""),
            "python_executable": status_data.get("python_executable", ""),
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
            "runtime_env": r.get("runtime_env", ""),
            "python_executable": r.get("python_executable", ""),
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
  # Default formal run: curated methods (proposed + ablation + baselines)
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
      --methods dec scdcc leiden \\
      --n_clusters 15 --seeds 42

  # Include pending/unverified (flagged as unverified in output)
  python scripts/run_formal_benchmark.py \\
      --data_path data/subsample_2k.h5ad \\
      --methods scgnn sccdcg \\
      --n_clusters 7 --allow_unverified --gpu 1

GPU Policy: GPU 0 and GPU 7 are FORBIDDEN.
Only the curated DEFAULT_FORMAL_METHODS list is included by default.
Substitute implementations are forbidden.
""",
    )

    parser.add_argument("--data_path", type=str, required=True,
                       help="Path to input h5ad file")
    parser.add_argument("--out_dir", type=str, default="results/formal",
                       help="Base output directory")
    parser.add_argument("--dataset_name", type=str, default=None,
                       help="Dataset name for summary tables (default: derived from data_path)")
    parser.add_argument("--n_clusters", type=str, required=True,
                       help="Number of clusters. Use 'auto' to infer from labels in .h5 or .h5ad.")
    parser.add_argument("--seeds", type=int, nargs="+", default=[42],
                       help="Random seeds (default: 42)")
    parser.add_argument("--epochs", type=int, default=200,
                       help="Training epochs for deep learning models")
    parser.add_argument("--pretrain_epochs", type=int, default=200,
                       help="Pretraining epochs for models with pretrain phase")
    parser.add_argument("--methods", type=str, nargs="+", default=None,
                       help=f"Method keys to run. Default formal list: {DEFAULT_FORMAL_METHODS}")
    parser.add_argument("--runtime_registry", type=str, default=None,
                       help="Path to runtime_registry.yaml for per-method Python resolution. "
                            "Default: envs/runtime_registry.yaml if exists.")
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

    # ── Auto-convert .h5 → .h5ad ─────────────────────────────
    dataset_name = args.dataset_name or os.path.splitext(os.path.basename(args.data_path))[0]

    data_path = os.path.abspath(args.data_path)
    n_clusters_raw = args.n_clusters

    if os.path.splitext(data_path)[1].lower() == ".h5":
        # .h5 always auto-converts; n_clusters comes from prepare_dataset
        data_path, n_clusters = auto_convert_h5(data_path, dataset_name)
        if n_clusters is None:
            print("ERROR: --n_clusters auto but prepare_dataset failed to infer n_clusters")
            sys.exit(1)
        print(f"  Resolved n_clusters={n_clusters} from .h5 conversion")
    else:
        # .h5ad: n_clusters may be 'auto' or a real int
        if n_clusters_raw == "auto":
            import scanpy as sc
            ad = sc.read_h5ad(data_path)
            # Use resolved_label if available, otherwise try common keys
            label_col = None
            for candidate in ["resolved_label", "cell_type", "Celltype",
                               "celltype", "cell_label", "label"]:
                if candidate in ad.obs.columns:
                    label_col = candidate
                    break
            if label_col is None:
                print(f"ERROR: Cannot auto-detect label column in {data_path}")
                print(f"  Available obs columns: {list(ad.obs.columns)}")
                sys.exit(1)
            n_clusters = len(ad.obs[label_col].unique())
            print(f"  Resolved n_clusters={n_clusters} from --n_clusters auto "
                  f"(label_col={label_col!r})")
        else:
            n_clusters = int(n_clusters_raw)

    # ── Load manifest ───────────────────────────────────────────
    manifest = load_manifest()

    # ── Load runtime registry ──────────────────────────────────
    runtime_registry_path = args.runtime_registry
    if runtime_registry_path is None:
        default_reg = PROJECT_ROOT / "envs" / "runtime_registry.yaml"
        if default_reg.exists():
            runtime_registry_path = str(default_reg)
    if runtime_registry_path:
        runtime_registry = load_runtime_registry(Path(str(runtime_registry_path)))
        print(f"  Runtime registry: {runtime_registry_path} ({len(runtime_registry)} entries)")
    else:
        runtime_registry = {}
        print("  Runtime registry: none (using default sys.executable)")

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
        selected = {}
        for k in DEFAULT_FORMAL_METHODS:
            if k not in manifest:
                print(f"WARNING: Default method key missing from manifest: {k!r}")
                continue
            method_info = manifest[k]
            allowed, reason = check_authenticity(k, method_info)
            if allowed:
                selected[k] = method_info
            else:
                print(f"WARNING: Default method {k!r} skipped ({reason})")

    if not selected:
        print("ERROR: No valid methods selected. Use --methods or check manifest.")
        sys.exit(2)

    # ── Git info ───────────────────────────────────────────────
    commit_sha, branch = get_git_info()

    # ── Summary of what will run ────────────────────────────────
    print("=" * 70)
    print("Formal Benchmark Runner")
    print("=" * 70)
    print(f"  Data:     {data_path}")
    print(f"  Dataset:  {dataset_name}")
    print(f"  Clusters: {n_clusters}")
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
        # Load runtime registry for dry run command display
        runtime_registry_local = load_runtime_registry(Path(str(runtime_registry_path))) if runtime_registry_path else {}
        for method_key, method_info in sorted(selected.items()):
            allowed, reason = check_authenticity(method_key, method_info)
            py_bin = resolve_python_executable(method_key, method_info, runtime_registry_local)
            for seed in args.seeds:
                extra_args = method_info.get("extra_args", None)
                cmd = build_command(
                    py_bin,
                    method_key, method_info, data_path,
                    base_out_dir / f"{method_key}__seed{seed}__dry",
                    n_clusters, args.epochs, args.pretrain_epochs,
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
                    skipped_runtime_env = get_runtime_env(method_info)
                    skipped_python = resolve_python_executable(method_key, method_info, runtime_registry)
                    _write_status(
                        out_dir, method_key, seed, "skipped",
                        None, 0.0, args.gpu, args.no_cuda,
                        datetime.now().isoformat(), datetime.now().isoformat(),
                        [], commit_sha, branch, reason, "",
                        dataset=dataset_name, n_clusters=n_clusters,
                        python_executable=skipped_python,
                        runtime_env=skipped_runtime_env,
                    )
                print()
                continue

            print(f"  Running {method_key}...")
            for seed in args.seeds:
                runtime_env = get_runtime_env(method_info)
                python_bin = resolve_python_executable(method_key, method_info, runtime_registry)
                result = run_method(
                    method_key, method_info,
                    data_path, base_out_dir,
                    n_clusters, args.epochs, args.pretrain_epochs,
                    seed, args.gpu, args.no_cuda,
                    args.allow_unverified,
                    commit_sha, branch,
                    verbose=args.verbose,
                    dry_run=args.dry_run,
                    dataset_name=dataset_name,
                    python_bin=python_bin,
                    runtime_env=runtime_env,
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
