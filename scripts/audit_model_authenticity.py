#!/usr/bin/env python3
"""
Model Authenticity Audit Script
================================
Per BDD Scenario 11: Automated model authenticity verification.

This script checks each migrated model for required core components
without understanding model semantics — using keyword and structural checks.

Exit codes:
    0  PASS   — All models pass or have documented warnings
    1  FAIL   — At least one model has FAIL status
    2  ERROR  — Script error (file not found, etc.)
"""

import os
import sys
import json
import re
from pathlib import Path
from typing import Dict, List, Tuple, Optional

# Project root
PROJECT_ROOT = Path(__file__).parent.parent.resolve()
METHODS_DIR = PROJECT_ROOT / "methods"

# ─────────────────────────────────────────────────────────────
# Model definitions: each entry specifies required keyword checks
# Format: { model_key: { "dir": ..., "required": [...], "forbidden": [...] } }
# ─────────────────────────────────────────────────────────────

MODEL_CHECKS: Dict[str, dict] = {
    # ── Deep Learning ────────────────────────────────────────
    "dec": {
        "name": "DEC",
        "dir": METHODS_DIR / "DeepLearning" / "dec",
        "entry": "run.py",
        "required": [
            ("Autoencoder", "Autoencoder class"),
            ("ClusteringLayer", "ClusteringLayer class"),
            ("target_distribution", "target_distribution function"),
            ("kl_div", "KL divergence loss"),
            ("KMeans", "KMeans initialization"),
            ("pretrain", "pretrain method"),
            ("clustering_layer", "clustering layer usage"),
        ],
        "forbidden": [
            ("--fake", "fake flag"),
            ("--dummy", "dummy flag"),
            ("--skip_model", "skip model flag"),
            ("--kmeans_only", "kmeans-only flag"),
        ],
        "auth_patterns": [
            r"class\s+Autoencoder",
            r"class\s+ClusteringLayer",
            r"def\s+target_distribution",
            r"kl_div",
        ],
    },

    "scdcc": {
        "name": "scDCC",
        "dir": METHODS_DIR / "DeepLearning" / "scDCC",
        "entry": "run.py",
        "required": [
            ("ZINBLoss", "ZINBLoss class"),
            ("MeanAct", "MeanAct activation"),
            ("DispAct", "DispAct activation"),
            ("soft_assign", "soft_assign method"),
            ("target_distribution", "target_distribution method"),
            ("_dec_mean", "ZINB decoder mean"),
            ("_dec_disp", "ZINB decoder dispersion"),
            ("_dec_pi", "ZINB decoder pi"),
        ],
        "forbidden": [
            ("--fake", "fake flag"),
            ("--dummy", "dummy flag"),
            ("--skip_model", "skip model flag"),
            ("--kmeans_only", "kmeans-only flag"),
        ],
        "auth_patterns": [
            r"ZINBLoss",
            r"MeanAct",
            r"DispAct",
            r"soft_assign",
            r"target_distribution",
            r"_dec_mean",
            r"_dec_disp",
            r"_dec_pi",
        ],
    },

    "scdsc": {
        "name": "scDSC (SDCN)",
        "dir": METHODS_DIR / "GNN" / "scDSC",
        "entry": "run.py",
        "required": [
            ("GNNLayer", "GNNLayer class"),
            ("ZINBLoss", "ZINBLoss"),
            ("cluster_layer", "cluster layer"),
            ("target_distribution", "target_distribution"),
            ("gnn_", "GNN layers (gnn_1, gnn_2, ...)"),
            ("SDCN", "SDCN model class"),
            ("_dec_mean", "ZINB decoder mean"),
        ],
        "forbidden": [
            ("--fake", "fake flag"),
            ("--dummy", "dummy flag"),
            ("--skip_model", "skip model flag"),
            ("--kmeans_only", "kmeans-only flag"),
        ],
        "auth_patterns": [
            r"GNNLayer",
            r"ZINBLoss",
            r"cluster_layer",
            r"target_distribution",
            r"SDCN",
            r"self\.gnn_\d",
        ],
    },

    "scdeepcluster": {
        "name": "scDeepCluster",
        "dir": METHODS_DIR / "DeepLearning" / "scDeepCluster",
        "entry": "run.py",
        "required": [
            ("SCDeepCluster", "SCDeepCluster class"),
            ("pretrain", "pretrain method"),
            ("fit", "fit method"),
            ("extract_feature", "extract_feature method"),
            ("code/", "code/ directory migrated"),
        ],
        "forbidden": [
            ("--fake", "fake flag"),
            ("--dummy", "dummy flag"),
            ("--kmeans_only", "kmeans-only flag"),
        ],
        "auth_patterns": [
            r"SCDeepCluster",
            r"pretrain",
            r"fit",
            r"extract_feature",
            r"code",
        ],
    },

    # ── Traditional ─────────────────────────────────────────
    "scanpy_standard": {
        "name": "ScanpyStandard",
        "dir": METHODS_DIR / "Traditional" / "ScanpyStandard",
        "entry": "run.py",
        "required": [
            ("normalize_total", "normalize_total"),
            ("highly_variable", "highly_variable_genes"),
            ("sc.tl.pca", "PCA"),
            ("sc.pp.neighbors", "neighbors graph"),
            ("sc.tl.umap", "UMAP"),
            ("sc.tl.leiden", "Leiden clustering"),
        ],
        "forbidden": [],
        "auth_patterns": [
            r"normalize_total",
            r"highly_variable",
            r"sc\.tl\.pca",
            r"sc\.pp\.neighbors",
            r"sc\.tl\.umap",
            r"sc\.tl\.leiden",
        ],
    },

    "leiden": {
        "name": "Leiden",
        "dir": METHODS_DIR / "Traditional" / "Leiden",
        "entry": "run.py",
        "required": [
            ("leidenalg", "leidenalg import"),
            ("leiden", "leiden function call"),
            ("RBConfigurationVertexPartition", "RBConfigurationVertexPartition"),
            ("igraph", "igraph import"),
        ],
        "forbidden": [
            ("--kmeans_only", "kmeans-only flag"),
        ],
        "auth_patterns": [
            r"leidenalg",
            r"find_partition",
            r"RBConfigurationVertexPartition",
            r"igraph",
        ],
    },

    "louvain": {
        "name": "Louvain",
        "dir": METHODS_DIR / "Traditional" / "Louvain",
        "entry": "run.py",
        "required": [
            ("louvain_communities", "louvain_communities"),
            ("networkx", "networkx import"),
        ],
        "forbidden": [
            ("--kmeans_only", "kmeans-only flag"),
        ],
        "auth_patterns": [
            r"louvain_communities",
            r"networkx",
        ],
    },

    "sc3": {
        "name": "sc3",
        "dir": METHODS_DIR / "Traditional" / "sc3",
        "entry": "run.py",
        "required": [
            ("consensus", "consensus matrix"),
            ("KMeans", "KMeans ensemble"),
            ("AgglomerativeClustering", "hierarchical clustering"),
        ],
        "forbidden": [
            ("--kmeans_only", "kmeans-only flag"),
        ],
        "auth_patterns": [
            r"consensus",
            r"KMeans",
            r"AgglomerativeClustering",
        ],
    },

    # ── GNN ─────────────────────────────────────────────────
    "scgnn": {
        "name": "scGNN",
        "dir": METHODS_DIR / "GNN" / "scGNN",
        "entry": "run.py",
        "required": [
            ("scGNN", "scGNN model"),
            ("graph", "graph construction"),
            ("cluster", "clustering module"),
            ("AE", "AE model component"),
            ("GAE", "GAE component"),
        ],
        "forbidden": [
            ("--fake", "fake flag"),
            ("--dummy", "dummy flag"),
        ],
        "auth_patterns": [
            r"scGNN",
            r"graph",
            r"cluster",
            r"autoencoder",
            r"gae",
        ],
    },

    "sccdcg": {
        "name": "scCDCG",
        "dir": METHODS_DIR / "GNN" / "scCDCG",
        "entry": "run.py",
        "required": [
            ("scCDCG", "scCDCG model"),
            ("AE", "autoencoder component"),
            ("cluster", "clustering module"),
            ("AE_GAT", "GAT-based autoencoder"),
            ("AE_NN", "NN-based autoencoder"),
        ],
        "forbidden": [
            ("--fake", "fake flag"),
            ("--dummy", "dummy flag"),
        ],
        "auth_patterns": [
            r"scCDCG",
            r"AE_GAT",
            r"AE_NN",
            r"ClusterAssignment",
            r"pdf_norm",
        ],
    },

    "attentionae_sc": {
        "name": "AttentionAE_sc",
        "dir": METHODS_DIR / "GNN" / "AttentionAE_sc",
        "entry": "run.py",
        "required": [
            ("attention", "attention mechanism"),
            ("autoencoder", "autoencoder model"),
        ],
        "forbidden": [
            ("--fake", "fake flag"),
            ("--dummy", "dummy flag"),
        ],
        "auth_patterns": [
            r"attention",
            r"autoencoder",
        ],
    },
}


def read_file_content(file_path: Path) -> str:
    """Read file content, return empty string if not found."""
    try:
        return file_path.read_text(encoding="utf-8", errors="ignore")
    except FileNotFoundError:
        return ""


def read_model_files(model_dir: Path, entry: str) -> Dict[str, str]:
    """Read all Python files in a model directory."""
    files = {}
    if not model_dir.exists():
        return files

    for py_file in model_dir.rglob("*.py"):
        rel = py_file.relative_to(model_dir)
        files[str(rel)] = read_file_content(py_file)

    return files


def check_required_keywords(content: str, required: List[Tuple[str, str]]) -> List[Tuple[str, bool]]:
    """Check that all required keywords appear in the content."""
    results = []
    for keyword, description in required:
        found = keyword.lower() in content.lower()
        results.append((description, found))
    return results


def check_forbidden_keywords(content: str, forbidden: List[Tuple[str, str]]) -> List[Tuple[str, bool]]:
    """Check that no forbidden patterns appear (True = violation found)."""
    results = []
    for keyword, description in forbidden:
        # Case-insensitive search
        found = keyword.lower() in content.lower()
        results.append((description, found))
    return results


def check_auth_patterns(content: str, patterns: List[str]) -> List[Tuple[str, bool]]:
    """Check regex patterns for structural components."""
    results = []
    for pattern in patterns:
        try:
            found = bool(re.search(pattern, content, re.IGNORECASE))
        except re.error:
            found = False
        results.append((pattern, found))
    return results


def check_gpu_default(model_dir: Path, entry: str) -> Tuple[bool, str]:
    """
    Check that --gpu does not default to 0.
    Returns (pass, detail).
    """
    entry_path = model_dir / entry
    content = read_file_content(entry_path)
    if not content:
        return True, "entry file not found, skipping GPU check"

    lines = content.split('\n')
    for line in lines:
        # Only match the same logical line (not across newlines)
        if '--gpu' in line.lower() or "'--gpu'" in line or '"--gpu"' in line:
            # Check if default=0 appears on the same line
            if re.search(r"['\"]?\s*--gpu\s*['\"]?\s*,.*?default\s*=\s*0\b", line, re.IGNORECASE):
                return False, f"--gpu default=0 found on line: {line.strip()[:120]}"
    return True, "GPU default != 0 or no --gpu argument"


def audit_model(model_key: str, config: dict) -> dict:
    """Audit a single model. Returns a dict with audit results."""
    result = {
        "model_key": model_key,
        "name": config["name"],
        "dir": str(config["dir"]),
        "status": "PASS",
        "warnings": [],
        "failures": [],
        "gpu_default_check": None,
        "missing_files": [],
    }

    model_dir = config["dir"]
    entry = config.get("entry", "run.py")

    # Check if directory exists
    if not model_dir.exists():
        result["status"] = "FAIL"
        result["failures"].append(f"Directory does not exist: {model_dir}")
        return result

    # Check if entry file exists
    entry_path = model_dir / entry
    if not entry_path.exists():
        result["status"] = "FAIL"
        result["failures"].append(f"Entry file does not exist: {entry_path}")
        return result

    # Read all Python files
    all_files = read_model_files(model_dir, entry)
    if not all_files:
        result["status"] = "FAIL"
        result["failures"].append("No Python files found in directory")
        return result

    # Combine all content for keyword checks
    combined_content = "\n".join(all_files.values())

    # Check required keywords
    required_results = check_required_keywords(combined_content, config.get("required", []))
    for desc, found in required_results:
        if not found:
            result["status"] = "FAIL"
            result["failures"].append(f"Required component missing: {desc}")

    # Check forbidden keywords
    forbidden_results = check_forbidden_keywords(combined_content, config.get("forbidden", []))
    for desc, found in forbidden_results:
        if found:
            result["status"] = "FAIL"
            result["failures"].append(f"Forbidden pattern found: {desc}")

    # Check auth regex patterns
    pattern_results = check_auth_patterns(combined_content, config.get("auth_patterns", []))
    missing_patterns = [p for p, found in pattern_results if not found]
    if missing_patterns:
        result["warnings"].append(f"Auth patterns not matched: {', '.join(missing_patterns)}")

    # Check GPU default
    gpu_pass, gpu_detail = check_gpu_default(model_dir, entry)
    result["gpu_default_check"] = {"pass": gpu_pass, "detail": gpu_detail}
    if not gpu_pass:
        result["warnings"].append(f"GPU default violation: {gpu_detail}")

    return result


def check_env_blocked(model_key: str, config: dict) -> dict:
    """Check if a model is environment-blocked (e.g., TensorFlow)."""
    result = {
        "model_key": model_key,
        "name": config["name"],
        "status": "ENV-GATED",
        "reason": "Environment-gated model",
        "detail": "",
    }

    model_dir = config["dir"]
    entry = model_dir / config.get("entry", "run.py")
    content = read_file_content(entry_path) if (entry_path := entry).exists() else ""

    # Check for TF/Keras import patterns
    tf_patterns = [
        r"import\s+tensorflow",
        r"from\s+tensorflow",
        r"import\s+keras",
        r"from\s+keras",
    ]
    for pattern in tf_patterns:
        if re.search(pattern, content, re.IGNORECASE):
            result["detail"] = "TensorFlow/Keras dependency detected"
            break

    return result


def main():
    import argparse

    parser = argparse.ArgumentParser(
        description="Audit model authenticity for all migrated methods."
    )
    parser.add_argument(
        "--model",
        type=str,
        default=None,
        help="Audit only a specific model key (e.g., 'dec', 'scdcc'). "
             "If not provided, audits all models.",
    )
    parser.add_argument(
        "--verbose", "-v",
        action="store_true",
        help="Show detailed output for each check",
    )
    parser.add_argument(
        "--json",
        action="store_true",
        help="Output results as JSON",
    )
    args = parser.parse_args()

    # Determine which models to audit
    if args.model:
        if args.model not in MODEL_CHECKS:
            print(f"ERROR: Unknown model key '{args.model}'")
            print(f"Available: {', '.join(sorted(MODEL_CHECKS.keys()))}")
            sys.exit(2)
        models_to_audit = {args.model: MODEL_CHECKS[args.model]}
    else:
        models_to_audit = MODEL_CHECKS

    results = []
    all_pass = True
    any_fail = False

    print("=" * 70)
    print("Model Authenticity Audit")
    print("=" * 70)
    print()

    for model_key, config in sorted(models_to_audit.items()):
        print(f"  Auditing: {config['name']} ({model_key})...")

        # Quick env-gated check: look for TF in the main entry file
        entry_path = config["dir"] / config.get("entry", "run.py")
        entry_content = read_file_content(entry_path)

        is_env_gated = False
        for tf_pattern in [r"import\s+tensorflow", r"from\s+tensorflow", r"import\s+keras", r"from\s+keras"]:
            if re.search(tf_pattern, entry_content, re.IGNORECASE):
                is_env_gated = True
                break

        if is_env_gated:
            print(f"    → ENV-GATED (TensorFlow/Keras dependency)")
            result = check_env_blocked(model_key, config)
        else:
            result = audit_model(model_key, config)

        results.append(result)

        status_icon = {"PASS": "✓", "FAIL": "✗", "ENV-GATED": "⊗", "WARN": "⚠"}.get(
            result["status"], "?"
        )
        status_str = f"  {status_icon} [{result['status']}] {config['name']}"
        if result.get("warnings"):
            status_str += f"  ({len(result['warnings'])} warning(s))"
        print(status_str)

        if args.verbose:
            if result.get("failures"):
                for f in result["failures"]:
                    print(f"      FAIL: {f}")
            if result.get("warnings"):
                for w in result["warnings"]:
                    print(f"      WARN: {w}")
            if result.get("gpu_default_check"):
                g = result["gpu_default_check"]
                print(f"      GPU:  {'PASS' if g['pass'] else 'FAIL'} — {g['detail']}")

        if result["status"] == "FAIL":
            any_fail = True
            all_pass = False
        elif result["status"] == "ENV-GATED":
            all_pass = False  # Env-gated is not a pass

        print()

    # ── Summary ───────────────────────────────────────────────
    print("=" * 70)
    print("Summary")
    print("=" * 70)

    status_counts = {"PASS": 0, "FAIL": 0, "ENV-GATED": 0, "WARN": 0}
    for r in results:
        s = r["status"]
        status_counts[s] = status_counts.get(s, 0) + 1

    total = len(results)
    print(f"  Total models audited: {total}")
    print(f"  PASS:        {status_counts.get('PASS', 0)}/{total}")
    print(f"  FAIL:        {status_counts.get('FAIL', 0)}/{total}")
    print(f"  ENV-GATED:   {status_counts.get('ENV-GATED', 0)}/{total}")

    if any_fail:
        print()
        print("  ⚠  Some models have authenticity issues.")
        print("     Models with FAIL status cannot enter the formal benchmark.")
        print("     Models with ENV-GATED status are blocked by environment dependencies.")

    if args.json:
        print()
        print(json.dumps({"results": results, "summary": status_counts}, indent=2, default=str))

    # ── Exit code ────────────────────────────────────────────
    if any_fail:
        sys.exit(1)
    else:
        sys.exit(0)


if __name__ == "__main__":
    main()
