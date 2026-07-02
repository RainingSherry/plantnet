#!/usr/bin/env python3
"""
Model Authenticity Audit Script
================================
Per BDD Scenario 7: Automated model authenticity verification.

This script checks each migrated model for required core components.
It supports two modes:

  1. Card-driven mode (default): reads checks from docs/model_core_cards/*.yaml
     Per BDD Scenario 7: model core definitions must come from cards.

  2. Legacy mode (--legacy): uses hardcoded MODEL_CHECKS.
     Kept for backward compatibility only; new models should use cards.

Scenarios covered:
  - Scenario 7:  core-card-driven audit
  - Scenario 8:  check core losses are used in training loop
  - Scenario 9:  check required_training_stages are preserved
  - Scenario 10: check label leakage
  - Scenario 11: check no OtherMode runtime dependency

Exit codes:
    0  PASS   — All models pass or have documented warnings
    1  FAIL   — At least one model has FAIL status
    2  ERROR  — Script error (file not found, etc.)
"""

import os
import sys
import json
import re
import yaml
from pathlib import Path
from typing import Dict, List, Tuple, Optional

# Project root
PROJECT_ROOT = Path(__file__).parent.parent.resolve()
METHODS_DIR = PROJECT_ROOT / "methods"
MANIFEST_PATH = METHODS_DIR / "method_manifest.yaml"
CARDS_DIR = PROJECT_ROOT / "docs" / "model_core_cards"


# ══════════════════════════════════════════════════════════════════════════════
# LEGACY MODEL_CHECKS — used only with --legacy flag (backward compatibility)
# New models should add entries to docs/model_core_cards/*.yaml
# ══════════════════════════════════════════════════════════════════════════════

LEGACY_MODEL_CHECKS: Dict[str, dict] = {
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
        "forbidden": [("--kmeans_only", "kmeans-only flag")],
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
        "forbidden": [("--kmeans_only", "kmeans-only flag")],
        "auth_patterns": [r"louvain_communities", r"networkx"],
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
        "forbidden": [("--kmeans_only", "kmeans-only flag")],
        "auth_patterns": [
            r"consensus",
            r"KMeans",
            r"AgglomerativeClustering",
        ],
    },
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
        "auth_patterns": [r"attention", r"autoencoder"],
    },
    "scmae": {
        "name": "scMAE",
        "dir": METHODS_DIR / "DeepLearning" / "scMAE",
        "entry": "run.py",
        "required": [
            ("AutoEncoder", "AutoEncoder class"),
            ("apply_scmae_noise", "scMAE mask noise"),
            ("loss_mask", "mask loss computation"),
            ("extract_embedding", "feature extraction"),
            ("KMeans", "KMeans clustering evaluation"),
            ("mask_prob", "mask probability"),
        ],
        "forbidden": [
            ("--fake", "fake flag"),
            ("--dummy", "dummy flag"),
            ("--skip_model", "skip model flag"),
            ("--kmeans_only", "kmeans-only flag"),
        ],
        "auth_patterns": [
            r"class\s+AutoEncoder",
            r"loss_mask",
            r"apply_scmae_noise",
            r"KMeans",
        ],
    },
    "neighbormix_scmae": {
        "name": "NeighborMix_scMAE",
        "dir": PROJECT_ROOT / "experimental_retired_models" / "NeighborMix_scMAE",
        "entry": "run.py",
        "required": [
            ("AutoEncoder", "AutoEncoder class"),
            ("apply_scmae_noise", "scMAE mask noise"),
            ("build_knn_distribution", "KNN distribution construction"),
            ("sample_mix", "neighbor-mix sampling"),
            ("pseudo_branch_enabled", "pseudo-cell branch control"),
            ("loss_mask", "mask loss computation"),
        ],
        "forbidden": [
            ("--fake", "fake flag"),
            ("--dummy", "dummy flag"),
            ("--skip_model", "skip model flag"),
            ("--kmeans_only", "kmeans-only flag"),
        ],
        "auth_patterns": [
            r"class\s+AutoEncoder",
            r"build_knn_distribution",
            r"sample_mix",
            r"loss_mask",
        ],
    },
    "nm_scmae_nomix": {
        "name": "nm_scmae_nomix (ablation)",
        "dir": PROJECT_ROOT / "experimental_retired_models" / "NeighborMix_scMAE",
        "entry": "run.py",
        "required": [
            ("--use_pseudo", "use_pseudo flag"),
            ("--pseudo_weight", "pseudo_weight flag"),
            ("--neighbor_k", "neighbor_k flag"),
            ("--mix_neighbors", "mix_neighbors flag"),
            ("--variant_name", "variant_name flag"),
        ],
        "manifest_check": {
            "key": "nm_scmae_nomix",
            "expected_extra_args": [
                "--use_pseudo", "false",
                "--pseudo_weight", "0",
                "--neighbor_k", "0",
                "--mix_neighbors", "0",
                "--variant_name", "nm_scmae_nomix",
            ],
        },
        "forbidden": [
            ("--fake", "fake flag"),
            ("--dummy", "dummy flag"),
            ("--skip_model", "skip model flag"),
            ("--kmeans_only", "kmeans-only flag"),
        ],
        "auth_patterns": [
            r"--use_pseudo",
            r"--pseudo_weight",
            r"--neighbor_k",
            r"--mix_neighbors",
            r"--variant_name",
        ],
    },
}


# ══════════════════════════════════════════════════════════════════════════════
# Utility functions
# ══════════════════════════════════════════════════════════════════════════════

def load_manifest() -> Dict[str, dict]:
    """Load the method manifest YAML keyed by method key."""
    try:
        with open(MANIFEST_PATH, "r", encoding="utf-8") as f:
            data = yaml.safe_load(f)
        return {m["key"]: m for m in data["methods"]}
    except Exception:
        return {}


def read_file_content(file_path: Path) -> str:
    """Read a file, return empty string if not found."""
    try:
        return file_path.read_text(encoding="utf-8", errors="ignore")
    except FileNotFoundError:
        return ""


def read_model_files(model_dir: Path) -> Dict[str, str]:
    """Read all Python files in a model directory."""
    files = {}
    if not model_dir.exists():
        return files
    for py_file in sorted(model_dir.rglob("*.py")):
        rel = py_file.relative_to(model_dir)
        files[str(rel)] = read_file_content(py_file)
    return files


def check_gpu_default(model_dir: Path, entry: str) -> Tuple[bool, str]:
    """Check --gpu default != 0 (BDD Scenario 13)."""
    entry_path = model_dir / entry
    content = read_file_content(entry_path)
    if not content:
        return True, "entry file not found"
    for line in content.split("\n"):
        if "--gpu" in line.lower():
            if re.search(r"['\"]?\s*--gpu\s*['\"]?\s*,.*?default\s*=\s*0\b", line, re.IGNORECASE):
                return False, f"--gpu default=0: {line.strip()[:100]}"
    return True, "GPU default != 0 or no --gpu argument"


# ══════════════════════════════════════════════════════════════════════════════
# Legacy audit functions (used with --legacy flag)
# ══════════════════════════════════════════════════════════════════════════════

def check_required_keywords(content: str, required: List[Tuple[str, str]]) -> List[Tuple[str, bool]]:
    results = []
    for keyword, description in required:
        found = keyword.lower() in content.lower()
        results.append((description, found))
    return results


def check_forbidden_keywords(content: str, forbidden: List[Tuple[str, str]]) -> List[Tuple[str, bool]]:
    results = []
    for keyword, description in forbidden:
        found = keyword.lower() in content.lower()
        results.append((description, found))
    return results


def check_auth_patterns(content: str, patterns: List[str]) -> List[Tuple[str, bool]]:
    results = []
    for pattern in patterns:
        try:
            found = bool(re.search(pattern, content, re.IGNORECASE))
        except re.error:
            found = False
        results.append((pattern, found))
    return results


def audit_model_legacy(model_key: str, config: dict) -> dict:
    """Legacy keyword-based audit."""
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

    if not model_dir.exists():
        result["status"] = "FAIL"
        result["failures"].append(f"Directory does not exist: {model_dir}")
        return result

    entry_path = model_dir / entry
    if not entry_path.exists():
        result["status"] = "FAIL"
        result["failures"].append(f"Entry file does not exist: {entry_path}")
        return result

    all_files = read_model_files(model_dir)
    if not all_files:
        result["status"] = "FAIL"
        result["failures"].append("No Python files found in directory")
        return result

    combined_content = "\n".join(all_files.values())

    required_results = check_required_keywords(combined_content, config.get("required", []))
    for desc, found in required_results:
        if not found:
            result["status"] = "FAIL"
            result["failures"].append(f"Required component missing: {desc}")

    forbidden_results = check_forbidden_keywords(combined_content, config.get("forbidden", []))
    for desc, found in forbidden_results:
        if found:
            result["status"] = "FAIL"
            result["failures"].append(f"Forbidden pattern found: {desc}")

    pattern_results = check_auth_patterns(combined_content, config.get("auth_patterns", []))
    missing_patterns = [p for p, found in pattern_results if not found]
    if missing_patterns:
        result["warnings"].append(f"Auth patterns not matched: {', '.join(missing_patterns)}")

    gpu_pass, gpu_detail = check_gpu_default(model_dir, entry)
    result["gpu_default_check"] = {"pass": gpu_pass, "detail": gpu_detail}
    if not gpu_pass:
        result["warnings"].append(f"GPU default violation: {gpu_detail}")

    mc = config.get("manifest_check")
    if mc:
        manifest = load_manifest()
        entry_in_manifest = manifest.get(mc.get("key", model_key), {})
        actual_args = entry_in_manifest.get("extra_args", [])
        expected_args = mc.get("expected_extra_args", [])
        if actual_args != expected_args:
            result["status"] = "FAIL"
            result["failures"].append(
                f"manifest extra_args mismatch: got {actual_args!r}, expected {expected_args!r}"
            )

    return result


# ══════════════════════════════════════════════════════════════════════════════
# Card-driven audit (Per BDD Scenarios 7-10)
# Model definitions come from docs/model_core_cards/*.yaml
# ══════════════════════════════════════════════════════════════════════════════

def load_cards() -> Dict[str, dict]:
    """Load all YAML core cards."""
    cards = {}
    if not CARDS_DIR.exists():
        return cards
    for fpath in sorted(CARDS_DIR.glob("*.yaml")):
        try:
            with open(fpath, "r", encoding="utf-8") as f:
                data = yaml.safe_load(f)
            key = data.get("model_key", fpath.stem)
            cards[key] = data
        except Exception as e:
            print(f"  WARNING: Could not parse {fpath}: {e}", file=sys.stderr)
    return cards


def _find_class(content: str, class_name: str) -> bool:
    return bool(re.search(rf"^class\s+{re.escape(class_name)}\s*[\(:]", content, re.MULTILINE))


def _find_method(content: str, method_name: str) -> bool:
    for p in [
        rf"^def\s+{re.escape(method_name)}\s*\(",
        rf"^\s+def\s+{re.escape(method_name)}\s*\(",
        rf"^async\s+def\s+{re.escape(method_name)}\s*\(",
    ]:
        if re.search(p, content, re.MULTILINE):
            return True
    return False


def _loss_in_training(content: str, search_patterns: List[str]) -> bool:
    """
    Check if any of search_patterns appears in a training loop.
    Each pattern is a regex checked (case-insensitive) inside fit/train/main blocks.
    """
    if not search_patterns:
        return True
    compiled = [re.compile(p, re.I) for p in search_patterns]
    lines = content.split("\n")
    in_training = False
    for line in lines:
        stripped = line.strip()
        if re.match(r"^\s*def\s+(fit|train|main)\s*\(", stripped):
            in_training = True
        elif in_training and re.match(r"^\s*def\s+", stripped):
            in_training = False
        if in_training:
            code_part = re.split(r"#", stripped)[0]
            for pattern in compiled:
                if pattern.search(code_part):
                    return True
    return False


def _check_label_leakage(content: str) -> Tuple[List[str], List[str]]:
    """
    Per BDD Scenario 1: distinguish HARD_LABEL_LEAKAGE from SOFT_LABEL_ACCESS.

    HARD (FAIL): ground truth Y used for model/epoch selection, checkpoint selection
                 e.g.  if acc > best: save_model()
    SOFT (WARN): ground truth Y accessed but only for progress printing, final metrics
                 e.g.  print(f"acc={cluster_acc(y, pred)}")
    """
    hard_violations = []
    soft_violations = []
    lines = content.split("\n")
    in_training = False

    # Track if we're inside an "if acc > best:" or "if nmi > best:" block
    # Only assignments inside these blocks count as HARD leakage
    in_acc_if_block = False
    in_nmi_if_block = False
    acc_if_indent = -1
    nmi_if_indent = -1

    def get_indent(line: str) -> int:
        stripped = line.lstrip()
        return len(line) - len(stripped)

    # ── SOFT patterns: Y used for printing / final metrics ───────────────────
    soft_patterns = [
        (r"eval_fn\s*\([^)]*Y\b",              "eval_fn(Y, ...) called in training loop"),
        (r"cluster_acc\s*\([^)]*Y\b",          "cluster_acc(y, ...) called in training loop"),
        (r"nmi_score\s*\([^)]*Y\b",            "nmi_score with Y called in training loop"),
        (r"ari_score\s*\([^)]*Y\b",            "ari_score with Y called in training loop"),
        (r"f1_score\s*\([^)]*Y\b",             "f1_score with Y called in training loop"),
        (r"accuracy_score\s*\([^)]*Y\b",        "accuracy_score with Y called in training loop"),
    ]

    # ── HARD standalone patterns: always hard regardless of context ───────────
    # These are one-liners that ARE the leakage (not inside an if-block)
    standalone_hard = [
        (r"pretrain_acc_max\s*=\s*acc\b",
         "pretrain_acc_max = acc (tracks best accuracy from labels)"),
        (r"acc_max\s*=\s*acc\b",
         "acc_max = acc (tracks best accuracy from labels)"),
        (r"best_nmi\s*=\s*nmi\b",
         "best_nmi = nmi (tracks best NMI from labels)"),
        (r"best_ari\s*=\s*ari\b",
         "best_ari = ari (tracks best ARI from labels)"),
        (r"if\s+acc\s*>=?\s*pretrain_acc_max",
         "if acc >= pretrain_acc_max (early stopping via label metric)"),
    ]

    soft_re = [(re.compile(p, re.I), d) for p, d in soft_patterns]
    standalone_hard_re = [(re.compile(p, re.I), d) for p, d in standalone_hard]

    for line in lines:
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue

        indent = get_indent(line)

        if re.match(r"^\s*def\s+(fit|train|main)\s*\(", stripped):
            in_training = True
            in_acc_if_block = False
            in_nmi_if_block = False
            acc_if_indent = -1
            nmi_if_indent = -1
        elif in_training and re.match(r"^\s*def\s+", stripped):
            in_training = False
            in_acc_if_block = False
            in_nmi_if_block = False
            continue

        if not in_training:
            continue

        code_part = re.split(r"#", stripped)[0]

        # ── Detect acc > best: / nmi > best: conditional blocks ─────────────
        # Enter block: `if ... acc ... > ... best`
        if re.search(r"if\s+.*\bac+[cs]\b.*\>", code_part, re.I):
            in_acc_if_block = True
            acc_if_indent = indent
        elif in_acc_if_block and indent <= acc_if_indent and not stripped.startswith(("elif", "else")):
            in_acc_if_block = False

        if re.search(r"if\s+.*\bnmi\b.*\>", code_part, re.I):
            in_nmi_if_block = True
            nmi_if_indent = indent
        elif in_nmi_if_block and indent <= nmi_if_indent and not stripped.startswith(("elif", "else")):
            in_nmi_if_block = False

        # ── HARD: inside acc > best: or nmi > best: block ──────────────────
        if in_acc_if_block or in_nmi_if_block:
            # Best state saved inside metric-conditional block = HARD leakage
            if re.search(r"best_embedding\s*=|best_y_pred\s*=", code_part):
                hard_violations.append(
                    f"HARD_LABEL_LEAKAGE [Y used: best state saved inside acc/nmi conditional]: {stripped[:100]}"
                )
            # Best model checkpoint saved inside this block = HARD leakage
            if re.search(r"torch\.save|model\.state_dict\(\)|best_model", code_part):
                hard_violations.append(
                    f"HARD_LABEL_LEAKAGE [Y used: model checkpoint saved inside acc/nmi conditional]: {stripped[:100]}"
                )
            continue  # skip soft check for lines inside hard blocks

        # ── SOFT: Y in training loop but not in acc > best: block ─────────
        for pattern, desc in soft_re:
            if pattern.search(code_part):
                soft_violations.append(f"SOFT_LABEL_ACCESS [{desc}]: {stripped[:100]}")
                break

        # ── STANDALONE HARD: patterns that are always hard ─────────────────
        for pattern, desc in standalone_hard_re:
            if pattern.search(code_part):
                hard_violations.append(f"HARD_LABEL_LEAKAGE [{desc}]: {stripped[:100]}")
                break

    return hard_violations, soft_violations


def _check_gpu_in_card(card: dict, target_path: Path) -> Tuple[bool, str]:
    policy = card.get("gpu_policy", "N/A")
    if policy == "N/A":
        return True, "N/A"
    entry = card.get("entry_file", "run.py")
    content = read_file_content(target_path / entry)
    if not content:
        return True, "entry not found"
    for line in content.split("\n"):
        if "--gpu" in line.lower():
            if re.search(r"['\"]?\s*--gpu\s*['\"]?\s*,.*?default\s*=\s*0\b", line, re.IGNORECASE):
                return False, f"--gpu default=0: {line.strip()[:100]}"
    return True, f"Policy={policy}"


def audit_card_driven(model_key: str, card: dict, manifest_entry: dict) -> dict:
    """
    Per BDD Scenarios 7-11: audit a model using its core card.

    Scenario 1:  HARD vs SOFT label leakage
    Scenario 7:  core-card-driven audit
    Scenario 8:  check core losses are used in training loop
    Scenario 9:  check required_training_stages are preserved
    Scenario 10: check label leakage (HARD/FAIL vs SOFT/WARN)
    Scenario 11: check no OtherMode runtime dependency
    """
    target_path_str = card.get("target_path", "")
    target_path = PROJECT_ROOT / target_path_str
    status = card.get("status", "UNKNOWN")
    result = {
        "model_key": model_key,
        "name": card.get("name", model_key),
        "status": "PASS",
        "source": "core_card",
        "card_version": card.get("card_version", "?"),
        "warnings": [],
        "failures": [],
        "info": [],
    }

    # ── Status-based shortcuts ───────────────────────────────────────────
    if status == "CORE-INCOMPLETE":
        result["status"] = "CORE-INCOMPLETE"
        result["info"].append("Placeholder — no core code defined")
        return result

    if status == "ENV-GATED":
        entry_path = target_path / card.get("entry_file", "run.py")
        content = read_file_content(entry_path)
        has_tf = bool(re.search(r"import\s+tensorflow|from\s+tensorflow", content, re.I))
        result["status"] = "ENV-GATED"
        result["info"].append(
            "TensorFlow/Keras dependency detected" if has_tf else "ENV-GATED (stub)"
        )
        return result

    # ── Load model files (including cross-directory siblings) ────────────────
    files = {}
    if target_path.exists():
        for py_file in sorted(target_path.rglob("*.py")):
            rel = str(py_file.relative_to(target_path))
            files[rel] = read_file_content(py_file)
        # Also load cross-directory sibling files referenced in core_source_files
        for rel_file in card.get("core_source_files", []):
            if rel_file.startswith(".."):
                cross_path = target_path / rel_file
                if cross_path.exists():
                    files[rel_file] = read_file_content(cross_path)
    all_content = "\n".join(files.values())

    # 1. Core source files exist
    for rel_file in card.get("core_source_files", []):
        if not (target_path / rel_file).exists():
            result["failures"].append(f"Core file missing: {rel_file}")

    # 2. Core classes present
    for cls in card.get("core_classes", []):
        fname = cls.get("file", "")
        cname = cls.get("name", "")
        content = files.get(fname, all_content)
        if not _find_class(content, cname):
            result["failures"].append(f"Core class '{cname}' not found in {fname or 'any file'}")
        for mp in cls.get("must_preserve", []):
            if isinstance(mp, str) and mp.lower() not in content.lower():
                result["failures"].append(f"  Missing: {mp}")

    # 3. Core functions present
    for fn in card.get("core_functions", []):
        fname = fn.get("file", "")
        fname2 = fn.get("name", "")
        if not fname2:
            sig = fn.get("signature", "")
            fname2 = sig.split("(")[0].replace("def ", "").strip()
        content = files.get(fname, all_content)
        if not _find_method(content, fname2):
            result["failures"].append(f"Core function '{fname2}' not found in {fname or 'any file'}")

    # 4. Core losses used in training (Scenario 3 + 8)
    for loss in card.get("core_losses", []):
        lname = loss.get("name", "")
        stage = loss.get("stage", "")
        # Scenario 3: prefer search_patterns list, fall back to name matching
        search_patterns = loss.get("search_patterns", [])
        if search_patterns:
            if not _loss_in_training(all_content, search_patterns):
                result["warnings"].append(
                    f"Core loss '{lname}' ({stage}) not used in training loop "
                    f"(search_patterns: {search_patterns})"
                )
        elif lname:
            # Fallback: use name as a single pattern
            if not _loss_in_training(all_content, [re.escape(lname)]):
                result["warnings"].append(
                    f"Core loss '{lname}' ({stage}) not used in training loop"
                )

    # 5. Required training stages (Scenario 9)
    for stage in card.get("required_training_stages", []):
        sname = stage.get("name", "")
        evidence = stage.get("evidence", [])
        missing_ev = [ev for ev in evidence if ev.lower() not in all_content.lower()]
        if missing_ev:
            result["warnings"].append(
                f"Training stage '{sname}' may be incomplete: missing [{', '.join(missing_ev)}]"
            )

    # 6. Forbidden changes not introduced
    for forb in card.get("forbidden_changes", []):
        desc = forb if isinstance(forb, str) else forb.get("name", str(forb))
        pattern = forb if isinstance(forb, str) else forb.get("pattern", forb.get("description", ""))
        if pattern and pattern.lower() in all_content.lower():
            result["failures"].append(f"Forbidden change: '{desc}'")

    # 7. Label leakage (BDD Scenario 1 + 10: HARD vs SOFT)
    hard_leaks, soft_leaks = _check_label_leakage(all_content)
    for v in hard_leaks:
        result["failures"].append(v)
    for v in soft_leaks:
        result["warnings"].append(v)

    # 8. GPU policy
    gpu_ok, gpu_detail = _check_gpu_in_card(card, target_path)
    if not gpu_ok:
        result["failures"].append(f"GPU policy violation: {gpu_detail}")

    if result["failures"]:
        result["status"] = "FAIL"

    return result


# ══════════════════════════════════════════════════════════════════════════════
# Main function
# ══════════════════════════════════════════════════════════════════════════════

def main():
    import argparse
    parser = argparse.ArgumentParser(
        description="Audit model authenticity for all migrated methods."
    )
    parser.add_argument("--model", type=str, default=None, help="Audit only a specific model key")
    parser.add_argument("-v", "--verbose", action="store_true", help="Show detailed output")
    parser.add_argument("--json", action="store_true", help="Output results as JSON")
    parser.add_argument(
        "--legacy",
        action="store_true",
        help="Use legacy hardcoded MODEL_CHECKS instead of core cards",
    )
    args = parser.parse_args()

    print("=" * 70)
    print("Model Authenticity Audit")
    print("=" * 70)

    if args.legacy:
        _run_legacy_mode(args)
        return

    # ── Card-driven mode (default) ────────────────────────────────────────
    print("Mode: CARD-DRIVEN (docs/model_core_cards/*.yaml)")
    print()
    manifest = load_manifest()
    cards = load_cards()
    print(f"Loaded {len(cards)} core cards, {len(manifest)} manifest entries.")
    print()

    models_to_check = sorted(set(list(cards.keys()) + list(manifest.keys())))
    results = []
    no_card_results = []  # models without cards (Scenarios 2+4)

    for model_key in models_to_check:
        card = cards.get(model_key)
        manifest_entry = manifest.get(model_key)

        if not card:
            # ── Scenario 2+4: no core card ───────────────────────────────
            manifest_auth = (manifest_entry or {}).get("authenticity", "UNKNOWN")
            is_placeholder = manifest_auth == "PLACEHOLDER"
            # Scenario 2: VERIFIED / default_in_formal without card = FAIL
            is_verified = (manifest_entry or {}).get("authenticity") == "VERIFIED"
            is_default = (manifest_entry or {}).get("default_in_formal") == True
            if is_verified or is_default:
                status = "FAIL"
                icon = "✗"
                reason = "manifest VERIFIED/default_in_formal requires a core card"
                if manifest_entry:
                    mname = manifest_entry.get("name", model_key)
                else:
                    mname = model_key
                print(f"  {icon} [{status:16}] {model_key}")
                if args.verbose:
                    print(f"    FAIL: {reason}")
                no_card_results.append({
                    "model_key": model_key, "name": mname, "status": status,
                    "reason": reason, "manifest_entry": manifest_entry,
                })
            elif is_placeholder:
                # Scenario 4: Foundation PLACEHOLDER → SKIPPED_PLACEHOLDER
                status = "SKIPPED_PLACEHOLDER"
                icon = "○"
                reason = (manifest_entry or {}).get("reason", "Foundation model placeholder")
                mname = (manifest_entry or {}).get("name", model_key)
                print(f"  {icon} [{status:16}] {model_key}")
                if args.verbose:
                    print(f"    INFO: {reason}")
                no_card_results.append({
                    "model_key": model_key, "name": mname, "status": status,
                    "reason": reason, "manifest_entry": manifest_entry,
                })
            else:
                # PENDING without card — informational only
                print(f"  ? [NO_CARD    ] {model_key} — no core card, skipping")
                no_card_results.append({
                    "model_key": model_key, "name": model_key,
                    "status": "NO_CARD", "reason": "no core card",
                    "manifest_entry": manifest_entry,
                })
            print()
            continue

        # ── Normal card-driven audit ─────────────────────────────────────
        print(f"  Checking: {card.get('name', model_key)} ({model_key})...", end=" ", flush=True)
        result = audit_card_driven(model_key, card, manifest_entry)
        results.append(result)

        icon = {
            "PASS": "✓", "FAIL": "✗", "WARN": "⚠",
            "ENV-GATED": "⊗", "CORE-INCOMPLETE": "○",
        }.get(result["status"], "?")

        print(f"{icon} [{result['status']:16}]  (card v{result.get('card_version','?')})")

        if args.verbose:
            for info in result.get("info", []):
                print(f"    INFO: {info}")
            for w in result.get("warnings", []):
                print(f"    WARN: {w}")
            for f in result.get("failures", []):
                print(f"    FAIL: {f}")

    # ── Scenario 2: manifest vs core card conflict detection ─────────────────
    manifest_conflicts = []
    for r in results + no_card_results:
        manifest_entry = r.get("manifest_entry") or {}
        card_status = r.get("status", "")
        manifest_auth = manifest_entry.get("authenticity", "")
        is_default = manifest_entry.get("default_in_formal", False)

        if manifest_auth == "VERIFIED" and card_status == "FAIL":
            manifest_conflicts.append(
                f"  CONFLICT: {r['model_key']} — manifest VERIFIED but core card is FAIL"
            )
        if is_default and card_status == "FAIL":
            manifest_conflicts.append(
                f"  CONFLICT: {r['model_key']} — manifest default_in_formal=true but core card is FAIL"
            )

    # ── Summary ────────────────────────────────────────────────────────────
    print()
    print("=" * 70)
    print("Summary")
    print("=" * 70)

    all_results = results + no_card_results
    counts = {}
    for r in all_results:
        s = r.get("status", "UNKNOWN")
        counts[s] = counts.get(s, 0) + 1
    total = len(all_results)
    for s, c in sorted(counts.items()):
        print(f"  {s:25}: {c}/{total}")

    if manifest_conflicts:
        print()
        print("  MANIFEST CONFLICTS (BDD Scenario 2):")
        for c in manifest_conflicts:
            print(c)

    if args.json:
        print()
        print(json.dumps({
            "mode": "card-driven",
            "results": results,
            "no_card_results": no_card_results,
            "manifest_conflicts": manifest_conflicts,
            "summary": counts,
        }, indent=2, default=str))

    # Exit code: 0 = all pass, 1 = any FAIL or manifest conflict
    has_fail = any(
        r.get("status") == "FAIL"
        for r in results + no_card_results
    )
    if has_fail or manifest_conflicts:
        print()
        print("  Some models have authenticity failures.")
        sys.exit(1)
    sys.exit(0)


def _run_legacy_mode(args):
    """Legacy keyword-based audit (backward compatibility)."""
    print("Mode: LEGACY (hardcoded MODEL_CHECKS)")
    print()
    print("NOTE: Prefer default card-driven mode. Core card definitions live in")
    print("      docs/model_core_cards/*.yaml")
    print()

    manifest = load_manifest()
    results = []
    any_fail = False

    for model_key, config in sorted(LEGACY_MODEL_CHECKS.items()):
        print(f"  Auditing: {config['name']} ({model_key})...", end=" ", flush=True)

        entry_path = config["dir"] / config.get("entry", "run.py")
        entry_content = read_file_content(entry_path)

        is_env_gated = False
        for tf_pattern in [
            r"import\s+tensorflow", r"from\s+tensorflow",
            r"import\s+keras", r"from\s+keras",
        ]:
            if re.search(tf_pattern, entry_content, re.IGNORECASE):
                is_env_gated = True
                break

        if is_env_gated:
            result = {
                "model_key": model_key,
                "name": config["name"],
                "status": "ENV-GATED",
                "detail": "TensorFlow/Keras dependency detected",
                "warnings": [],
                "failures": [],
            }
            print(f"⊗ [ENV-GATED]")
        else:
            result = audit_model_legacy(model_key, config)
            icon = "✓" if result["status"] == "PASS" else "✗"
            extra = f" ({len(result['warnings'])} warning(s))" if result.get("warnings") else ""
            print(f"{icon} [{result['status']:9}]{extra}")

            if args.verbose:
                for f in result.get("failures", []):
                    print(f"    FAIL: {f}")
                for w in result.get("warnings", []):
                    print(f"    WARN: {w}")
                if result.get("gpu_default_check"):
                    g = result["gpu_default_check"]
                    print(f"    GPU:  {'PASS' if g['pass'] else 'FAIL'} — {g['detail']}")

            if result["status"] == "FAIL":
                any_fail = True

        results.append(result)

    print()
    print("=" * 70)
    print("Summary")
    print("=" * 70)
    counts = {"PASS": 0, "FAIL": 0, "ENV-GATED": 0}
    for r in results:
        s = r["status"]
        counts[s] = counts.get(s, 0) + 1
    total = len(results)
    for s, c in sorted(counts.items()):
        print(f"  {s:20}: {c}/{total}")

    if args.json:
        print()
        print(json.dumps({"mode": "legacy", "results": results, "summary": counts}, indent=2, default=str))

    sys.exit(1 if any_fail else 0)


if __name__ == "__main__":
    main()
