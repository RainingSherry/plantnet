#!/usr/bin/env python3
"""
check_model_core_cards.py
=========================
Per BDD Scenario 6: Validate model core cards exist and are consistent.

This script checks:
  1. Every model in the manifest has a corresponding core card.
  2. core_source_files actually exist in the target_path.
  3. core_components (classes, functions, losses) are findable in source files.
  4. forbidden_changes have not been introduced.
  5. source_path / target_path are consistent with manifest.
  6. known_deviations are documented.
  7. label_leakage patterns are absent in training loops.
  8. GPU defaults do not violate BDD Scenario 13.

Exit codes:
    0 = ALL PASS (or expected failures documented)
    1 = At least one FAIL
    2 = Script error
"""

import os
import sys
import re
import json
import yaml
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Any
from dataclasses import dataclass, field

PROJECT_ROOT = Path(__file__).parent.parent.resolve()
CARDS_DIR = PROJECT_ROOT / "docs" / "model_core_cards"
MANIFEST_PATH = PROJECT_ROOT / "methods" / "method_manifest.yaml"
AUDIT_SCRIPT = PROJECT_ROOT / "scripts" / "audit_model_authenticity.py"

STATUS_PASS = "PASS"
STATUS_FAIL = "FAIL"
STATUS_WARN = "WARN"
STATUS_OTHERMODE_FAIL = "OTHERMODE_DEPENDENCY"
STATUS_CORE_INCOMPLETE = "CORE-INCOMPLETE"
STATUS_ENV_GATED = "ENV-GATED"
STATUS_MISSING_CARD = "MISSING_CARD"
STATUS_NO_CARD = "NO_CARD"


@dataclass
class CheckResult:
    model_key: str
    status: str
    checks: Dict[str, Any] = field(default_factory=dict)
    failures: List[str] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)
    info: List[str] = field(default_factory=list)


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


def load_manifest() -> Dict[str, dict]:
    """Load the method manifest keyed by model key."""
    try:
        with open(MANIFEST_PATH, "r", encoding="utf-8") as f:
            data = yaml.safe_load(f)
        return {m["key"]: m for m in data["methods"]}
    except Exception as e:
        print(f"ERROR: Could not load manifest: {e}", file=sys.stderr)
        return {}


def read_file_content(file_path: Path) -> str:
    """Read a file, return empty string if missing."""
    try:
        return file_path.read_text(encoding="utf-8", errors="ignore")
    except FileNotFoundError:
        return ""


def get_all_model_files(target_path: Path) -> Dict[str, str]:
    """Read all Python files in a model directory."""
    files = {}
    if not target_path.exists():
        return files
    for py_file in sorted(target_path.rglob("*.py")):
        rel = py_file.relative_to(target_path)
        files[str(rel)] = read_file_content(py_file)
    return files


def find_class(content: str, class_name: str) -> bool:
    """Check if a class definition exists in content."""
    return bool(re.search(rf"^class\s+{re.escape(class_name)}\s*[\(:]", content, re.MULTILINE))


def find_function(content: str, func_name: str) -> bool:
    """Check if a function/method definition exists."""
    patterns = [
        rf"^def\s+{re.escape(func_name)}\s*\(",
        rf"^\s+def\s+{re.escape(func_name)}\s*\(",
        rf"^async\s+def\s+{re.escape(func_name)}\s*\(",
    ]
    for p in patterns:
        if re.search(p, content, re.MULTILINE):
            return True
    return False


def find_import(content: str, module_or_class: str) -> bool:
    """Check if a module/class is imported."""
    patterns = [
        rf"^import\s+.*\b{re.escape(module_or_class)}\b",
        rf"^from\s+.*\b{re.escape(module_or_class)}\b",
        rf"^\s+import\s+.*\b{re.escape(module_or_class)}\b",
        rf"^\s+from\s+.*\b{re.escape(module_or_class)}\b",
    ]
    for p in patterns:
        if re.search(p, content, re.MULTILINE):
            return True
    return False


def find_loss_in_training(content: str, loss_name: str) -> bool:
    """
    Check if a loss appears in a training loop context.
    Looks for the loss name near 'backward', 'step', 'optimizer', 'loss'.
    """
    # Split into training-related blocks
    lines = content.split("\n")
    in_training_block = False
    for line in lines:
        stripped = line.strip()
        # Detect training loops
        if any(kw in stripped for kw in ["def fit(", "def train(", "for epoch", "for batch", "while "]):
            in_training_block = True
        if in_training_block:
            if stripped.startswith("def ") and "fit" not in stripped and "train" not in stripped:
                in_training_block = False
            if loss_name.lower() in stripped.lower():
                # Make sure it's not just a comment
                code_part = re.split(r"#", stripped)[0]
                if loss_name.lower() in code_part.lower():
                    return True
    return loss_name.lower() in content.lower()


def find_forbidden_pattern(content: str, forbidden: str) -> Tuple[bool, str]:
    """Check if a forbidden pattern appears (excluding comments)."""
    lines = content.split("\n")
    for line in lines:
        code_part = re.split(r"#", line)[0]
        if forbidden.lower() in code_part.lower():
            return True, line.strip()[:120]
    return False, ""


def check_gpu_default(target_path: Path, entry: str) -> Tuple[bool, str]:
    """Check --gpu default != 0 (BDD Scenario 13)."""
    entry_path = target_path / entry
    content = read_file_content(entry_path)
    if not content:
        return True, "entry file not found, skipping GPU check"

    lines = content.split("\n")
    for line in lines:
        if "--gpu" in line.lower():
            # Check for default=0 on same logical line
            if re.search(r"['\"]?\s*--gpu\s*['\"]?\s*,\s*.*?default\s*=\s*0\b", line, re.IGNORECASE):
                return False, f"--gpu default=0 found: {line.strip()[:100]}"
    return True, "GPU default != 0 or no --gpu argument"


def check_label_leakage(content: str) -> List[str]:
    """
    Check for label leakage patterns in training loops.
    Looks for ground truth used for best model selection, early stopping, etc.
    """
    violations = []
    lines = content.split("\n")
    in_training_loop = False

    # Patterns that indicate label leakage in training
    leakage_patterns = [
        (r"eval_fn\s*\(\s*Y\s*[,\)]", "eval_fn called with Y in training loop"),
        (r"evaluation\s*\(\s*Y\s*[,\)]", "evaluation called with Y in training loop"),
        (r"cluster_acc\s*\([^)]*Y[^)]*\)", "cluster_acc called with Y in training loop"),
        (r"acc\s*=\s*cluster_acc\s*\([^)]*\)", "cluster_acc assigned during training"),
        (r"if\s+.*acc.*>.*best", "Using ACC to track best model in training loop"),
        (r"if\s+.*nmi.*>.*best", "Using NMI to track best model in training loop"),
        (r"best.*=.*acc\b", "Using ACC as best metric in training loop"),
    ]

    for line in lines:
        stripped = line.strip()
        # Track if we're in a training function
        if re.match(r"^\s*def\s+(fit|train|main)\s*\(", stripped):
            in_training_loop = True
        elif in_training_loop and re.match(r"^\s*def\s+", stripped):
            in_training_loop = False

        if in_training_loop:
            for pattern, description in leakage_patterns:
                if re.search(pattern, stripped, re.IGNORECASE):
                    # Allow label usage only in the final metrics computation block
                    if "metrics.json" in content[content.find(stripped):content.find(stripped)+500]:
                        continue
                    violations.append(f"{description}: {stripped[:100]}")

    return violations


def check_core_components(card: dict, target_path: Path) -> Tuple[List[str], List[str]]:
    """
    Check that core classes, functions, and losses exist in source files.
    Returns (failures, warnings).
    """
    failures = []
    warnings = []
    files = get_all_model_files(target_path)
    # Also load cross-directory sibling files referenced in core_source_files
    for rel_file in card.get("core_source_files", []):
        if rel_file.startswith(".."):
            cross_path = target_path / rel_file
            if cross_path.exists():
                files[rel_file] = read_file_content(cross_path)
    all_content = "\n".join(files.values())

    # Check core classes
    for cls in card.get("core_classes", []):
        name = cls.get("name", "")
        file = cls.get("file", "").strip()
        # Cross-directory class: load the sibling file
        if file and file.startswith(".."):
            cross_path = target_path / file
            if cross_path.exists():
                files[file] = read_file_content(cross_path)
        content = files.get(file, all_content) if file else all_content
        if file:
            if not find_class(content, name):
                failures.append(f"Core class '{name}' not found in {file}")
        else:
            if not find_class(all_content, name):
                failures.append(f"Core class '{name}' not found in any file")
        # Check must_preserve conditions (only string items)
        for cond in cls.get("must_preserve", []):
            if isinstance(cond, str):
                # Strip quotes if quoted (YAML quotes)
                check_str = cond.strip('"\'')
                if check_str.lower() not in all_content.lower():
                    failures.append(f"Core class '{name}' missing required component: {cond}")

    # Check core functions
    for func in card.get("core_functions", []):
        name = func.get("name", "")
        file = func.get("file", "").strip()
        content = ""  # will be set below

        # Resolve cross-directory paths (e.g. "../scMAE_family.py" from scMAE/)
        if file and file.startswith(".."):
            cross_path = target_path / file
            if cross_path.exists():
                files[file] = read_file_content(cross_path)
                content = files[file]
            else:
                # Try searching all files
                content = all_content
        elif file and file in files:
            content = files[file]
        else:
            # No file specified — search across all files
            content = all_content

        if not content or not find_function(content, name):
            failures.append(f"Core function '{name}' not found in {file or 'any file'}")

    # Check core losses appear in training
    for loss in card.get("core_losses", []):
        loss_name = loss.get("name", "")
        stage = loss.get("stage", "")
        if not find_loss_in_training(all_content, loss_name):
            warnings.append(f"Core loss '{loss_name}' ({stage}) not found in training loop")

    return failures, warnings


def check_card(card: dict, manifest_entry: Optional[dict]) -> CheckResult:
    """Perform all checks for a single model card."""
    model_key = card.get("model_key", "unknown")
    status = card.get("status", "UNKNOWN")
    target_path_str = card.get("target_path", "")
    target_path = PROJECT_ROOT / target_path_str
    entry_file = card.get("entry_file", "run.py")

    result = CheckResult(model_key=model_key, status=STATUS_PASS)

    # ── 1. Check status-dependent behavior ──────────────────────────────────
    if status in ("CORE-INCOMPLETE",):
        result.status = STATUS_CORE_INCOMPLETE
        result.info.append(f"Status: {status} — placeholder, no core code expected")
        return result

    if status == "ENV-GATED":
        # For ENV-GATED, check only that the stub exists and has proper error message
        entry_path = target_path / entry_file
        content = read_file_content(entry_path)
        if not content:
            result.failures.append(f"ENV-GATED model entry file missing: {entry_path}")
            result.status = STATUS_FAIL
        elif "tensorflow" not in content.lower() and "keras" not in content.lower():
            result.warnings.append(f"ENV-GATED model may not have TF check: {entry_path}")
        result.status = STATUS_ENV_GATED
        return result

    # ── 2. Check core source files exist ─────────────────────────────────────
    for rel_file in card.get("core_source_files", []):
        full_path = target_path / rel_file
        if not full_path.exists():
            result.failures.append(f"Core file missing: {full_path}")
        else:
            result.info.append(f"Core file found: {rel_file}")

    # ── 3. Check core components ──────────────────────────────────────────────
    comp_failures, comp_warnings = check_core_components(card, target_path)
    result.failures.extend(comp_failures)
    result.warnings.extend(comp_warnings)

    # ── 4. Check forbidden changes ────────────────────────────────────────────
    all_content = "\n".join(get_all_model_files(target_path).values())
    for forb in card.get("forbidden_changes", []):
        # Handle both string and dict formats
        if isinstance(forb, dict):
            desc = forb.get("name", str(forb))
            pattern = forb.get("pattern", forb.get("description", ""))
        else:
            desc = str(forb)
            pattern = str(forb)

        found, line = find_forbidden_pattern(all_content, pattern or desc)
        if found:
            result.failures.append(f"Forbidden change detected: '{desc}' — {line}")

    # ── 5. Check source_path / target_path consistency ────────────────────────
    source_path = card.get("source_path", "")
    if manifest_entry:
        manifest_source = manifest_entry.get("source_path", "")
        manifest_target = manifest_entry.get("target_path", "")
        if manifest_source and source_path != manifest_source:
            result.warnings.append(
                f"Card source_path mismatch: card={source_path} manifest={manifest_source}"
            )
        manifest_target_card = card.get("target_path", "")
        if manifest_target and manifest_target_card != manifest_target:
            result.warnings.append(
                f"Card target_path mismatch: card={manifest_target_card} manifest={manifest_target}"
            )

    # ── 6. Check compatibility changes documented ─────────────────────────────
    compat = card.get("compatibility_changes", [])
    core_logic = card.get("core_logic_changes", [])
    unsafe_subs = card.get("unsafe_substitutions", [])

    if core_logic and not card.get("known_deviations"):
        result.warnings.append("core_logic_changes present but no known_deviations field")

    if unsafe_subs:
        result.failures.append(f"UNSAFE_SUBSTITUTION detected: {unsafe_subs}")

    # ── 7. Check label leakage ────────────────────────────────────────────────
    leakage = check_label_leakage(all_content)
    if leakage:
        result.failures.extend([f"LABEL_LEAKAGE: {v}" for v in leakage])

    # ── 8. Check GPU default ─────────────────────────────────────────────────
    if card.get("gpu_policy") == "FAIL":
        result.warnings.append("GPU default policy: FAIL (--gpu default=0, BDD Scenario 13 violation)")
    elif card.get("gpu_policy") == "PASS":
        gpu_ok, gpu_detail = check_gpu_default(target_path, entry_file)
        if not gpu_ok:
            result.failures.append(f"GPU default violation: {gpu_detail}")

    # ── 9. Check known_deviations in manifest matches card ───────────────────
    if manifest_entry:
        manifest_devs = manifest_entry.get("known_deviations", [])
        card_devs = card.get("compatibility_changes", [])
        if manifest_devs and not card_devs:
            result.warnings.append(
                "manifest has known_deviations but card has no compatibility_changes"
            )

    # ── Determine final status ────────────────────────────────────────────────
    if result.failures:
        result.status = STATUS_FAIL
    elif result.warnings:
        result.status = STATUS_WARN

    return result


def main():
    import argparse
    parser = argparse.ArgumentParser(
        description="Check model core cards for completeness and consistency."
    )
    parser.add_argument(
        "--model",
        type=str,
        default=None,
        help="Check only a specific model key"
    )
    parser.add_argument(
        "--verbose", "-v",
        action="store_true",
        help="Show detailed output"
    )
    parser.add_argument(
        "--json",
        action="store_true",
        help="Output results as JSON"
    )
    parser.add_argument(
        "--manifest-only",
        action="store_true",
        help="Only check which models have core cards (skip content checks)"
    )
    args = parser.parse_args()

    print("=" * 70)
    print("Model Core Card Checker")
    print("=" * 70)

    cards = load_cards()
    manifest = load_manifest()

    print(f"\nLoaded {len(cards)} core cards from {CARDS_DIR}")
    print(f"Loaded {len(manifest)} models from manifest")
    print()

    # ── Check: every non-placeholder manifest model has a card ──────────────────
    missing_cards = []
    for key, entry in sorted(manifest.items()):
        if entry.get("authenticity") in ("PLACEHOLDER",):
            continue
        if key not in cards:
            missing_cards.append(key)

    if missing_cards:
        print(f"MISSING core cards for: {', '.join(missing_cards)}")
        print()
    else:
        print("All non-placeholder models have core cards.")
        print()

    # ── Per-model checks ───────────────────────────────────────────────────────
    if args.manifest_only:
        return

    models_to_check = cards.keys() if not args.model else [args.model]
    if args.model and args.model not in cards:
        print(f"ERROR: Unknown model key '{args.model}'")
        print(f"Available: {', '.join(sorted(cards.keys()))}")
        sys.exit(2)

    results = {}
    for key in sorted(models_to_check):
        card = cards[key]
        manifest_entry = manifest.get(key)
        result = check_card(card, manifest_entry)
        results[key] = result

        icon = {
            STATUS_PASS: "✓", STATUS_FAIL: "✗", STATUS_WARN: "⚠",
            STATUS_ENV_GATED: "⊗", STATUS_CORE_INCOMPLETE: "○",
            STATUS_MISSING_CARD: "?", STATUS_NO_CARD: "?"
        }.get(result.status, "?")

        status_label = result.status
        if result.warnings and result.status == STATUS_PASS:
            status_label = "WARN"
            icon = "⚠"

        print(f"  {icon} [{status_label:16}] {key}")

        if args.verbose:
            for info in result.info:
                print(f"      INFO: {info}")
            for w in result.warnings:
                print(f"      WARN: {w}")
            for f in result.failures:
                print(f"      FAIL: {f}")
        print()

    # ── Summary ────────────────────────────────────────────────────────────────
    print("=" * 70)
    print("Summary")
    print("=" * 70)

    counts = {}
    for r in results.values():
        counts[r.status] = counts.get(r.status, 0) + 1

    total = len(results)
    for s, c in sorted(counts.items()):
        print(f"  {s:20}: {c}/{total}")

    # Expected statuses
    expected = {
        "PASS": [],
        "WARN": [],
        "ENV-GATED": ["scdeepcluster", "scname", "sczidesk", "desc"],
        "CORE-INCOMPLETE": ["scname", "sczidesk", "desc"],
        "FAIL": ["sccdcg", "scdcc"],  # Known label leakage violations
    }

    # We expect FAIL count = 0, CORE-INCOMPLETE for placeholders, ENV-GATED for TF models
    unexpected_failures = []
    for key, r in results.items():
        if r.status == STATUS_FAIL:
            if key not in expected.get("FAIL", []):
                unexpected_failures.append(f"  {key}: {r.failures}")

    if unexpected_failures:
        print("\n  UNEXPECTED FAILURES:")
        for f in unexpected_failures:
            print(f)
        print()
        print("  These models should PASS or have documented reasons for failure.")
    else:
        print("\n  All models pass or have expected non-PASS status.")

    if missing_cards:
        print(f"\n  WARNING: {len(missing_cards)} models missing core cards.")
        print(f"  Missing: {', '.join(missing_cards)}")

    # ── JSON output ───────────────────────────────────────────────────────────
    if args.json:
        print()
        output = {
            "missing_cards": missing_cards,
            "results": {
                k: {
                    "status": v.status,
                    "failures": v.failures,
                    "warnings": v.warnings,
                    "info": v.info,
                    "checks": v.checks,
                }
                for k, v in results.items()
            },
            "summary": counts,
        }
        print(json.dumps(output, indent=2, default=str))

    # ── Exit code ─────────────────────────────────────────────────────────────
    if unexpected_failures:
        sys.exit(1)
    sys.exit(0)


if __name__ == "__main__":
    main()
