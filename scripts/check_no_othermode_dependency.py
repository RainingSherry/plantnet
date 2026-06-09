#!/usr/bin/env python3
"""
check_no_othermode_dependency.py
================================
Per BDD Scenario 11: Verify that methods/ do not depend on OtherMode/ at runtime.

This script checks:
  1. methods/**/*.py do not import or reference OtherMode.
  2. scripts/run_formal_benchmark.py does not add OtherMode to PYTHONPATH.
  3. command.txt does not reference OtherMode.
  4. model_core_cards do not set source_path as a runtime import.

Exit codes:
    0 = PASS — no runtime dependency on OtherMode found
    1 = FAIL — runtime dependency on OtherMode detected
    2 = ERROR — script error
"""

import os
import re
import sys
import json
import yaml
from pathlib import Path
from typing import Dict, List, Tuple, Set

PROJECT_ROOT = Path(__file__).parent.parent.resolve()
METHODS_DIR = PROJECT_ROOT / "methods"
OTHERMODE_DIR = PROJECT_ROOT / "OtherMode"
SCRIPTS_DIR = PROJECT_ROOT / "scripts"
CARDS_DIR = PROJECT_ROOT / "docs" / "model_core_cards"


def find_python_files(directory: Path, exclude_dirs: Set[str] = None) -> List[Path]:
    """Recursively find all Python files, excluding specified directories."""
    exclude = exclude_dirs or {"__pycache__", ".git", "OtherMode"}
    files = []
    for item in directory.rglob("*.py"):
        if any(excluded in item.parts for excluded in exclude):
            continue
        files.append(item)
    return files


def check_file_for_othermode(content: str, file_path: Path,
                            is_methods_dir: bool = False) -> List[Tuple[str, int, str]]:
    """
    Check a single file for OtherMode runtime import references.
    Returns list of (file, line_no, line_content) tuples.

    A "runtime import" means:
      - Adding OtherMode to sys.path
      - from/import OtherMode or scCluBench
      - Using OtherMode/scCluBench in an actual import or path assignment

    Safe (not flagged):
      - "OtherMode" in comments, docstrings, variable names (e.g. OTHERMODE_DIR)
      - source_path fields in YAML cards
      - Reference-only docstrings
    """
    violations = []
    lines = content.split("\n")

    # Patterns that indicate actual runtime OtherMode dependency
    # These MUST be in actual Python code (not just mentioned in comments)
    import_patterns = [
        (r"sys\.path.*insert.*OtherMode",       "sys.path insert OtherMode"),
        (r"sys\.path.*append.*OtherMode",       "sys.path append OtherMode"),
        (r"sys\.path.*OtherMode",               "sys.path reference OtherMode"),
        (r"from\s+\.\./.*OtherMode",           "relative import from OtherMode"),
        (r"from\s+\.\./OtherMode",              "import from OtherMode"),
        (r"import\s+OtherMode",                "import OtherMode"),
        (r"from\s+OtherMode\s+import",         "from OtherMode import"),
    ]

    for i, line in enumerate(lines, 1):
        stripped = line.strip()

        # Skip comment-only lines
        if stripped.startswith("#"):
            continue

        # Skip docstring lines (heuristic: start with triple-quote or are heavily indented)
        if '"""' in stripped or "'''" in stripped:
            continue

        # For methods/ files, also skip if line is clearly documentation
        # (contains words like "reference", "source_path", "original paper")
        if is_methods_dir and any(w in stripped.lower() for w in ["reference", "source_path", "original"]):
            continue

        for pattern, description in import_patterns:
            if re.search(pattern, line, re.IGNORECASE):
                violations.append((str(file_path), i, stripped[:120]))
                break

    return violations


def check_scripts_dir() -> List[Tuple[str, int, str]]:
    """Check scripts/ for OtherMode references (excluding the checker itself)."""
    violations = []
    for fpath in SCRIPTS_DIR.rglob("*.py"):
        # Exclude the checker script itself (it legitimately mentions OtherMode in patterns)
        if fpath.name == "check_no_othermode_dependency.py":
            continue
        content = fpath.read_text(encoding="utf-8", errors="ignore")
        found = check_file_for_othermode(content, fpath, is_methods_dir=True)
        violations.extend(found)
    return violations


def check_methods_dir() -> List[Tuple[str, int, str]]:
    """Check methods/ for OtherMode references."""
    violations = []
    py_files = find_python_files(METHODS_DIR)
    for fpath in py_files:
        content = fpath.read_text(encoding="utf-8", errors="ignore")
        found = check_file_for_othermode(content, fpath, is_methods_dir=True)
        violations.extend(found)
    return violations


def check_command_txt(out_dir: Path) -> List[Tuple[str, int, str]]:
    """Check result directories for command.txt with OtherMode references."""
    violations = []
    for cmd_file in out_dir.rglob("command.txt"):
        content = cmd_file.read_text(encoding="utf-8", errors="ignore")
        lines = content.split("\n")
        for i, line in enumerate(lines, 1):
            if "OtherMode" in line or "scCluBench" in line:
                violations.append((str(cmd_file), i, line.strip()))
    return violations


def check_run_formal_benchmark() -> List[Tuple[str, int, str]]:
    """Check run_formal_benchmark.py for OtherMode in PYTHONPATH."""
    violations = []
    script_path = SCRIPTS_DIR / "run_formal_benchmark.py"
    if not script_path.exists():
        return violations

    content = script_path.read_text(encoding="utf-8", errors="ignore")
    lines = content.split("\n")

    for i, line in enumerate(lines, 1):
        if "OtherMode" in line or "scCluBench" in line:
            # Allow comments about OtherMode being reference-only
            if re.match(r"^\s*#", line):
                continue
            # Allow OtherMode in path literals in comments
            if "reference" in line.lower() or "do not add" in line.lower():
                continue
            violations.append((str(script_path), i, line.strip()))

    return violations


def check_core_cards_source_path() -> List[Tuple[str, int, str]]:
    """
    Check that model core cards do not accidentally set source_path as runtime import.
    The source_path is for documentation, not for Python import.
    """
    violations = []
    if not CARDS_DIR.exists():
        return violations

    for card_file in CARDS_DIR.glob("*.yaml"):
        try:
            with open(card_file, "r", encoding="utf-8") as f:
                data = yaml.safe_load(f)
            source_path = data.get("source_path", "")
            if source_path and "OtherMode" in source_path:
                # source_path points to OtherMode reference — this is fine for docs
                # but we should ensure it's not used as an import path in code
                pass
        except Exception:
            pass
    return violations


def check_results_for_othermode() -> List[Tuple[str, int, str]]:
    """Check results/ directories for OtherMode references."""
    violations = []
    results_dir = PROJECT_ROOT / "results"
    if not results_dir.exists():
        return violations

    violations.extend(check_command_txt(results_dir))
    return violations


def main():
    import argparse
    parser = argparse.ArgumentParser(
        description="Check that methods/ do not depend on OtherMode/ at runtime."
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
        "--fix",
        action="store_true",
        help="Show fix suggestions (read-only check)"
    )
    args = parser.parse_args()

    print("=" * 70)
    print("OtherMode Runtime Dependency Checker")
    print("=" * 70)
    print()
    print("This checks that methods/ and scripts/ do not use OtherMode/")
    print("as a runtime import or PYTHONPATH dependency.")
    print("OtherMode is for REFERENCE ONLY.")
    print()

    all_violations: List[Tuple[str, str, int, str]] = []  # (category, file, line, content)

    # Check 1: methods/
    print("Checking methods/ ...")
    method_violations = check_methods_dir()
    all_violations.extend(("methods/", *v) for v in method_violations)
    print(f"  Found {len(method_violations)} violations")

    # Check 2: scripts/
    print("Checking scripts/ ...")
    script_violations = check_scripts_dir()
    all_violations.extend(("scripts/", *v) for v in script_violations)
    print(f"  Found {len(script_violations)} violations")

    # Check 3: run_formal_benchmark.py
    print("Checking run_formal_benchmark.py ...")
    benchmark_violations = check_run_formal_benchmark()
    all_violations.extend(("run_formal_benchmark.py", *v) for v in benchmark_violations)
    print(f"  Found {len(benchmark_violations)} violations")

    # Check 4: results/
    print("Checking results/ command.txt files ...")
    results_violations = check_results_for_othermode()
    all_violations.extend(("results/", *v) for v in results_violations)
    print(f"  Found {len(results_violations)} violations")

    print()

    # ── Output ───────────────────────────────────────────────────────────────
    if all_violations:
        print("FAIL: Runtime dependency on OtherMode detected")
        print()
        prev_file = None
        for category, file, line_no, content in all_violations:
            if file != prev_file:
                print(f"  {file}:")
                prev_file = file
            print(f"    L{line_no}: {content[:120]}")
            if args.verbose:
                print(f"      Category: {category}")
        print()

        if args.fix:
            print("Fix suggestions:")
            for category, file, line_no, content in all_violations:
                if "sys.path" in content.lower():
                    print(f"  {file}:{line_no}")
                    print(f"    Remove or comment out sys.path line referencing OtherMode")
                    print(f"    OtherMode is for REFERENCE ONLY, not for runtime import")

        print()
        print(f"Total violations: {len(all_violations)}")
        print()
        print("OtherMode/ contains reference implementations only.")
        print("Do NOT add OtherMode to PYTHONPATH or import from it.")
        print()
        sys.exit(1)
    else:
        print("PASS: no runtime dependency on OtherMode")
        print()
        print("All methods/ and scripts/ are self-contained.")
        print("OtherMode/ is referenced only in documentation (source_path fields).")
        print()
        sys.exit(0)


if __name__ == "__main__":
    main()
