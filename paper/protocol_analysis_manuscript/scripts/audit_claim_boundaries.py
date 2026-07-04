#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import re
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
MANUSCRIPT_ROOT = SCRIPT_DIR.parent
GENERATED_ROOT = MANUSCRIPT_ROOT / "generated"
MAIN_TEX = MANUSCRIPT_ROOT / "main.tex"
README = MANUSCRIPT_ROOT / "README.md"
CHECKLIST = MANUSCRIPT_ROOT / "REPRODUCIBILITY_CHECKLIST.md"


FORBIDDEN_PATTERNS = {
    "positive_caam_method_claim": re.compile(r"\bCAAM[- ]scMAE\s+(improves|outperforms|achieves|is\s+validated)", re.IGNORECASE),
    "advmask_positive_claim": re.compile(r"\bAdvMask\s+(improves|outperforms|boosts|is\s+better|is\s+validated)", re.IGNORECASE),
    "axial_positive_claim": re.compile(r"\bAxial\b[^.\n]{0,80}\b(improves|outperforms|boosts|rescues|is\s+validated)", re.IGNORECASE),
    "synergy_claim": re.compile(r"\bsynergy\b|\bsynergistic\b", re.IGNORECASE),
    "publication_validation_claim": re.compile(r"\b(publication[- ]level|validated|validation)\s+(evidence|claim|method|result)", re.IGNORECASE),
    "unknown_k_overclaim": re.compile(r"\bfully\s+unsupervised\b|\bunknown[- ]K\b", re.IGNORECASE),
}

NEGATION_TERMS = (
    "not",
    "no ",
    "does not",
    "do not",
    "cannot",
    "unsupported",
    "not support",
    "not supported",
    "not evidence",
    "unsafe",
    "future",
    "should not",
    "must not",
    "without",
    "lacks",
    "negative",
)

REQUIRED_PHRASES = (
    "development evidence only",
    "Validation has not been run",
    "Known-\\(K\\) K-means ARI remains a useful development diagnostic",
    "does not support a positive CAAM-scMAE method paper",
    "not a new positive claim",
)


def load_json(path: Path) -> dict:
    with open(path, "r", encoding="utf-8") as handle:
        return json.load(handle)


def read_text(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def write_json(path: Path, data: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, sort_keys=True), encoding="utf-8")


def generated_file_count() -> int:
    return sum(1 for path in GENERATED_ROOT.rglob("*") if path.is_file())


def csv_rows(path: Path) -> list[dict]:
    with open(path, "r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def mean(values: list[float]) -> float:
    return sum(values) / len(values)


def condition_summary(rows: list[dict]) -> dict[str, dict[str, float]]:
    out: dict[str, dict[str, float]] = {}
    for condition in sorted({row["condition"] for row in rows}):
        group = [row for row in rows if row["condition"] == condition]
        out[condition] = {
            "n_runs": float(len(group)),
            "known_k_ari_mean": mean([float(row["kmeans_known_k.ari"]) for row in group]),
            "fixed_leiden_ari_mean": mean([float(row["leiden_fixed.ari"]) for row in group]),
            "wall_clock_seconds_mean": mean([float(row["wall_clock_seconds"]) for row in group]),
            "gpu_delta_mib_mean": mean([float(row["peak_gpu_memory_delta_mib"]) for row in group]),
        }
    return out


def is_negated(text: str, start: int, end: int) -> bool:
    window = text[max(0, start - 140) : min(len(text), end + 140)].lower()
    return any(term in window for term in NEGATION_TERMS)


def scan_forbidden_claims(text: str) -> list[dict]:
    findings = []
    for name, pattern in FORBIDDEN_PATTERNS.items():
        for match in pattern.finditer(text):
            if is_negated(text, match.start(), match.end()):
                continue
            line = text.count("\n", 0, match.start()) + 1
            excerpt = text[max(0, match.start() - 80) : min(len(text), match.end() + 80)].replace("\n", " ")
            findings.append({"name": name, "line": line, "match": match.group(0), "excerpt": excerpt})
    return findings


def find_missing_required_phrases(text: str) -> list[str]:
    return [phrase for phrase in REQUIRED_PHRASES if phrase not in text]


def extract_count_mentions(text: str, unit: str) -> list[int]:
    pattern = re.compile(rf"\b(\d+)\s+{re.escape(unit)}\b", re.IGNORECASE)
    return [int(match.group(1)) for match in pattern.finditer(text)]


def consistency_findings(texts: dict[str, str], evidence: dict) -> list[dict]:
    findings = []
    expected_files = int(evidence["generated_file_count"])
    expected_resource_runs = int(evidence["instrumented_resource_smoke"]["n_runs"])
    for name, text in texts.items():
        for count in extract_count_mentions(text, "files"):
            if count != expected_files:
                findings.append(
                    {
                        "name": "generated_file_count_mismatch",
                        "file": name,
                        "observed_text_count": count,
                        "expected_count": expected_files,
                    }
                )
        if "instrumented resource" in text or "resource smoke" in text:
            for count in extract_count_mentions(text, "instrumented resource-smoke runs"):
                if count != expected_resource_runs:
                    findings.append(
                        {
                            "name": "instrumented_resource_run_count_mismatch",
                            "file": name,
                            "observed_text_count": count,
                            "expected_count": expected_resource_runs,
                        }
                    )
        stale_phrases = (
            "single-dataset resource smoke",
            "one development dataset only",
            "one development dataset",
            "single-dataset smoke",
        )
        for phrase in stale_phrases:
            if phrase in text:
                findings.append({"name": "stale_resource_scope_phrase", "file": name, "phrase": phrase})
    return findings


def build_evidence() -> dict:
    protocol_manifest = load_json(GENERATED_ROOT / "artifact_manifest.json")
    resource_manifest = load_json(GENERATED_ROOT / "instrumented_resource_smoke/artifact_manifest.json")
    resource_rows = csv_rows(GENERATED_ROOT / "instrumented_resource_smoke/instrumented_resource_smoke.csv")
    mechanism_rows = csv_rows(GENERATED_ROOT / "mechanism_decisions/mechanism_decisions.csv")
    feature_space = load_json(GENERATED_ROOT / "feature_space_smoke/feature_space_smoke.json")
    mask_ratio = load_json(GENERATED_ROOT / "mask_ratio_smoke/mask_ratio_smoke.json")
    summary = condition_summary(resource_rows)
    return {
        "generated_file_count": generated_file_count(),
        "phase13_runs": int(protocol_manifest["phase13_runs"]),
        "phase14_runs": int(protocol_manifest["phase14_runs"]),
        "attention_runs": int(protocol_manifest["attention_runs"]),
        "instrumented_resource_smoke": {
            "n_runs": len(resource_rows),
            "datasets": [item["name"] for item in resource_manifest["datasets"]],
            "seed": int(resource_manifest["seed"]),
            "epochs": int(resource_manifest["epochs"]),
            "conditions": resource_manifest["conditions"],
            "gpu_memory_sampling_mode": resource_manifest["gpu_memory_sampling_mode"],
            "condition_summary": summary,
        },
        "mechanism_decisions": [
            {
                "study": row["study"],
                "mechanism": row["mechanism"],
                "decision": row["decision"],
                "claim_boundary": row["claim_boundary"],
            }
            for row in mechanism_rows
        ],
        "feature_space_smoke": feature_space,
        "mask_ratio_smoke": mask_ratio,
    }


def write_markdown(path: Path, audit: dict) -> None:
    resource = audit["evidence"]["instrumented_resource_smoke"]
    summary = resource["condition_summary"]
    mechanism_decisions = audit["evidence"]["mechanism_decisions"]
    feature_space = audit["evidence"]["feature_space_smoke"]
    feature_rows = {row["role"]: row for row in feature_space["rows"]}
    feature_deltas = feature_space["deltas"]
    mask_ratio = audit["evidence"]["mask_ratio_smoke"]
    mask_rows = mask_ratio["aggregate_rows"]
    mask_summary = mask_ratio["summary"]
    lines = [
        "# Claim Boundary Audit",
        "",
        "Status: generated manuscript-support audit. This is not validation evidence.",
        "",
        "## Evidence Snapshot",
        "",
        f"- Generated artifact files: {audit['evidence']['generated_file_count']}",
        f"- Phase 13 runs: {audit['evidence']['phase13_runs']}",
        f"- Phase 14 runs: {audit['evidence']['phase14_runs']}",
        f"- Attention/context smoke runs: {audit['evidence']['attention_runs']}",
        f"- Instrumented resource-smoke runs: {resource['n_runs']}",
        f"- Feature-space smoke runs: {feature_space['n_runs']}",
        f"- Mask-ratio smoke runs: {mask_ratio['n_runs']}",
        f"- Resource-smoke datasets: {', '.join(resource['datasets'])}",
        f"- Resource-smoke seed/epochs: {resource['seed']} / {resource['epochs']}",
        f"- GPU memory mode: {resource['gpu_memory_sampling_mode']}",
        "",
        "## Instrumented Resource Summary",
        "",
        "| condition | known-K ARI mean | fixed-Leiden ARI mean | wall sec mean | GPU delta MiB mean |",
        "|---|---:|---:|---:|---:|",
    ]
    for condition in ("control", "advmask", "axial", "mlp_parammatched"):
        row = summary[condition]
        lines.append(
            f"| {condition} | {row['known_k_ari_mean']:.6f} | {row['fixed_leiden_ari_mean']:.6f} | "
            f"{row['wall_clock_seconds_mean']:.2f} | {row['gpu_delta_mib_mean']:.1f} |"
        )

    lines.extend(
        [
            "",
            "## Generated Mechanism Decisions",
            "",
            "| study | mechanism | decision | claim boundary |",
            "|---|---|---|---|",
        ]
    )
    for row in mechanism_decisions:
        lines.append(f"| {row['study']} | {row['mechanism']} | {row['decision']} | {row['claim_boundary']} |")

    lines.extend(
        [
            "",
            "## Feature-Space Smoke Boundary",
            "",
            "| role | genes | student params | known-K ARI | fixed-Leiden ARI |",
            "|---|---:|---:|---:|---:|",
        ]
    )
    for role in ("hvg2000", "full_gene_stress"):
        row = feature_rows[role]
        lines.append(
            f"| {role} | {row['n_genes']} | {row['student_trainable_params']} | "
            f"{row['kmeans_known_k_ari']:.6f} | {row['leiden_fixed_ari']:.6f} |"
        )
    lines.extend(
        [
            f"- Full-gene minus HVG known-K ARI: {feature_deltas['full_minus_hvg_known_k_ari']:.6f}.",
            f"- Full-gene minus HVG fixed-Leiden ARI: {feature_deltas['full_minus_hvg_leiden_fixed_ari']:.6f}.",
            f"- Full-gene/HVG student parameter ratio: {feature_deltas['full_over_hvg_parameter_ratio']:.2f}x.",
            "- This is a single-dataset dense-MLP development smoke, not validation and not a universal claim against full-gene modeling.",
            "",
            "## Mask-Ratio Smoke Boundary",
            "",
            "| mask ratio | runs | masked effective mean | global effective mean | known-K ARI mean | known-K ARI std | fixed-Leiden ARI mean | fixed-Leiden ARI std |",
            "|---:|---:|---:|---:|---:|---:|---:|---:|",
        ]
    )
    for row in mask_rows:
        lines.append(
            f"| {row['mask_ratio']:.1f} | {row['n_runs']} | {row['effective_corruption_rate_masked_mean']:.6f} | "
            f"{row['global_effective_change_rate_estimate_mean']:.6f} | {row['kmeans_known_k_ari_mean']:.6f} | "
            f"{row['kmeans_known_k_ari_std']:.6f} | {row['leiden_fixed_ari_mean']:.6f} | {row['leiden_fixed_ari_std']:.6f} |"
        )
    lines.extend(
        [
            f"- Best mean known-K ARI in this smoke: mask ratio {mask_summary['best_known_k_mask_ratio']:.1f}.",
            f"- Best mean fixed-Leiden ARI in this smoke: mask ratio {mask_summary['best_leiden_mask_ratio']:.1f}.",
            f"- Known-K ARI mean range: {mask_summary['known_k_mean_range']:.6f}; fixed-Leiden ARI mean range: {mask_summary['leiden_mean_range']:.6f}.",
            f"- Maximum known-K seed standard deviation: {mask_summary['max_known_k_seed_std']:.6f}; maximum fixed-Leiden seed standard deviation: {mask_summary['max_leiden_seed_std']:.6f}.",
            "- This is a single-dataset mask-ratio sensitivity diagnostic, not a tuning decision for validation.",
            "",
            "## Supported Claims",
            "",
            "- Development evidence supports a protocol-analysis route rather than a positive CAAM method route.",
            "- HVG 2000 remains the current dense-MLP feature-space default under development evidence.",
            "- Nominal mask ratio can change global effective perturbation without producing monotonic downstream clustering behavior.",
            "- Mask-ratio figures are development diagnostics and must be read with the table-level claim boundary.",
            "- scMAE-style shuffle remains the safest corruption choice to carry forward in the frozen development protocol.",
            "- Effective corruption diagnostics and downstream clustering quality can diverge.",
            "- AdvMask is active but unsupported as a main clustering-improvement mechanism under current development evidence.",
            "- The current Axial encoder is unsupported as a rescue path and should not be treated as evidence against all attention mechanisms.",
            "- Parameter-matched MLP controls are necessary before interpreting architecture effects.",
            "",
            "## Unsupported Claims",
            "",
            "- CAAM-scMAE improves clustering.",
            "- AdvMask improves clustering.",
            "- The current Axial encoder improves clustering.",
            "- Axial plus AdvMask has synergy.",
            "- Development evidence is publication-level validation.",
            "- Known-K ARI proves fully unsupervised unknown-K clustering quality.",
            "- HVG 2000 is validated by the feature-space smoke.",
            "- Full-gene or gene-token feature-space models are generally inferior.",
            "- A different mask ratio should be selected from the single-dataset smoke.",
            "- Mask ratio 0.6 is generally best because it has the highest fixed-Leiden mean in this smoke.",
            "- The frozen validation mask ratio should be changed based on this smoke.",
            "",
            "## Audit Findings",
            "",
        ]
    )
    if audit["status"] == "pass":
        lines.append("No blocking claim-boundary findings.")
    else:
        for finding in audit["findings"]:
            lines.append(f"- {finding}")
    lines.extend(
        [
            "",
            "## Safe Next-Step Boundary",
            "",
            "Validation is still not run. Any validation pass must use the frozen protocol and must not tune corruption, mask policy, architecture, loss, clustering resolution, or manuscript claims.",
            "",
        ]
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")


def build_audit() -> dict:
    texts = {
        "main.tex": read_text(MAIN_TEX),
        "README.md": read_text(README),
        "REPRODUCIBILITY_CHECKLIST.md": read_text(CHECKLIST),
    }
    evidence = build_evidence()
    findings = []
    forbidden = scan_forbidden_claims(texts["main.tex"])
    if forbidden:
        findings.extend({"type": "forbidden_claim", **item} for item in forbidden)
    missing = find_missing_required_phrases(texts["main.tex"])
    if missing:
        findings.extend({"type": "missing_required_phrase", "phrase": phrase} for phrase in missing)
    findings.extend({"type": "consistency", **item} for item in consistency_findings(texts, evidence))
    return {
        "status": "pass" if not findings else "fail",
        "findings": findings,
        "evidence": evidence,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Audit manuscript claim boundaries against generated CAAM evidence.")
    parser.add_argument("--output-dir", type=Path, default=GENERATED_ROOT / "claim_audit")
    args = parser.parse_args()
    audit = build_audit()
    write_json(args.output_dir / "claim_boundary_audit.json", audit)
    write_markdown(args.output_dir / "claim_boundary_audit.md", audit)
    print(f"claim_audit_status={audit['status']}")
    print(f"findings={len(audit['findings'])}")
    return 0 if audit["status"] == "pass" else 1


if __name__ == "__main__":
    raise SystemExit(main())
