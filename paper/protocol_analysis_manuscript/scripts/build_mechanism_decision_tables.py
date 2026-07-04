#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
MANUSCRIPT_ROOT = SCRIPT_DIR.parent
GENERATED_ROOT = MANUSCRIPT_ROOT / "generated"


def read_csv(path: Path) -> list[dict]:
    with open(path, "r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        raise ValueError(f"No rows to write for {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = []
    for row in rows:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def write_json(path: Path, data: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, sort_keys=True), encoding="utf-8")


def mean(values: list[float]) -> float:
    return sum(values) / len(values)


def latex_escape(value: object) -> str:
    return str(value).replace("_", r"\_")


def fmt(value: float) -> str:
    return f"{float(value):.6f}"


def phase13_decision(rows: list[dict]) -> tuple[list[dict], dict]:
    by_corruption: dict[str, list[dict]] = {}
    for row in rows:
        by_corruption.setdefault(row["corruption"], []).append(row)

    summaries = []
    for corruption, group in sorted(by_corruption.items()):
        summaries.append(
            {
                "study": "phase13_corruption_triad",
                "mechanism": corruption,
                "n_dataset_groups": len(group),
                "known_k_ari_mean": mean([float(row["kmeans_known_k.ari.mean"]) for row in group]),
                "fixed_leiden_ari_mean": mean([float(row["leiden_fixed.ari.mean"]) for row in group]),
                "effective_corruption_mean": mean([float(row["effective_corruption_rate.mean"]) for row in group]),
                "decision": "carry_forward" if corruption == "scmae_shuffle" else "diagnostic_only",
                "claim_boundary": "development evidence only; not validation",
            }
        )
    control = next(row for row in summaries if row["mechanism"] == "scmae_shuffle")
    best_known = max(summaries, key=lambda row: row["known_k_ari_mean"])
    best_leiden = max(summaries, key=lambda row: row["fixed_leiden_ari_mean"])
    evidence = {
        "recommended_corruption": "scmae_shuffle",
        "best_known_k_corruption": best_known["mechanism"],
        "best_fixed_leiden_corruption": best_leiden["mechanism"],
        "scmae_shuffle_known_k_mean": control["known_k_ari_mean"],
        "scmae_shuffle_fixed_leiden_mean": control["fixed_leiden_ari_mean"],
    }
    return summaries, evidence


def paired_delta(rows: list[dict], group_key: str, condition_key: str, baseline: str, candidate: str, metric: str) -> list[dict]:
    grouped: dict[str, dict[str, dict]] = {}
    for row in rows:
        grouped.setdefault(row[group_key], {})[row[condition_key]] = row
    deltas = []
    for group, values in sorted(grouped.items()):
        if baseline not in values or candidate not in values:
            raise ValueError(f"Missing {baseline}/{candidate} pair for {group}")
        deltas.append(
            {
                "group": group,
                "baseline": baseline,
                "candidate": candidate,
                "metric": metric,
                "baseline_value": float(values[baseline][metric]),
                "candidate_value": float(values[candidate][metric]),
                "delta": float(values[candidate][metric]) - float(values[baseline][metric]),
            }
        )
    return deltas


def phase14_decision(rows: list[dict]) -> tuple[list[dict], dict]:
    known = paired_delta(rows, "dataset", "variant", "control", "advmask", "kmeans_known_k.ari.mean")
    leiden = paired_delta(rows, "dataset", "variant", "control", "advmask", "leiden_fixed.ari.mean")
    known_mean = mean([row["delta"] for row in known])
    leiden_mean = mean([row["delta"] for row in leiden])
    positive_known = sum(1 for row in known if row["delta"] > 0)
    positive_leiden = sum(1 for row in leiden if row["delta"] > 0)
    grad_rows = [row for row in rows if row["variant"] == "advmask"]
    grad_positive = all(float(row["generator_grad_norm.mean"]) > 0 for row in grad_rows)
    decision = {
        "study": "phase14_advmask_triage",
        "mechanism": "advmask",
        "known_k_delta_mean": known_mean,
        "fixed_leiden_delta_mean": leiden_mean,
        "known_k_positive_datasets": positive_known,
        "fixed_leiden_positive_datasets": positive_leiden,
        "generator_grad_positive": grad_positive,
        "decision": "drop_or_downgrade",
        "claim_boundary": "active generator is not evidence of useful clustering representation",
    }
    return [decision], {"known_k_deltas": known, "fixed_leiden_deltas": leiden}


def attention_decision(rows: list[dict]) -> tuple[list[dict], dict]:
    known_axial = paired_delta(rows, "dataset", "role", "control", "axial", "kmeans_known_k.ari")
    known_param = paired_delta(rows, "dataset", "role", "control", "mlp_parammatched", "kmeans_known_k.ari")
    leiden_axial = paired_delta(rows, "dataset", "role", "control", "axial", "leiden_fixed.ari")
    decision = {
        "study": "attention_context_smoke",
        "mechanism": "current_axial_encoder",
        "known_k_delta_vs_control_mean": mean([row["delta"] for row in known_axial]),
        "fixed_leiden_delta_vs_control_mean": mean([row["delta"] for row in leiden_axial]),
        "known_k_positive_datasets": sum(1 for row in known_axial if row["delta"] > 0),
        "known_k_parammatched_delta_mean": mean([row["delta"] for row in known_param]),
        "decision": "do_not_use_as_rescue_path",
        "claim_boundary": "not evidence that all attention mechanisms fail",
    }
    return [decision], {"axial_known_k_deltas": known_axial, "axial_fixed_leiden_deltas": leiden_axial, "parammatched_known_k_deltas": known_param}


def resource_decision(rows: list[dict]) -> tuple[list[dict], dict]:
    by_condition: dict[str, list[dict]] = {}
    for row in rows:
        by_condition.setdefault(row["condition"], []).append(row)
    summaries = []
    for condition, group in sorted(by_condition.items()):
        summaries.append(
            {
                "study": "instrumented_resource_smoke",
                "mechanism": condition,
                "n_runs": len(group),
                "known_k_ari_mean": mean([float(row["kmeans_known_k.ari"]) for row in group]),
                "fixed_leiden_ari_mean": mean([float(row["leiden_fixed.ari"]) for row in group]),
                "wall_clock_seconds_mean": mean([float(row["wall_clock_seconds"]) for row in group]),
                "gpu_delta_mib_mean": mean([float(row["peak_gpu_memory_delta_mib"]) for row in group]),
                "decision": "diagnostic_control" if condition in {"control", "mlp_parammatched"} else "cost_without_stable_gain",
                "claim_boundary": "resource smoke only; not submission-scale runtime benchmark",
            }
        )
    known_advmask = paired_delta(rows, "dataset", "condition", "control", "advmask", "kmeans_known_k.ari")
    leiden_advmask = paired_delta(rows, "dataset", "condition", "control", "advmask", "leiden_fixed.ari")
    known_axial = paired_delta(rows, "dataset", "condition", "control", "axial", "kmeans_known_k.ari")
    return summaries, {
        "advmask_known_k_deltas": known_advmask,
        "advmask_fixed_leiden_deltas": leiden_advmask,
        "axial_known_k_deltas": known_axial,
    }


def write_latex_table(path: Path, rows: list[dict]) -> None:
    lines = [
        r"\begin{table}[ht]",
        r"\centering",
        r"\caption{Generated mechanism-decision table from current development evidence. Decisions are development-route decisions, not validation claims.}",
        r"\label{tab:generated-mechanism-decisions}",
        r"\scriptsize",
        r"\begin{tabular}{llll}",
        r"\toprule",
        r"Study & Mechanism & Decision & Claim boundary \\",
        r"\midrule",
    ]
    for row in rows:
        lines.append(
            " & ".join(
                [
                    latex_escape(row["study"]),
                    latex_escape(row["mechanism"]),
                    latex_escape(row["decision"]),
                    latex_escape(row["claim_boundary"]),
                ]
            )
            + r" \\"
        )
    lines.extend([r"\bottomrule", r"\end{tabular}", r"\end{table}", ""])
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")


def write_markdown(path: Path, rows: list[dict], evidence: dict) -> None:
    lines = [
        "# Mechanism Decision Matrix",
        "",
        "Status: generated development-evidence decision matrix. This is not validation evidence.",
        "",
        "## Decisions",
        "",
        "| study | mechanism | decision | claim boundary |",
        "|---|---|---|---|",
    ]
    for row in rows:
        lines.append(f"| {row['study']} | {row['mechanism']} | {row['decision']} | {row['claim_boundary']} |")
    lines.extend(
        [
            "",
            "## Key Paired Deltas",
            "",
            f"- Phase 14 AdvMask mean known-K ARI delta: {evidence['phase14']['decision'][0]['known_k_delta_mean']:.6f}.",
            f"- Phase 14 AdvMask mean fixed-Leiden ARI delta: {evidence['phase14']['decision'][0]['fixed_leiden_delta_mean']:.6f}.",
            f"- Attention smoke Axial mean known-K ARI delta vs control: {evidence['attention']['decision'][0]['known_k_delta_vs_control_mean']:.6f}.",
            f"- Resource smoke AdvMask mean known-K ARI delta vs control: {mean([row['delta'] for row in evidence['resource']['deltas']['advmask_known_k_deltas']]):.6f}.",
            f"- Resource smoke Axial mean known-K ARI delta vs control: {mean([row['delta'] for row in evidence['resource']['deltas']['axial_known_k_deltas']]):.6f}.",
            "",
            "## Safe Route",
            "",
            "Continue as protocol-analysis / diagnostic paper unless a new mechanism later passes a fresh development gate. Do not claim AdvMask, current Axial, full CAAM synergy, publication validation, or unknown-K clustering superiority from current evidence.",
            "",
        ]
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")


def build_outputs(args: argparse.Namespace) -> None:
    phase13_rows = read_csv(args.phase13_summary)
    phase14_rows = read_csv(args.phase14_summary)
    attention_rows = read_csv(args.attention_runs)
    resource_rows = read_csv(args.resource_runs)

    phase13_rows_out, phase13_evidence = phase13_decision(phase13_rows)
    phase14_rows_out, phase14_evidence = phase14_decision(phase14_rows)
    attention_rows_out, attention_evidence = attention_decision(attention_rows)
    resource_rows_out, resource_evidence = resource_decision(resource_rows)

    decision_rows = phase13_rows_out + phase14_rows_out + attention_rows_out + resource_rows_out
    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    write_csv(output_dir / "mechanism_decisions.csv", decision_rows)
    write_latex_table(output_dir / "mechanism_decisions.tex", decision_rows)

    evidence = {
        "claim_scope": "development evidence only; no validation or sealed test data",
        "phase13": phase13_evidence,
        "phase14": {"decision": phase14_rows_out, "deltas": phase14_evidence},
        "attention": {"decision": attention_rows_out, "deltas": attention_evidence},
        "resource": {"decision": resource_rows_out, "deltas": resource_evidence},
        "safe_route": "protocol_analysis",
        "unsafe_claims": [
            "AdvMask improves clustering",
            "current Axial improves clustering",
            "Axial plus AdvMask synergy",
            "development evidence validates a publication-ready method",
            "known-K ARI proves unknown-K clustering quality",
        ],
    }
    write_json(output_dir / "mechanism_decisions.json", evidence)
    write_markdown(output_dir / "mechanism_decisions.md", decision_rows, evidence)


def main() -> int:
    parser = argparse.ArgumentParser(description="Build CAAM mechanism-decision tables from generated development evidence.")
    parser.add_argument("--phase13-summary", type=Path, default=GENERATED_ROOT / "data/phase13_corruption_summary.csv")
    parser.add_argument("--phase14-summary", type=Path, default=GENERATED_ROOT / "data/phase14_advmask_summary.csv")
    parser.add_argument("--attention-runs", type=Path, default=GENERATED_ROOT / "data/attention_context_smoke_runs.csv")
    parser.add_argument("--resource-runs", type=Path, default=GENERATED_ROOT / "instrumented_resource_smoke/instrumented_resource_smoke.csv")
    parser.add_argument("--output-dir", type=Path, default=GENERATED_ROOT / "mechanism_decisions")
    args = parser.parse_args()
    build_outputs(args)
    print(f"Wrote mechanism decision tables to {args.output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
