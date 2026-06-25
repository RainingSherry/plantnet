#!/usr/bin/env python3
from __future__ import annotations

import argparse
import copy
import csv
import json
import sys
from pathlib import Path
from typing import Any

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[4]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from methods.DeepLearning.CAAM_scMAE.data.donor_candidates import DonorCandidateProvider
from methods.DeepLearning.CAAM_scMAE.data.preprocessing import load_caam_data
from methods.DeepLearning.CAAM_scMAE.registry import DEFAULT_CONFIG


STATUS_PASS = "pass"
STATUS_BUDGET = "blocked_by_budget_deficit"
STATUS_PARAM = "blocked_by_parameter_match"
STATUS_RUNTIME = "blocked_by_runtime"
STATUS_OTHER = "blocked_by_other"


def _linear_params(in_dim: int, out_dim: int) -> int:
    return int(in_dim) * int(out_dim) + int(out_dim)


def _layernorm_params(dim: int) -> int:
    return 2 * int(dim)


def _mha_params(embed_dim: int) -> int:
    d = int(embed_dim)
    return (3 * d * d) + (3 * d) + (d * d) + d


def _shared_student_params(n_genes: int, latent_dim: int, conditioning: str) -> int:
    g = int(n_genes)
    z = int(latent_dim)
    mask_head = _linear_params(z, g)
    decoder_in = z if conditioning == "none" else z + g
    decoder_hidden = max(z, g)
    decoder = _linear_params(decoder_in, decoder_hidden) + _linear_params(decoder_hidden, g)
    return mask_head + decoder


def _mlp_encoder_params(n_genes: int, latent_dim: int, hidden_dim: int) -> int:
    g = int(n_genes)
    z = int(latent_dim)
    h = int(hidden_dim)
    return (
        _linear_params(g, h)
        + _layernorm_params(h)
        + _linear_params(h, z)
        + _layernorm_params(z)
        + _linear_params(z, z)
        + _layernorm_params(z)
    )


def _gene_axis_block_params(token_dim: int) -> int:
    d = int(token_dim)
    ffn = _linear_params(d, 4 * d) + _linear_params(4 * d, d)
    return _layernorm_params(d) + _mha_params(d) + _layernorm_params(d) + ffn


def _cell_axis_params(token_dim: int) -> int:
    d = int(token_dim)
    projections = 4 * _linear_params(d, d)
    ffn = _linear_params(d, 4 * d) + _linear_params(4 * d, d)
    return projections + _layernorm_params(d) + _layernorm_params(d) + ffn


def _axial_encoder_params(config: dict[str, Any]) -> int:
    d = int(config["axial"]["token_dim"])
    z = int(config["model"]["latent_dim"])
    m = int(config["axial"]["n_gene_modules"])
    layers = int(config["axial"]["gene_attention_layers"])
    tokenizer = _linear_params(1, d) + _linear_params(d, d) + (m * d)
    gene_axis = layers * _gene_axis_block_params(d)
    cell_axis = _cell_axis_params(d)
    projection = _linear_params(d, z) + _layernorm_params(z)
    return tokenizer + gene_axis + cell_axis + projection


def estimate_student_params(n_genes: int, config: dict[str, Any], encoder_type: str, mlp_hidden_dim: int | None = None) -> int:
    z = int(config["model"]["latent_dim"])
    shared = _shared_student_params(n_genes, z, str(config["model"]["decoder_mask_conditioning"]))
    if encoder_type == "axial":
        encoder = _axial_encoder_params(config)
    elif encoder_type == "mlp":
        hidden = int(mlp_hidden_dim if mlp_hidden_dim is not None else config["model"]["mlp_hidden_dim"])
        encoder = _mlp_encoder_params(n_genes, z, hidden)
    else:
        raise ValueError(f"Unknown encoder_type={encoder_type!r}")
    return int(shared + encoder)


def infer_param_matched_hidden_dim_estimate(n_genes: int, config: dict[str, Any]) -> dict[str, Any]:
    target = estimate_student_params(n_genes, config, "axial")

    def count_for(hidden: int) -> int:
        return estimate_student_params(n_genes, config, "mlp", hidden)

    low = 1
    high = max(4, int(config["model"]["mlp_hidden_dim"]))
    while count_for(high) < target and high < 65536:
        low = high
        high *= 2
    for _ in range(24):
        mid = (low + high) // 2
        if mid <= low:
            break
        if count_for(mid) < target:
            low = mid
        else:
            high = mid

    best_hidden = high
    best_count = count_for(high)
    best_gap = abs(best_count - target) / max(1, target)
    for hidden in range(max(1, low - 64), high + 65):
        count = count_for(hidden)
        gap = abs(count - target) / max(1, target)
        if gap < best_gap:
            best_hidden = hidden
            best_count = count
            best_gap = gap
    return {
        "axial_student_params": int(target),
        "mlp_parammatched_hidden_dim": int(best_hidden),
        "mlp_parammatched_student_params": int(best_count),
        "relative_param_gap": float(best_gap),
        "can_parameter_match": bool(best_gap <= 0.05),
    }


def estimate_eligibility(
    *,
    x: np.ndarray,
    candidates: np.ndarray,
    mask_ratio: float,
    atol: float,
    rtol: float,
    seed: int,
    max_cells: int,
    gene_chunk_size: int,
) -> dict[str, Any]:
    rng = np.random.default_rng(seed)
    n_cells, n_genes = x.shape
    sample_n = min(int(max_cells), int(n_cells))
    if sample_n <= 0:
        raise ValueError("No cells available for eligibility estimation.")
    sampled_cells = rng.choice(np.arange(n_cells, dtype=np.int64), size=sample_n, replace=False)
    eligible_counts = np.zeros(sample_n, dtype=np.int64)
    total_positions = sample_n * n_genes
    candidate_pool_size = int(candidates.shape[1])
    for start in range(0, n_genes, int(gene_chunk_size)):
        end = min(n_genes, start + int(gene_chunk_size))
        width = end - start
        slots = rng.integers(0, candidate_pool_size, size=(sample_n, width), endpoint=False)
        donor_indices = np.take_along_axis(candidates[sampled_cells], slots, axis=1)
        genes = np.arange(start, end, dtype=np.int64)[None, :]
        replacement = x[donor_indices, genes]
        original = x[sampled_cells, start:end]
        eligible = ~np.isclose(replacement, original, atol=atol, rtol=rtol)
        eligible_counts += eligible.sum(axis=1).astype(np.int64)
    budget = int(round(float(mask_ratio) * n_genes))
    deficit = np.maximum(0, budget - eligible_counts)
    return {
        "estimated_cells": int(sample_n),
        "estimated_gene_count": int(n_genes),
        "estimated_eligibility_rate": float(eligible_counts.sum() / max(1, total_positions)),
        "estimated_budget_per_cell": int(budget),
        "estimated_budget_deficit_rate": float(np.mean(deficit > 0)),
        "estimated_budget_deficit_mean": float(deficit.mean()),
        "estimated_min_eligible_per_cell": int(eligible_counts.min()),
        "estimated_median_eligible_per_cell": float(np.median(eligible_counts)),
    }


def _fallback_rates(levels: dict[str, int], n_cells: int) -> dict[str, float]:
    denom = max(1, int(n_cells))
    return {str(k): float(v) / denom for k, v in levels.items()}


def preflight_one(
    data_path: str | Path,
    *,
    config: dict[str, Any] | None = None,
    max_runtime_genes: int = 5000,
    max_runtime_params: int = 250_000_000,
    max_estimate_cells: int = 512,
    gene_chunk_size: int = 2048,
) -> dict[str, Any]:
    cfg = copy.deepcopy(config or DEFAULT_CONFIG)
    cfg["benchmark_mode"] = True
    cfg["preprocessing"]["input_mode"] = "log1p"
    cfg["preprocessing"]["n_top_genes"] = 2000
    cfg["preprocessing"]["scale_input"] = False
    path = Path(data_path)
    try:
        bundle = load_caam_data(
            str(path),
            input_mode="log1p",
            target_sum=float(cfg["preprocessing"]["target_sum"]),
            n_top_genes=int(cfg["preprocessing"]["n_top_genes"]),
            scale_input=False,
            benchmark_mode=True,
            seed=int(cfg.get("seed", 0)),
        )
        x = np.asarray(bundle.x, dtype=np.float32)
        n_cells, n_genes = x.shape
        donor = DonorCandidateProvider(
            x,
            bundle.batch_code,
            bundle.library_size,
            bundle.zero_ratio,
            candidate_pool_size=int(cfg["corruption"]["candidate_pool_size"]),
            library_size_bins=int(cfg["corruption"]["library_size_bins"]),
            zero_ratio_bins=int(cfg["corruption"]["zero_ratio_bins"]),
            atol=float(cfg["mask"]["changed_tolerance_abs"]),
            rtol=float(cfg["mask"]["changed_tolerance_rel"]),
            seed=int(cfg.get("seed", 0)),
        )
        eligibility = estimate_eligibility(
            x=x,
            candidates=donor.candidates,
            mask_ratio=float(cfg["mask"]["ratio"]),
            atol=float(cfg["mask"]["changed_tolerance_abs"]),
            rtol=float(cfg["mask"]["changed_tolerance_rel"]),
            seed=int(cfg.get("seed", 0)),
            max_cells=max_estimate_cells,
            gene_chunk_size=gene_chunk_size,
        )
        params = infer_param_matched_hidden_dim_estimate(n_genes, cfg)
        fallback_levels = donor.stats.get("fallback_levels", {})
        max_budget_deficit = float(cfg["corruption"]["max_budget_deficit_fraction"])
        strict_effective_budget = bool(cfg["corruption"].get("strict_effective_budget", False))
        runtime_blocked = bool(n_genes > int(max_runtime_genes) or params["axial_student_params"] > int(max_runtime_params))
        if strict_effective_budget and eligibility["estimated_budget_deficit_rate"] > max_budget_deficit:
            status = STATUS_BUDGET
            reason = (
                f"estimated_budget_deficit_rate={eligibility['estimated_budget_deficit_rate']:.4f} "
                f"exceeds max_budget_deficit_fraction={max_budget_deficit:.4f}"
            )
        elif not params["can_parameter_match"]:
            status = STATUS_PARAM
            reason = f"relative_param_gap={params['relative_param_gap']:.4f} exceeds 0.0500"
        elif runtime_blocked:
            status = STATUS_RUNTIME
            reason = (
                f"n_genes={n_genes} or axial_student_params={params['axial_student_params']} "
                f"exceeds quick thresholds n_genes<={max_runtime_genes}, params<={max_runtime_params}"
            )
        else:
            status = STATUS_PASS
            reason = "quick ablation preflight passed"
        report = {
            "data_path": str(path),
            "dataset_name": path.stem,
            "quick_ablation_status": status,
            "status_reason": reason,
            "n_cells": int(n_cells),
            "n_genes": int(n_genes),
            "sparsity": float((x <= 0.0).mean()),
            "library_size": {
                "min": float(np.min(bundle.library_size)),
                "mean": float(np.mean(bundle.library_size)),
                "median": float(np.median(bundle.library_size)),
                "max": float(np.max(bundle.library_size)),
            },
            "zero_ratio": {
                "min": float(np.min(bundle.zero_ratio)),
                "mean": float(np.mean(bundle.zero_ratio)),
                "median": float(np.median(bundle.zero_ratio)),
                "max": float(np.max(bundle.zero_ratio)),
            },
            "donor_fallback_levels": {str(k): int(v) for k, v in fallback_levels.items()},
            "donor_fallback_rates": _fallback_rates(fallback_levels, n_cells),
            **eligibility,
            **params,
            "runtime_thresholds": {
                "max_runtime_genes": int(max_runtime_genes),
                "max_runtime_params": int(max_runtime_params),
            },
            "preprocessing": {
                "benchmark_mode": True,
                "input_mode": "log1p",
                "n_top_genes": int(cfg["preprocessing"]["n_top_genes"]),
                "scale_input": False,
                "strict_effective_budget": strict_effective_budget,
            },
        }
        return report
    except Exception as exc:
        return {
            "data_path": str(path),
            "dataset_name": path.stem,
            "quick_ablation_status": STATUS_OTHER,
            "status_reason": f"{type(exc).__name__}: {exc}",
        }


SUMMARY_COLUMNS = [
    "dataset_name",
    "data_path",
    "quick_ablation_status",
    "status_reason",
    "n_cells",
    "n_genes",
    "sparsity",
    "library_size_mean",
    "zero_ratio_mean",
    "donor_fallback_matched",
    "donor_fallback_batch",
    "donor_fallback_global",
    "estimated_eligibility_rate",
    "estimated_budget_deficit_rate",
    "axial_student_params",
    "mlp_parammatched_hidden_dim",
    "relative_param_gap",
    "can_parameter_match",
]


def _summary_row(report: dict[str, Any]) -> dict[str, Any]:
    levels = report.get("donor_fallback_levels", {})
    return {
        "dataset_name": report.get("dataset_name"),
        "data_path": report.get("data_path"),
        "quick_ablation_status": report.get("quick_ablation_status"),
        "status_reason": report.get("status_reason"),
        "n_cells": report.get("n_cells"),
        "n_genes": report.get("n_genes"),
        "sparsity": report.get("sparsity"),
        "library_size_mean": report.get("library_size", {}).get("mean") if isinstance(report.get("library_size"), dict) else None,
        "zero_ratio_mean": report.get("zero_ratio", {}).get("mean") if isinstance(report.get("zero_ratio"), dict) else None,
        "donor_fallback_matched": levels.get("matched"),
        "donor_fallback_batch": levels.get("batch"),
        "donor_fallback_global": levels.get("global"),
        "estimated_eligibility_rate": report.get("estimated_eligibility_rate"),
        "estimated_budget_deficit_rate": report.get("estimated_budget_deficit_rate"),
        "axial_student_params": report.get("axial_student_params"),
        "mlp_parammatched_hidden_dim": report.get("mlp_parammatched_hidden_dim"),
        "relative_param_gap": report.get("relative_param_gap"),
        "can_parameter_match": report.get("can_parameter_match"),
    }


def write_outputs(reports: list[dict[str, Any]], output_dir: str | Path) -> tuple[Path, Path]:
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    report_path = out / "preflight_report.json"
    summary_path = out / "preflight_summary.csv"
    report_path.write_text(json.dumps({"datasets": reports}, indent=2), encoding="utf-8")
    with summary_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=SUMMARY_COLUMNS)
        writer.writeheader()
        for report in reports:
            writer.writerow(_summary_row(report))
    return report_path, summary_path


def run_preflight(paths: list[str], output_dir: str | Path, **kwargs) -> tuple[list[dict[str, Any]], Path, Path]:
    reports = [preflight_one(path, **kwargs) for path in paths]
    report_path, summary_path = write_outputs(reports, output_dir)
    return reports, report_path, summary_path


def main() -> int:
    parser = argparse.ArgumentParser(description="Preflight CAAM-scMAE datasets for quick internal ablation.")
    parser.add_argument("data_paths", nargs="+", help="One or more h5ad dataset paths.")
    parser.add_argument("--output_dir", default="results/CAAM_scMAE_preflight")
    parser.add_argument("--max_runtime_genes", type=int, default=5000)
    parser.add_argument("--max_runtime_params", type=int, default=250_000_000)
    parser.add_argument("--max_estimate_cells", type=int, default=512)
    parser.add_argument("--gene_chunk_size", type=int, default=2048)
    args = parser.parse_args()
    reports, report_path, summary_path = run_preflight(
        args.data_paths,
        args.output_dir,
        max_runtime_genes=args.max_runtime_genes,
        max_runtime_params=args.max_runtime_params,
        max_estimate_cells=args.max_estimate_cells,
        gene_chunk_size=args.gene_chunk_size,
    )
    print(f"Wrote {report_path}")
    print(f"Wrote {summary_path}")
    for report in reports:
        print(f"{report['dataset_name']}: {report['quick_ablation_status']} - {report['status_reason']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
