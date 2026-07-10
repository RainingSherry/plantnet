from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from .config import DATASETS, PAPER_ROOT, SEEDS, VARIANTS
from .run_baseline import BASELINES


DISPLAY = {
    "nomix": "NoMix",
    "random_mix": "RandomMix",
    "fixed": "scVICAR-F",
    "topology_edge_only": "T-Edge",
    "topology_gate_only": "T-Gate",
    "topology_full": "scVICAR-T",
}


def require_primary(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path)
    formal = frame[frame["execution_mode"] == "formal"].copy()
    expected = len(DATASETS) * len(VARIANTS) * len(SEEDS)
    keys = formal[["dataset", "variant", "seed"]].drop_duplicates()
    if len(formal) != expected or len(keys) != expected:
        raise RuntimeError(f"Refusing manuscript assets from partial primary results: {len(formal)}/{expected}")
    if formal[["ari", "nmi", "acc", "f1_macro"]].isna().any().any():
        raise ValueError("Primary result matrix contains missing manuscript metrics")
    return formal


def write_confirmatory(formal: pd.DataFrame, contrasts: pd.DataFrame, manuscript: Path, tables: Path) -> None:
    required = {"F_vs_NoMix", "T_vs_NoMix", "T_vs_F"}
    ari = contrasts[contrasts["metric"] == "ari"].set_index("contrast")
    if not required.issubset(ari.index):
        raise ValueError("ARI contrast table is incomplete")
    seed_mean = formal.groupby(["dataset", "variant"], as_index=False)[["ari", "nmi", "f1_macro"]].mean()
    overall = seed_mean.groupby("variant")[["ari", "nmi", "f1_macro"]].agg(["mean", "std"]).reindex(VARIANTS)
    lines = [
        r"\begin{table}[t]",
        r"\centering",
        r"\caption{Matched-backbone confirmatory performance. Values are means across dataset-level seed averages; dispersion is the standard deviation across datasets.}",
        r"\label{tab:confirmatory_overall}",
        r"\begin{tabular}{lccc}",
        r"\toprule",
        r"Method & ARI & NMI & Macro-F1 \\",
        r"\midrule",
    ]
    for variant, row in overall.iterrows():
        lines.append(
            f"{DISPLAY[variant]} & {row[('ari', 'mean')]:.3f} ({row[('ari', 'std')]:.3f}) & "
            f"{row[('nmi', 'mean')]:.3f} ({row[('nmi', 'std')]:.3f}) & "
            f"{row[('f1_macro', 'mean')]:.3f} ({row[('f1_macro', 'std')]:.3f}) \\\\"
        )
    lines.extend([r"\bottomrule", r"\end{tabular}", r"\end{table}"])
    tables.mkdir(parents=True, exist_ok=True)
    (tables / "confirmatory_overall.tex").write_text("\n".join(lines) + "\n", encoding="utf-8")

    def sentence(key: str, left: str, right: str) -> str:
        row = ari.loc[key]
        return (
            f"{left} versus {right} yielded a mean paired ARI difference of {row['mean_delta']:.3f} "
            f"(dataset bootstrap 95\\% CI {row['ci95_low']:.3f} to {row['ci95_high']:.3f}; "
            f"wins/ties/losses {int(row['wins'])}/{int(row['ties'])}/{int(row['losses'])}; "
            f"Holm-adjusted permutation p={row['holm_p']:.4f})."
        )

    fragment = "\n\n".join([
        r"\input{../tables/generated/confirmatory_overall}",
        sentence("F_vs_NoMix", "scVICAR-F", "NoMix"),
        sentence("T_vs_NoMix", "scVICAR-T", "NoMix"),
        sentence("T_vs_F", "scVICAR-T", "scVICAR-F")
        + " The latter is a clean-graph comparison and is not interpreted as a requirement that topology adaptation dominate on every dataset.",
    ])
    (manuscript / "generated/confirmatory_results.tex").write_text(fragment + "\n", encoding="utf-8")


def write_downstream(dataset_metrics: pd.DataFrame, manuscript: Path, tables: Path) -> None:
    expected = len(DATASETS) * len(VARIANTS) * 3
    if len(dataset_metrics) != expected:
        raise RuntimeError(f"Refusing partial downstream manuscript assets: {len(dataset_metrics)}/{expected}")
    selected = [
        ("marker", None, "recovery_recovery_at_100", "Recovery@100"),
        ("marker", None, "annotation_f1_macro", "Marker annotation macro-F1"),
        ("linear_probe", 0.1, "probe_f1_macro", r"10\% probe macro-F1"),
        ("linear_probe", 0.3, "probe_f1_macro", r"30\% probe macro-F1"),
    ]
    rows = []
    for task, fraction, metric, label in selected:
        subset = dataset_metrics[dataset_metrics["task"] == task]
        subset = subset[subset["label_fraction"].isna()] if fraction is None else subset[np.isclose(subset["label_fraction"], fraction)]
        rows.append((label, subset.groupby("variant")[metric].mean()))
    lines = [
        r"\begin{table}[t]",
        r"\centering",
        r"\caption{Downstream utility after averaging split seeds and then model seeds within each dataset.}",
        r"\label{tab:downstream_overall}",
        r"\begin{tabular}{lccc}",
        r"\toprule",
        r"Endpoint & NoMix & scVICAR-F & scVICAR-T \\",
        r"\midrule",
    ]
    for label, means in rows:
        lines.append(f"{label} & {means['nomix']:.3f} & {means['fixed']:.3f} & {means['topology_full']:.3f} \\\\")
    lines.extend([r"\bottomrule", r"\end{tabular}", r"\end{table}"])
    (tables / "downstream_overall.tex").write_text("\n".join(lines) + "\n", encoding="utf-8")
    fragment = (
        r"\input{../tables/generated/downstream_overall}" + "\n"
        + "Table~\\ref{tab:downstream_overall} reports all prespecified downstream endpoints after split-level and model-seed aggregation. "
        + "Paired effects, confidence intervals, and Holm-adjusted tests are provided in the generated contrast CSV; no best seed or split is selected.\n"
    )
    (manuscript / "generated/downstream_results.tex").write_text(fragment, encoding="utf-8")


def write_leiden(overall: pd.DataFrame, contrasts: pd.DataFrame, manuscript: Path, tables: Path) -> None:
    if len(overall) != len(VARIANTS) or set(overall["variant"]) != set(VARIANTS):
        raise RuntimeError("Refusing incomplete fixed-Leiden overall results")
    if set(overall["n_datasets"]) != {len(DATASETS)}:
        raise RuntimeError("Fixed-Leiden overall results do not use six dataset units")
    confirmatory = contrasts[contrasts["metric"] == "ari"].set_index("contrast")
    required = {"F_vs_NoMix", "T_vs_NoMix", "T_vs_F"}
    if len(contrasts) != 12 or not required.issubset(confirmatory.index):
        raise RuntimeError("Fixed-Leiden contrast family is incomplete")
    indexed = overall.set_index("variant").reindex(VARIANTS)
    lines = [
        r"\begin{table}[t]", r"\centering",
        r"\caption{Fixed-resolution Leiden sensitivity analysis. Values are means (SD) across the six dataset-level seed averages; resolution is 1.0 and is never label-tuned.}",
        r"\label{tab:leiden_overall}", r"\begin{tabular}{lccc}", r"\toprule",
        r"Method & ARI & NMI & Macro-F1 \\ ", r"\midrule",
    ]
    for variant, row in indexed.iterrows():
        lines.append(
            f"{DISPLAY[variant]} & {row['ari_mean']:.3f} ({row['ari_sd']:.3f}) & "
            f"{row['nmi_mean']:.3f} ({row['nmi_sd']:.3f}) & "
            f"{row['f1_macro_mean']:.3f} ({row['f1_macro_sd']:.3f}) \\\\"
        )
    lines.extend([r"\bottomrule", r"\end{tabular}", r"\end{table}"])
    (tables / "leiden_overall.tex").write_text("\n".join(lines) + "\n", encoding="utf-8")

    def effect(key: str, left: str, right: str) -> str:
        row = confirmatory.loc[key]
        return (
            f"{left} versus {right}: mean paired ARI difference {row['mean_delta']:.3f} "
            f"(hierarchical-bootstrap 95\\% CI {row['ci95_low']:.3f} to {row['ci95_high']:.3f}; "
            f"wins/ties/losses {int(row['wins'])}/{int(row['ties'])}/{int(row['losses'])}; "
            f"Holm-adjusted permutation p={row['holm_p']:.4f})."
        )

    fragment = "\n\n".join([
        r"\input{../tables/generated/leiden_overall}",
        "Fixed-resolution Leiden changed absolute scores and some method ordering but did not support a universal adaptive advantage. "
        + effect("F_vs_NoMix", "scVICAR-F", "NoMix"),
        effect("T_vs_NoMix", "scVICAR-T", "NoMix"),
        effect("T_vs_F", "scVICAR-T", "scVICAR-F"),
    ])
    (manuscript / "generated/leiden_results.tex").write_text(fragment + "\n", encoding="utf-8")


def write_stress(stress_runs: pd.DataFrame, manuscript: Path, tables: Path) -> None:
    if len(stress_runs) != 126:
        raise RuntimeError(f"Refusing partial stress manuscript assets: {len(stress_runs)}/126")
    current = stress_runs[stress_runs["estimator"] == "current"].copy()
    expected = 3 * 2 * 3 * 5
    if len(current[["dataset", "variant", "seed", "contamination"]].drop_duplicates()) != expected:
        raise RuntimeError("Stress contamination grid is incomplete")
    seed_mean = current.groupby(
        ["dataset", "variant", "contamination"], as_index=False
    ).mean(numeric_only=True)
    overall = seed_mean.groupby(["variant", "contamination"], as_index=False).mean(numeric_only=True)
    pivot = overall.pivot(index="contamination", columns="variant", values="ari")
    clean = pivot.loc[0.0]
    degradation = pivot.subtract(clean, axis=1)
    lines = [
        r"\begin{table}[t]", r"\centering",
        r"\caption{Graph-contamination stress results. ARI values average model seeds within dataset and then the three stress-test datasets. Degradation is relative to each variant's clean graph.}",
        r"\label{tab:stress_overall}", r"\begin{tabular}{rrrrr}", r"\toprule",
        r"Injected edges (\%) & F ARI & T ARI & F degradation & T degradation \\ ", r"\midrule",
    ]
    for contamination in sorted(pivot.index):
        lines.append(
            f"{int(round(100 * contamination))} & {pivot.loc[contamination, 'fixed']:.3f} & "
            f"{pivot.loc[contamination, 'topology_full']:.3f} & "
            f"{degradation.loc[contamination, 'fixed']:.3f} & "
            f"{degradation.loc[contamination, 'topology_full']:.3f} \\\\"
        )
    lines.extend([r"\bottomrule", r"\end{tabular}", r"\end{table}"])
    (tables / "stress_overall.tex").write_text("\n".join(lines) + "\n", encoding="utf-8")

    final = max(pivot.index)
    difference = degradation.loc[final, "topology_full"] - degradation.loc[final, "fixed"]
    diagnostic = seed_mean[(seed_mean["variant"] == "topology_full") & (seed_mean["contamination"] > 0)]
    auroc = float(diagnostic["affinity_same_edge_auroc"].mean())
    gate = float(diagnostic["gate_purity_spearman"].mean())
    estimator = stress_runs[stress_runs["contamination"] == 0].groupby(
        ["dataset", "variant", "estimator"], as_index=False
    )["ari"].mean()
    estimator_means = estimator.groupby("estimator")["ari"].mean()
    fragment = "\n\n".join([
        r"\input{../tables/generated/stress_overall}",
        (
            f"At 100\\% injected cross-class edges, mean ARI degradation was "
            f"{degradation.loc[final, 'fixed']:.3f} for scVICAR-F and "
            f"{degradation.loc[final, 'topology_full']:.3f} for scVICAR-T; "
            f"the descriptive T--F degradation difference was {difference:.3f} "
            f"(positive values indicate slower degradation for T)."
        ),
        (
            f"Across nonzero contamination levels, the topology-informed affinity had mean same-edge AUROC {auroc:.3f}, "
            f"and the mean gate--purity Spearman association was {gate:.3f}. "
            "These are mechanism diagnostics, not independent biological replicates."
        ),
        (
            "At zero injected contamination, descriptive ARI means for the implemented, uniform-sample, and full-neighborhood estimators were "
            f"{estimator_means.get('current', np.nan):.3f}, {estimator_means.get('uniform_sample', np.nan):.3f}, and "
            f"{estimator_means.get('full', np.nan):.3f}, respectively."
        ),
    ])
    (manuscript / "generated/stress_results.tex").write_text(fragment + "\n", encoding="utf-8")


def write_baselines(baseline: pd.DataFrame, formal: pd.DataFrame, manuscript: Path, tables: Path) -> None:
    expected = len(DATASETS) * len(BASELINES) * len(SEEDS)
    if len(baseline) != expected or len(baseline[["dataset", "method", "seed"]].drop_duplicates()) != expected:
        raise RuntimeError(f"Refusing partial baseline manuscript assets: {len(baseline)}/{expected}")
    external_dataset = baseline.groupby(["dataset", "method"], as_index=False)[
        ["ari", "nmi", "f1_macro", "runtime_seconds"]
    ].mean()
    matched_dataset = formal[formal["variant"].isin(["nomix", "fixed", "topology_full"])].groupby(
        ["dataset", "variant"], as_index=False
    )[["ari", "nmi", "f1_macro", "runtime_seconds"]].mean().rename(columns={"variant": "method"})
    combined = pd.concat([external_dataset, matched_dataset], ignore_index=True)
    order = [*BASELINES, "nomix", "fixed", "topology_full"]
    display = {
        "pca_kmeans": "PCA+KMeans", "scmae": "Original scMAE", "scvi": "scVI",
        "scdcc": "scDCC", "scdeepcluster": "scDeepCluster", "scrcl": "scRCL",
        "nomix": "Matched NoMix", "fixed": "scVICAR-F", "topology_full": "scVICAR-T",
    }
    overall = combined.groupby("method")[["ari", "nmi", "f1_macro", "runtime_seconds"]].agg(["mean", "std"]).reindex(order)
    lines = [
        r"\begin{table*}[t]", r"\centering",
        r"\caption{External baselines and principal matched variants. Values are means (SD) across six dataset-level seed averages. scDCC, scDeepCluster, and scRCL use known $K$ during their published clustering-oriented training; other rows use $K$ only for post-hoc KMeans.}",
        r"\label{tab:external_baselines}", r"\begin{tabular}{lrrrr}", r"\toprule",
        r"Method & ARI & NMI & Macro-F1 & Runtime (s) \\ ", r"\midrule",
    ]
    for method, row in overall.iterrows():
        lines.append(
            f"{display[method]} & {row[('ari', 'mean')]:.3f} ({row[('ari', 'std')]:.3f}) & "
            f"{row[('nmi', 'mean')]:.3f} ({row[('nmi', 'std')]:.3f}) & "
            f"{row[('f1_macro', 'mean')]:.3f} ({row[('f1_macro', 'std')]:.3f}) & "
            f"{row[('runtime_seconds', 'mean')]:.1f} ({row[('runtime_seconds', 'std')]:.1f}) \\\\"
        )
    lines.extend([r"\bottomrule", r"\end{tabular}", r"\end{table*}"])
    (tables / "external_baselines.tex").write_text("\n".join(lines) + "\n", encoding="utf-8")
    (manuscript / "generated/baseline_results.tex").write_text(
        r"\input{../tables/generated/external_baselines}" + "\n"
        + "External methods are reported under their frozen published adapters and are not used to attribute backbone-matched gains. "
        + "Every row averages all three seeds within dataset before the six-dataset summary; no historical value or best seed is substituted. "
        + "Wall-clock times are descriptive because CPU/GPU frameworks and interpreters differ across published implementations.\n",
        encoding="utf-8",
    )


def validate_optional_complete(baseline: Path, stress: Path, require_all: bool) -> None:
    if baseline.is_file():
        frame = pd.read_csv(baseline)
        expected = len(DATASETS) * len(BASELINES) * len(SEEDS)
        keys = frame[["dataset", "method", "seed"]].drop_duplicates()
        if require_all and (len(frame) != expected or len(keys) != expected):
            raise RuntimeError(f"External baseline matrix is partial: {len(frame)}/{expected}")
    elif require_all:
        raise FileNotFoundError(baseline)
    if stress.is_file():
        frame = pd.read_csv(stress)
        if require_all and len(frame) != 126:
            raise RuntimeError(f"Stress matrix is partial: {len(frame)}/126")
    elif require_all:
        raise FileNotFoundError(stress)


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate manuscript prose/tables only from complete frozen aggregates")
    parser.add_argument("--run-master", type=Path, default=PAPER_ROOT / "experiments/protocol_v1/run_master.csv")
    parser.add_argument("--contrasts", type=Path, default=PAPER_ROOT / "tables/confirmatory_contrasts_ari.csv")
    parser.add_argument("--baseline-master", type=Path, default=PAPER_ROOT / "experiments/baselines_v1/run_master.csv")
    parser.add_argument("--stress-runs", type=Path, default=PAPER_ROOT / "experiments/stress_v1/stress_runs.csv")
    parser.add_argument("--downstream-metrics", type=Path, default=PAPER_ROOT / "experiments/downstream_v1/dataset_variant_metrics.csv")
    parser.add_argument(
        "--leiden-overall", type=Path,
        default=PAPER_ROOT / "experiments/leiden_fixed_v1/aggregate/variant_overall_metrics.csv",
    )
    parser.add_argument(
        "--leiden-contrasts", type=Path,
        default=PAPER_ROOT / "experiments/leiden_fixed_v1/aggregate/contrasts.csv",
    )
    parser.add_argument("--require-all", action="store_true")
    args = parser.parse_args()
    formal = require_primary(args.run_master)
    manuscript = PAPER_ROOT / "manuscript"
    tables = PAPER_ROOT / "tables/generated"
    (manuscript / "results_macros.tex").write_text(
        "% Generated from the complete frozen primary aggregate.\n"
        "\\newcommand{\\FormalDatasetCount}{6}\n"
        "\\newcommand{\\FormalSeedCount}{3}\n"
        "\\newcommand{\\FormalRunCount}{108}\n"
        "\\newcommand{\\FormalResultStatus}{primary confirmatory matrix complete; robustness and downstream validation in progress}\n"
        "\\newcommand{\\DevelopmentDatasetCount}{16}\n",
        encoding="utf-8",
    )
    write_confirmatory(formal, pd.read_csv(args.contrasts), manuscript, tables)
    if args.leiden_overall.is_file() and args.leiden_contrasts.is_file():
        write_leiden(
            pd.read_csv(args.leiden_overall), pd.read_csv(args.leiden_contrasts), manuscript, tables
        )
    elif args.require_all:
        raise FileNotFoundError("Complete fixed-Leiden aggregate is required")
    validate_optional_complete(args.baseline_master, args.stress_runs, args.require_all)
    if args.baseline_master.is_file():
        baseline = pd.read_csv(args.baseline_master)
        expected_baselines = len(DATASETS) * len(BASELINES) * len(SEEDS)
        if len(baseline) == expected_baselines:
            write_baselines(baseline, formal, manuscript, tables)
        elif args.require_all:
            raise RuntimeError(f"External baseline matrix is partial: {len(baseline)}/{expected_baselines}")
    if args.stress_runs.is_file():
        write_stress(pd.read_csv(args.stress_runs), manuscript, tables)
    elif args.require_all:
        raise FileNotFoundError(args.stress_runs)
    if args.downstream_metrics.is_file():
        write_downstream(pd.read_csv(args.downstream_metrics), manuscript, tables)
    elif args.require_all:
        raise FileNotFoundError(args.downstream_metrics)
    print("Generated manuscript assets from complete frozen aggregates")


if __name__ == "__main__":
    main()
