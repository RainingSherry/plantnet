# APA-scMAE v2 ablation on processed_scmae

## Run summary

- Branch: `apa-v2-objective`
- Data root: `/data/luolie/biopipeline/scCluBench/data/processed_scmae`
- Datasets: `Limb_Muscle`, `Quake_10x_Spleen`, `Guo`, `Macosko`, `Tosches`, `Young`
- Seeds: `42`, `2024`, `3407`
- Variants: A random Student, B VICReg, C EMA Teacher, D prototype consistency, E full v2 with Generator
- Runs: `90/90` complete
- Failures: `0`
- `metrics_by_run.csv`: `180` rows, covering `kmeans_known_k` and `leiden_fixed`
- `metrics_summary.csv`: `60` groups

## Main conclusion

The full v2 Generator variant is not consistently dominant across evaluation modes. It improves over the no-generator prototype variant on most `kmeans_known_k` paired comparisons, but this does not hold reliably for `leiden_fixed`, where several datasets show no stable E-over-D gain.

The safest paper conclusion is:

> The representation-driven Student objectives are the main effective component; the current Generator is not yet the consistently dominant contributor.

## Key paired E vs D findings

- `kmeans_known_k`: E beats D by paired ARI on `Limb_Muscle`, `Quake_10x_Spleen`, `Guo`, `Tosches`, and `Young`, but not on `Macosko`.
- `leiden_fixed`: E beats D clearly on `Limb_Muscle`; results are mixed or worse on `Quake_10x_Spleen`, `Guo`, `Macosko`, and `Tosches`; `Young` has a small ARI gain but not a paired NMI gain.
- If reporting one balanced interpretation, emphasize that Generator benefits are dataset- and clustering-method-dependent.

## Files

- `tables/apa_v2_ablation_metrics_by_run.csv`: per-run metrics, always two rows per run.
- `tables/apa_v2_ablation_metrics_summary.csv`: mean/std summary by dataset, variant, and cluster method.
- `tables/apa_v2_ablation_paper_table.csv`: combined ARI/NMI wide paper table.
- `tables/apa_v2_ablation_paper_table_ari.csv`: ARI-only paper table.
- `tables/apa_v2_ablation_paper_table_nmi.csv`: NMI-only paper table.
- `tables/apa_v2_ablation_E_vs_D_paired.csv`: paired seed comparison for E vs D.
- `metadata/run_environment.json`: branch, commit, Python, data root, datasets, variants, and start time.
- `metadata/task_manifest.json`: complete task manifest with commands and resolved paths.
- `metadata/failures.json`: failure list; expected to be `[]`.
- `code/apa_v2_ablation_processed_scmae_paper.py`: driver script used for this run.

## Notes

- The report package intentionally excludes run directories, embeddings, checkpoints, h5ad files, and other large artifacts.
- Raw benchmark artifacts remain under `benchmark_outputs/APA_scMAE_v2_ablation_processed_scmae_paper/` and should not be committed unless explicitly requested.
