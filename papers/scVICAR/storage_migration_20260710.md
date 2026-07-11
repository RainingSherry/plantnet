# Local storage migration (2026-07-10)

Large run products and regenerable staging files were moved from the `/home`
workspace to the matching paths under
`/data/luolie/biopipeline/dimension-reduction/plantnet`.  The original paths are
symbolic links, so existing scripts keep the same interfaces.

Before each source directory was replaced, an incremental `rsync -a --partial`
copy was completed and a dry-run comparison reported no source differences.
The links were then checked as readable directories before the local backups
were removed.

| Workspace link | Data target |
|---|---|
| `experiment_reports` | `/data/luolie/biopipeline/dimension-reduction/plantnet/experiment_reports` |
| `benchmark_outputs` | `/data/luolie/biopipeline/dimension-reduction/plantnet/benchmark_outputs` |
| `experimental_retired_models/scMAEs/runs` | `/data/luolie/biopipeline/dimension-reduction/plantnet/experimental_retired_models/scMAEs/runs` |
| `result/scmae_retired_methods_20260706_full` | `/data/luolie/biopipeline/dimension-reduction/plantnet/result/scmae_retired_methods_20260706_full` |
| `result/scmae_nm_rg_contrast_param_search_20260707_full` | `/data/luolie/biopipeline/dimension-reduction/plantnet/result/scmae_nm_rg_contrast_param_search_20260707_full` |
| `result/scmae_unified_comparison_benchmark_20260710` | `/data/luolie/biopipeline/dimension-reduction/plantnet/result/scmae_unified_comparison_benchmark_20260710` |
| `papers/scVICAR/.staging` | `/data/luolie/biopipeline/dimension-reduction/plantnet/papers/scVICAR/.staging` |

The migration increased free space on the `/home` filesystem from approximately
52 GiB to 86 GiB.  Paper source, aggregate tables, manifests, and final figures
remain in the workspace.  The pre-existing disconnected path
`result/scmae_all_methods_20260705_full` was not modified.
