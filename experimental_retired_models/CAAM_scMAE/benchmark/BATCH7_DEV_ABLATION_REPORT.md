# CAAM-scMAE Batch 7 Development Ablation Report

## Scope

- Preflight only used label-free expression, batch code, library size, and zero ratio.
- No model mechanism, loss, donor logic, mask logic, or formal benchmark registration was changed.
- No development ablation was launched for datasets blocked by preflight.
- Do not claim `synergy_confirmed` from this batch.

## Preflight Results

| Dataset | Status | Reason | n_cells | n_genes | Estimated eligibility rate | Estimated budget deficit rate | Param-match status |
| --- | --- | --- | ---: | ---: | ---: | ---: | --- |
| `data/processed/Quake_Smart-seq2_Lung.h5ad` | `blocked_by_budget_deficit` | estimated budget deficit exceeds 1% threshold | 1676 | 23341 | 0.1737 | 1.0000 | pass, hidden_dim=26, gap=0.000008 |
| `data/其他/Mouse_Pancreas_1.h5ad` | `blocked_by_budget_deficit` | estimated budget deficit exceeds 1% threshold | 1886 | 14878 | 0.1702 | 1.0000 | pass, hidden_dim=41, gap=0.000009 |
| `data/processed_scmae/Limb_Muscle.h5ad` | `blocked_by_budget_deficit` | estimated budget deficit exceeds 1% threshold | 3909 | 23341 | 0.0962 | 1.0000 | pass, hidden_dim=26, gap=0.000008 |

## Ablation Results

No development dataset passed preflight, so no Batch 7 quick internal ablation was run.

| Dataset | claim_status | delta_AB | full_minus_axial | full_minus_advmask | full_minus_parammatched_mlp |
| --- | --- | --- | --- | --- | --- |
| `data/processed/Quake_Smart-seq2_Lung.h5ad` | not_run_preflight_blocked | N/A | N/A | N/A | N/A |
| `data/其他/Mouse_Pancreas_1.h5ad` | not_run_preflight_blocked | N/A | N/A | N/A | N/A |
| `data/processed_scmae/Limb_Muscle.h5ad` | not_run_preflight_blocked | N/A | N/A | N/A | N/A |

## Output Locations

- Preflight report: `results/CAAM_scMAE_preflight_batch7/preflight_report.json`
- Preflight summary: `results/CAAM_scMAE_preflight_batch7/preflight_summary.csv`

These result files are local artifacts and are not committed.
