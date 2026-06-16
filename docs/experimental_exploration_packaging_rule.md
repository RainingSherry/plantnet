# Experimental Exploration Packaging Rule

Effective date: 2026-06-16

This project often contains exploratory method code, ablations, and benchmark outputs that are useful for future research but should not be mixed into production method paths without context. From now on, every exploratory experiment that is worth uploading to GitHub should be archived as a self-contained package.

## Required Rule

For each exploratory experiment, create:

```text
experiment_packages/<experiment_id>/
```

The package must include:

```text
README.md
code/
results/summaries/
results/run_artifacts/
notes/
MANIFEST.files
MANIFEST.sha256
```

Use a stable experiment id with date:

```text
<topic>_<YYYYMMDD>
```

Example:

```text
experiment_packages/neighbormix_stochastic_regularization_20260616/
```

## What To Include

Code:

```text
code/
```

Include snapshots of all scripts required to rerun or summarize the experiment. Prefer copying only the touched exploratory files rather than the whole repository.

Results:

```text
results/summaries/
```

Include final tables and interpretation files, at minimum:

```text
main_results.csv
group_summary.csv
neighbor_diagnostics.csv
interpretation.md
```

Per-run evidence:

```text
results/run_artifacts/<dataset>/<method>/seed<seed>/
```

Include lightweight files that prove how each row was produced:

```text
args.json
eval_fixed.csv
eval_metrics.json
metrics.json
summary.json
neighbor_diagnostics*.json
neighbor_diagnostics*.csv
training_history.json
per_class_neighbor_purity.csv
```

Analysis notes:

```text
notes/
```

Include the experiment README, interpretation, decisions, and any caveats.

## What Not To Include In Git Packages

Do not commit large binary outputs to normal Git packages:

```text
*.npy
*.npz
*.h5
*.h5ad
*.pt
*.pth
*.tar
stdout.log
stderr.log
```

If large artifacts are essential, keep them in the external `results/` storage and reference their paths from the package README. Use Git LFS only when explicitly chosen for a specific artifact class.

## Required Interpretation

Do not upload only code or only tables. Every package must contain a short interpretation with this structure:

```text
Finding 1:
...

Evidence:
...

Interpretation:
...

Next action:
...
```

If the result is negative, write a negative conclusion directly. Do not frame weak or failed evidence as success.

## Git Workflow

1. Run experiments under `results/experiments/<experiment_id>/`.
2. Generate summaries under `results/experiments/<experiment_id>/summaries/`.
3. Copy lightweight evidence into `experiment_packages/<experiment_id>/`.
4. Add `MANIFEST.files` and `MANIFEST.sha256`.
5. Commit only the relevant exploratory code, package, and documentation.
6. Push the branch to GitHub.

This rule is meant to keep future exploratory code reviewable, reproducible, and separated from large runtime artifacts.
