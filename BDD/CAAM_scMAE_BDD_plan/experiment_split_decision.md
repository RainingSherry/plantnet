# CAAM-scMAE Experiment Split Decision

## Smoke Datasets

- `data/smoke/subsample_2k.h5ad`: formal smoke only; may be regenerated locally and must not be treated as a benchmark result.
- `data/smoke/SRP182008_200x500.h5ad`: tiny runner/debug smoke only; not for model claims.

## Development Datasets

- `data/processed/Quake_Smart-seq2_Lung.h5ad`: allowed for debugging, failure analysis, and small internal ablation.
- `data/其他/Mouse_Pancreas_1.h5ad`: allowed for debugging, failure analysis, and small internal ablation.
- `data/processed_scmae/Limb_Muscle.h5ad`: optional medium development dataset if runtime allows.

## Validation Datasets

- `data/processed_scmae/Young.h5ad`
- `data/processed_scmae/Baron.h5ad`
- `data/SRP182008.h5ad`

## Sealed Test Datasets

- `data/CRA007122.h5ad`
- `data/SRP224648.h5ad`
- `data/SRP235541.h5ad`
- `data/SRP145013.h5ad`
- `data/SRP309176.h5ad`
- `data/processed_benchmark/PRJNA895163.h5ad`

Do not tune model structure, losses, donor logic, mask logic, or hyperparameters based on sealed test results.
