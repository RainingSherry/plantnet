# Rank 20: MultiMAE Targets scMAE

This folder is an independent scMAE-family adaptation of:

- Paper: MultiMAE: Multi-modal Multi-task Masked Autoencoders
- Official code: https://github.com/EPFL-VILAB/MultiMAE
- Local source snapshot: `../external_sources/rank20_multimae`
- Commit: `66910f5b5ba236f5e731883db85fe4f24ee01106`

## scMAE Adaptation

The benchmark datasets are single-modality scRNA-seq matrices, so this implementation does not invent unavailable ADT/ATAC modalities. Instead it maps the MultiMAE multi-task design onto three RNA-internal tasks:

- `expr`: scaled expression patch values.
- `rank`: per-cell gene rank percentiles, patchified on the same HVGs.
- `stat`: nonnegative patch statistics: mean, standard deviation, and zero fraction.

The model keeps the MultiMAE core mechanics:

- one input adapter per task;
- task and position embeddings;
- Dirichlet sampling of a fixed number of visible tokens across tasks;
- a shared Transformer encoder that receives only visible tokens plus global tokens;
- task-specific decoders that unshuffle the full token sequence, insert learned mask tokens, and cross-attend from task queries to encoder tokens;
- loss only on masked task tokens.

Mask semantics are fixed as `1 = removed/masked task token to reconstruct`, matching the original MultiMAE convention.

## Outputs

`run.py` writes the standard independent scMAEs artifacts:

- `args.json`
- `source_manifest.json`
- `dataset_profile.json`
- `preprocess_config.json`
- `training_history.json`
- `model_checkpoint.pth`
- `embedding_final.npy`
- `labels.npy`
- `gene_names.npy`
- `embedding.h5`
- `metrics.json` and `eval_fixed.csv` when evaluation is enabled

## Fairness Notes

Preprocessing remains the shared scMAE-family protocol through `methods.DeepLearning.scMAE_family`: HVG selection, normalization, optional scaling, and KMeans known-k evaluation are not changed. The auxiliary `stat` task uses the same selected genes with `scale_input=False` only to compute nonnegative patch statistics.
