# Rank 21: AudioMAE Local-Window scMAE

This folder is an independent scMAE-family adaptation of:

- Paper: AudioMAE: Masked Autoencoders that Listen
- Official code: https://github.com/facebookresearch/AudioMAE
- Local source snapshot: `../external_sources/rank21_audiomae`
- Commit: `bd60e29651285f80d32a6405082835ad26e6f19f`

## scMAE Adaptation

AudioMAE is designed for 2D audio spectrograms. The scRNA benchmark has a single expression vector per cell, so this implementation maps the expression vector to a padded 2D gene grid before patch embedding. It does not use audio checkpoints or audio feature extraction.

The adapted model keeps the core AudioMAE structure:

- Conv2d patch embedding over a 2D grid;
- high-ratio masking;
- optional 2D row/column masking analogous to time/frequency masking;
- MAE encoder over visible patches plus a class token;
- decoder unshuffle with learned mask tokens;
- local-window decoder attention with alternating shifted windows;
- reconstruction loss only on removed patches.

Mask semantics are fixed as `1 = removed/masked patch`.

## Fairness Notes

Preprocessing and evaluation stay inside the shared scMAE-family protocol. Padded gene-grid positions are excluded from the reconstruction loss denominator, so loss normalization uses only real gene values inside masked patches.

## Outputs

`run.py` writes the standard independent scMAEs artifacts: `args.json`, `source_manifest.json`, `training_history.json`, `model_checkpoint.pth`, `embedding_final.npy`, `labels.npy`, `gene_names.npy`, `embedding.h5`, and `eval_fixed.csv`/`metrics.json` when evaluation is enabled.
