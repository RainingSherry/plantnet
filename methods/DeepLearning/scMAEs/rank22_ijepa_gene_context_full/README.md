# Rank 22: I-JEPA Gene-Context scMAE

This folder is an independent scMAE-family adaptation of:

- Paper: I-JEPA: Self-Supervised Learning from Images with Joint-Embedding Predictive Architecture
- Official code: https://github.com/facebookresearch/ijepa
- Local source snapshot: `../external_sources/rank22_ijepa`
- Commit: `52c1ae95d05f743e000e8f10a1f3a79b10cff048`

## scMAE Adaptation

I-JEPA predicts target-region latent representations from a surrounding context representation, using an EMA target encoder and a predictor network with mask tokens. The scRNA input is a 1D expression vector, so this implementation splits genes into contiguous patches instead of image patches.

The adapted model keeps the core I-JEPA structure:

- online context encoder over expression patches;
- EMA target encoder over the full expression vector;
- non-overlapping context and target patch masks;
- predictor with learned target tokens and positional embeddings;
- SmoothL1 latent loss only on target patches;
- auxiliary scMAE reconstruction and target-mask prediction heads for benchmark-compatible masked-expression training.

Mask semantics are fixed as `1 = target patch predicted/reconstructed from context`.

## Fairness Notes

No image checkpoint, ImageNet data, or external pretrained representation is used. Preprocessing, KMeans evaluation, seed handling, and output files stay inside the shared scMAE-family protocol. The method-specific architecture and losses are contained in this directory.

## Outputs

`run.py` writes the standard independent scMAEs artifacts: `args.json`, `source_manifest.json`, `training_history.json`, `model_checkpoint.pth`, `embedding_final.npy`, `labels.npy`, `gene_names.npy`, `embedding.h5`, and `eval_fixed.csv`/`metrics.json` when evaluation is enabled.
