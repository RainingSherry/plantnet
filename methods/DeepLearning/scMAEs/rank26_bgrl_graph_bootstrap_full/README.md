# Rank 26: BGRL Graph Bootstrap scMAE

This folder is an independent scMAE-family adaptation of:

- Paper: BGRL: Large-Scale Representation Learning on Graphs via Bootstrapping
- Official code: https://github.com/nerdslab/bgrl
- Local source snapshot: `../external_sources/rank26_bgrl`
- Commit: `60f9f19ad0598f9163ad70ebbde3e7297760901e`

## scMAE Adaptation

BGRL learns node representations from two augmented graph views with an online encoder, an EMA target encoder, and an online predictor. The scRNA benchmark has one cell-expression vector per node, so each mini-batch is converted to a batch-local KNN cell graph and two stochastic graph views are generated.

The adapted model keeps the core BGRL structure:

- online graph encoder and independently initialized EMA target encoder;
- MLP predictor on the online branch;
- two augmented views with edge dropout and feature masking;
- symmetric cosine bootstrap loss from online predictions to stopped-gradient target embeddings;
- dependency-free dense GraphSAGE-style graph encoder over batch-local cell graphs;
- auxiliary masked-expression reconstruction and mask prediction heads for the scMAE benchmark protocol.

Mask semantics are fixed as `1 = expression feature set to zero by BGRL view corruption`.

## Fairness Notes

No pretrained graph weights are used. The official code depends on PyTorch Geometric; this implementation ports the BGRL architecture and augmentations into dense PyTorch operations so it can run inside the existing benchmark environment without changing dependencies, preprocessing, or evaluation.

## Outputs

`run.py` writes the standard independent scMAEs artifacts: `args.json`, `source_manifest.json`, `training_history.json`, `model_checkpoint.pth`, `embedding_final.npy`, `labels.npy`, `gene_names.npy`, `embedding.h5`, and `eval_fixed.csv`/`metrics.json` when evaluation is enabled.
