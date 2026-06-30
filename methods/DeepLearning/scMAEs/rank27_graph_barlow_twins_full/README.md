# Rank 27: Graph Barlow Twins scMAE

This folder is an independent scMAE-family adaptation of:

- Paper: Graph Barlow Twins
- Official code: https://github.com/pbielak/graph-barlow-twins
- Local source snapshot: `../external_sources/rank27_graph_barlow_twins`
- Commit: `ec62580aa89bf3f0d20c92e7549031deedc105ab`

## scMAE Adaptation

Graph Barlow Twins trains graph encoders with two augmented views and a Barlow Twins cross-correlation loss. The scRNA benchmark has one expression vector per cell, so each mini-batch becomes a KNN cell graph and two stochastic graph views are generated.

The adapted model keeps the core Graph Barlow structure:

- two graph views with edge dropout and feature masking;
- shared graph encoder for both views;
- projection head before the Barlow objective;
- batch-normalized cross-correlation matrix;
- diagonal invariance loss and off-diagonal redundancy reduction;
- auxiliary masked-expression reconstruction and mask prediction heads for the scMAE benchmark protocol.

Mask semantics are fixed as `1 = expression feature zeroed by graph-view corruption`.

## Fairness Notes

No pretrained graph model or external labels are used. The official code depends on PyTorch Geometric; this implementation ports the Graph Barlow objective and graph augmentations into dense PyTorch operations so preprocessing, dependencies, and KMeans evaluation remain unchanged.

## Outputs

`run.py` writes the standard independent scMAEs artifacts: `args.json`, `source_manifest.json`, `training_history.json`, `model_checkpoint.pth`, `embedding_final.npy`, `labels.npy`, `gene_names.npy`, `embedding.h5`, and `eval_fixed.csv`/`metrics.json` when evaluation is enabled.
