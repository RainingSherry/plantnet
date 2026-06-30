# Rank 28: Graphormer Bias scMAE

This folder is an independent scMAE-family adaptation of:

- Paper: Graphormer
- Official code: https://github.com/microsoft/Graphormer
- Local source snapshot: `../external_sources/rank28_graphormer`
- Commit: `a04573c40705fb174db261bb746a8258d00992f5`

## scMAE Adaptation

Graphormer represents graph structure inside Transformer attention using centrality encodings, shortest-path spatial encodings, edge encodings, and a graph-token virtual distance. The scRNA benchmark has one expression vector per cell, so each mini-batch is converted into a KNN cell graph and cells are treated as graph nodes.

The adapted model keeps the core Graphormer structure:

- graph token prepended to node tokens;
- node feature projection plus in/out degree centrality embeddings;
- shortest-path spatial position embeddings per attention head;
- edge-type attention bias;
- graph-token virtual-distance bias;
- multi-head self-attention with per-head graph attention bias added before softmax;
- auxiliary masked-expression reconstruction and graph reconstruction losses for the scMAE benchmark protocol.

Mask semantics are fixed as `1 = expression feature replaced by another cell's expression value`.

## Fairness Notes

No pretrained Graphormer weights are used. The official implementation depends on Fairseq and graph molecule preprocessing; this implementation ports the structural attention-bias mechanism into dense PyTorch operations over batch-local cell graphs. Shared scMAE-family preprocessing, seed handling, and KMeans evaluation are unchanged.

## Outputs

`run.py` writes the standard independent scMAEs artifacts: `args.json`, `source_manifest.json`, `training_history.json`, `model_checkpoint.pth`, `embedding_final.npy`, `labels.npy`, `gene_names.npy`, `embedding.h5`, and `eval_fixed.csv`/`metrics.json` when evaluation is enabled.
