# RG_NeighborMix_scMAE

Reliability-Gated NeighborMix-scMAE is an isolated proposed variant that tests whether fixed NeighborMix can be made less data-structure dependent.

Core form:

```text
x_i_prime = (1 - g_i) x_i + g_i sum_j a_ij x_j
```

- `a_ij` is an analytic edge reliability weight.
- `g_i` is an unsupervised node-level gate.
- Phase 1 supports `none`, `fixed`, `mutual`, `reliability`, `random`, and `far`.
- `attention` is reserved for phase 2 and requires `--allow_attention_phase2 true`.

The original `methods/DeepLearning/NeighborMix_scMAE/` and sibling `RA_NeighborMix_scMAE/` directories are not modified by this method.

## Main Command

```bash
python -m methods.DeepLearning.RG_NeighborMix_scMAE.run \
  --data_path data/processed/Macosko.h5ad \
  --save_dir results/RG_NeighborMix_scMAE/Macosko/seed42/reliability \
  --dataset_name Macosko \
  --method_name RG_NeighborMix_scMAE \
  --variant_name rg_nm_v1_reliability \
  --mix_mode reliability \
  --gate_mode topology \
  --edge_reliability_mode sim_mutual_snn_distance \
  --gpu 1
```

## Canonical Outputs

```text
embedding_final.npy
embeddings_base.npy
labels.npy
gene_names.npy
embedding.h5
metrics.json
args.json
summary.json
training_history.json
model.pt
```

Additional diagnostics include `neighbor_indices.npy`, `edge_reliability.npy`, `node_gate.npy`, `edge_weight_summary.json`, `gate_summary.json`, `embedding_geometry_summary.csv`, and `per_cell_type_metrics.csv`.
