# CutAware_NeighborMix_scMAE

This is an experimental migration of graph-structured single-cell clustering ideas into the NeighborMix/scMAE family.

The central change is conceptual:

```text
NeighborMix old question:
  How should a cell be mixed with its neighbors?

CutAware question:
  Which graph edges should be preserved inside clusters, and which edges should be cut across clusters?
```

The default route therefore does not deepen a GNN and does not add attention as the main mechanism. It keeps the scMAE expression encoder and adds a cluster head trained by a mini-batch cut-informed loss plus Sinkhorn-balanced self-training.

## Why This Exists

The RG NeighborMix experiments showed stable gates and non-collapsed edge weights, but the mechanism did not beat random controls robustly. That suggests the failure is not just an implementation problem. The training target still treats neighbor perturbation as something the model should reconstruct away, so it has no explicit reason to preserve cluster boundaries.

scCluBench makes the same distinction for graph methods: graph methods have structural awareness, but GNN-style message passing often causes over-smoothing and embedding collapse. scCDCG is the exception because it uses cut-informed graph embedding rather than ordinary neighborhood smoothing.

## Migrated Ideas

| Source method | What it solves | What is migrated here | What is deliberately not migrated |
| --- | --- | --- | --- |
| scGNN | Learns an AE embedding, builds/prunes a cell graph, then uses graph learning | PCA/embedding KNN graph, optional edge pruning, optional graph refresh | Full GNN message passing and LTMG |
| scDSC | Adds clustering supervision to AE + graph learning | Explicit cluster head and clustering loss | ZINB decoder and GNN propagation |
| AttentionAE-sc | Fuses denoising and graph-topology embeddings | `canm_attention_fusion_probe` diagnostic variant | Attention as the main solution |
| scCDCG/DCGC | Reframes graph use as cut-informed clustering with OT-balanced assignment | Mini-batch normalized-cut surrogate, Sinkhorn-balanced assignment, and cut-reweighted NeighborMix edges | Full published architecture |

## Variants

```text
canm_diagnostic_only
  scMAE training only. Saves graph cut, assignment, and embedding collapse diagnostics.

canm_cut_ot
  Default. scMAE + cut-aware graph clustering loss + Sinkhorn-balanced self-training.

canm_mix_plus_cut
  NeighborMix expression perturbation + cut/OT losses. This tests whether mixing still helps once the boundary objective is explicit.

canm_cut_reweighted_mix
  NeighborMix expression perturbation with candidate cross-cluster edges explicitly downweighted before mixing. This is the direct network-level migration of the scCDCG idea into NeighborMix.

canm_attention_fusion_probe
  Diagnostic AttentionAE-sc-style fusion probe. It is intentionally not the main route because attention does not by itself cut bad edges.
```

## Example

Do not use GPU 0 or 7. GPU 1 is the default.

```bash
python methods/DeepLearning/CutAware_NeighborMix_scMAE/run.py \
  --data_path data/processed_scmae/Macosko.h5ad \
  --save_dir results/experimental/cutaware_nm/Macosko/seed42/canm_cut_ot \
  --dataset_name Macosko \
  --variant_name canm_cut_ot \
  --seed 42 \
  --n_clusters 12 \
  --gpu 1 \
  --epochs 80 \
  --batch_size 256 \
  --neighbor_k 10 \
  --cut_weight 0.2 \
  --ot_weight 0.1 \
  --no_save_h5ad
```

## Outputs

The runner writes the same evaluation style as scMAE-family methods:

```text
eval_fixed.csv
eval_metrics.json
metrics.json
summary.json
embedding_final.npy
training_history.json
```

Additional cut-aware diagnostics:

```text
cut_diagnostics.json
embedding_similarity_diagnostics.json
cluster_probs.npy
neighbor_graph_profile.json
edge_weight_summary.json
```

Key diagnostics to compare against NeighborMix:

```text
ARI / ACC / NMI
cut_diagnostics.ncut_surrogate
cut_diagnostics.cluster_mass_min / cluster_mass_max
embedding_similarity_diagnostics.fraction_cosine_gt_0p9
embedding_geometry_summary.between_within_ratio
rare_cell_effect_summary
```

## First Ablation

Run these on the same datasets and seeds:

```text
scMAE baseline
NeighborMix_scMAE or RG/RC NeighborMix best current variant
canm_diagnostic_only
canm_cut_ot
canm_mix_plus_cut
canm_cut_reweighted_mix
```

Decision rule:

```text
If canm_cut_ot > NeighborMix on negative-transfer datasets:
  The problem was the use of neighbor information, not absence of neighbor information.

If canm_mix_plus_cut > canm_cut_ot:
  Expression mixing can remain, but only after a boundary-aware objective is present.

If canm_cut_reweighted_mix > canm_mix_plus_cut:
  Cutting candidate cross-cluster edges before NeighborMix is safer than adding an unstable cut/OT gradient.

If canm_diagnostic_only already shows high collapse or high bad-cut diagnostics:
  Graph construction/pruning must be fixed before training with cut loss.
```

## References

- scGNN: https://www.nature.com/articles/s41467-021-22197-x and https://github.com/juexinwang/scGNN
- scDSC: https://doi.org/10.1093/bib/bbac018
- AttentionAE-sc: https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1011641 and https://github.com/LiShenghao813/AttentionAE-sc
- scCDCG: https://arxiv.org/abs/2404.06167 and https://github.com/XPgogogo/scCDCG
- scCluBench: https://arxiv.org/abs/2512.02471
