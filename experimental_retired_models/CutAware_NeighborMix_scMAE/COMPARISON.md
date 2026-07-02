# Graph Methods To NeighborMix Migration

## Core Diagnosis

NeighborMix and most GNN-style scRNA clustering methods share one risky assumption:

```text
If two cells are neighbors, making their representations closer is usually helpful.
```

That assumption breaks at cluster boundaries. A KNN graph inevitably contains cross-type edges, especially in transitional, rare-cell, and imbalanced datasets. If the training objective treats all neighbor edges as smoothing signals, the model can erase the exact boundaries needed for clustering.

CutAware_NeighborMix_scMAE changes the graph role:

```text
The graph is not a smoother.
The graph is a cut constraint over soft cluster assignments.
```

## Method-by-method Comparison

### scGNN

scGNN solves sparsity and cell-cell structure by combining feature autoencoding, graph autoencoding, iterative graph construction/pruning, and LTMG-based expression-state regularization.

Transferable parts:

```text
AE first, graph second
Build graph from learned/denoised embedding rather than raw expression only
Prune weak/noisy graph edges
Optionally refresh the graph during training
```

Implementation in this folder:

```text
build_pca_knn_graph(...)
build_embedding_knn_graph(...)
--edge_prune_quantile
--graph_refresh_interval
```

Not transferred:

```text
GNN message passing
LTMG
cluster/imputation autoencoder loop
```

Reason: the current failure mode is not lack of graph propagation; it is insufficient boundary discrimination.

### scDSC

scDSC addresses target mismatch by jointly training a ZINB-based AE, GNN structure module, and clustering supervision.

Transferable parts:

```text
Do not rely on reconstruction alone
Add an explicit clustering-aware objective
Train expression representation and graph/clustering objective together
```

Implementation in this folder:

```text
CutAwareAutoEncoder.cluster_head
graph_cut_loss(...)
ot_self_training_loss(...)
```

Not transferred:

```text
ZINB decoder
Full GNN structure module
Mutual-supervised GNN pipeline
```

Reason: scMAE-family already has a stable masked reconstruction target; the missing piece is a boundary-aware clustering target.

### AttentionAE-sc

AttentionAE-sc fuses denoising and topological embeddings with multi-head attention.

Transferable parts:

```text
Keep expression-denoising representation separate from graph-derived representation
Measure whether graph assignment agrees with expression assignment
Use graph information conditionally, not blindly
```

Implementation in this folder:

```text
canm_attention_fusion_probe
attention_fusion_probe(...)
```

Not transferred as main mechanism:

```text
Full attention fusion backbone
Full GAE branch
```

Reason: attention can reweight graph information, but it does not guarantee that bad cross-cluster edges are cut. It is useful as a probe, not as the first fix.

### scCDCG / DCGC

scCDCG reframes graph use as cut-informed graph embedding with OT-guided self-supervised clustering. This is the closest match to the NeighborMix failure mode.

Transferable parts:

```text
Use graph cut rather than graph smoothing
Train soft cluster assignments directly
Use balanced OT-style assignments to avoid clustering collapse
Track embedding distinguishability
```

Implementation in this folder:

```text
graph_cut_loss(...)
sinkhorn_balanced_assignment(...)
ot_self_training_loss(...)
apply_cluster_cut_reweight(...)
CutAwareAutoEncoder.edge_gate_scores(...)
cut_diagnostics.json
embedding_similarity_diagnostics.json
```

Not transferred:

```text
Full scCDCG/DCGC architecture
High-order cut-informed encoder stack
Dataset-specific training tricks
```

Reason: the first experiment should isolate whether cut/OT fixes NeighborMix negative transfer before adding a larger architecture.

## Expected Comparison

| Variant | Neighbor information use | Boundary signal | Collapse guard | Main question |
| --- | --- | --- | --- | --- |
| scMAE | None | None | KMeans only | How good is expression-only masked AE? |
| NeighborMix | Expression mixing | None | None | Does neighbor perturbation help? |
| RG NeighborMix | Reliability-weighted mixing | Indirect | Gate diagnostics only | Does reliability avoid bad mixing? |
| canm_diagnostic_only | Graph diagnostics only | None in training | Measured only | Is the graph already unsafe? |
| canm_cut_ot | Graph cut objective | Explicit | Sinkhorn balance | Does cut-aware graph use beat mixing? |
| canm_mix_plus_cut | Mixing + cut objective | Explicit | Sinkhorn balance | Is mixing still useful after boundary control? |
| canm_cut_reweighted_mix | Mixing after cross-edge downweighting | Explicit in edge weights | Avoids direct cut-gradient collapse | Does scCDCG-style cutting fix NeighborMix itself? |
| canm_gated_cut_mix | Learnable gated mixing after cross-edge downweighting | Explicit plus learned edge reliability | Gate prior/entropy/cluster consistency | Can attention/gating recover useful neighbors after cut priors? |
| canm_gated_cut_warm | Warm-started gated mixing after cross-edge downweighting | Explicit plus learned edge reliability | Delayed gate, no cluster gate loss | Does a conservative gate help without amplifying early pseudo-cluster noise? |
| canm_attention_fusion_probe | Assignment fusion probe | Indirect | Sinkhorn/cut optional | Is attention disagreement informative? |

## Recommended First Run

Use the Phase-2 representative datasets where RG NeighborMix was already inspected:

```text
SRP182008
Melanoma_5K
Macosko
Tosches
Wang
```

Run only seed 42 first:

```text
canm_diagnostic_only
canm_cut_ot
canm_mix_plus_cut
canm_cut_reweighted_mix
```

Then compare against:

```text
results/formal/rg_phase2_sensitivity_e80/rg_phase2_all_sweeps_raw.csv
results/formal/rg_phase1_allseeds_e80/rg_phase1_allseeds_summary.csv
```

## Interpretation Rules

```text
canm_cut_ot beats canm_diagnostic_only:
  Cut/OT contributes beyond passive graph diagnostics.

canm_cut_ot beats NeighborMix/RG on negative datasets:
  Neighbor information should be used as a cut constraint, not as expression averaging.

canm_mix_plus_cut beats canm_cut_ot:
  Expression mixing is not fundamentally wrong, but it needs boundary-aware constraints.

canm_mix_plus_cut underperforms canm_cut_ot:
  NeighborMix perturbation is likely adding harmful invariance pressure.

canm_cut_reweighted_mix beats canm_mix_plus_cut:
  The useful migration is to cut/downweight bad NeighborMix edges, not to add another global clustering loss.

cluster_mass_max is high or fraction_cosine_gt_0p9 is high:
  The method is collapsing despite OT; increase OT weight, lower cut weight, or improve graph pruning.
```
