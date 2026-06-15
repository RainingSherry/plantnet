# CutAware NeighborMix Two-Change Package

Created: 2026-06-15T13:19:59
Workspace: `/home/luolie/biopipeline/dimension-reduction/plantnet`

This package is intended for GPT-side follow-up analysis. It contains the code and lightweight results for the two NeighborMix/scMAE graph-structure modifications.

## Change Groups

1. `change_1_direct_cut_ot`
   - Variants: `canm_diagnostic_only`, `canm_cut_ot`, `canm_mix_plus_cut`, `canm_cut_ot_warm`, `canm_mix_plus_cut_warm`, `canm_attention_fusion_probe`.
   - Idea: migrate scGNN/scDSC/AttentionAE-sc/scCDCG concepts as explicit graph diagnostics, cluster head, cut loss, OT self-training, and attention fusion probe.

2. `change_2_cut_neighbor_mix`
   - Variants: `canm_cut_reweighted_mix`, plus any available gated cut variants (`canm_gated_cut_mix`, `canm_gated_cut_warm`).
   - Idea: inject the scCDCG cut-informed idea directly into NeighborMix edge weights, so candidate cross-cluster edges are downweighted before expression mixing.

## Contents

- `code/methods/DeepLearning/CutAware_NeighborMix_scMAE/`: current method implementation and configs.
- `results/tables/`: compact CSV comparison tables.
- `results/lightweight_runs/`: per-run JSON/CSV/log artifacts copied from the full local output.
- `results/full_result_pointers/omitted_large_artifacts.json`: list of omitted binary artifacts and their original local paths.
- `previous_experiment_report_snapshot/`: previously curated lightweight report snapshot, if available.
- `git/`: commit/status context at package time.

## Large Artifact Policy

The package intentionally omits `*.pt`, `*.pth`, `*.npy`, `*.npz`, `*.h5`, and `*.h5ad` files. Full local results remain under:

`results/experimental/cutaware_neighbormix_20260615`

This keeps the package small enough for GPT analysis while preserving exact local pointers for embeddings, model checkpoints, and binary graph arrays.

## Current Variant Summary

```text
                variant  n_runs  mean_ari  min_ari  max_ari
canm_cut_reweighted_mix       5  0.632087 0.373287 0.939177
   canm_diagnostic_only       5  0.605034 0.410052 0.744993
     canm_gated_cut_mix       5  0.572877 0.292344 0.938741
    canm_gated_cut_warm       5  0.571910 0.290779 0.939177
       canm_cut_ot_warm       5  0.281111 0.000000 0.785268
 canm_mix_plus_cut_warm       5  0.281013 0.000000 0.785268
            canm_cut_ot       5  0.256310 0.000000 0.733063
      canm_mix_plus_cut       5  0.255616 0.000000 0.734670
```
