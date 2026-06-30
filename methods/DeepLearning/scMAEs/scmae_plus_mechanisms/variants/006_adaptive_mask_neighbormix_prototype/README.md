# 006 Adaptive Mask NeighborMix Prototype

This variant is the planned second-batch combination from the mechanism-search
objective. It keeps the original scMAE backbone as the main task and combines
three local regularizers:

```text
scMAE + variance-adaptive replacement mask + first reliable NeighborMix + delayed DEC prototype
```

All auxiliary losses warm up from zero. The scMAE reconstruction and mask
prediction losses remain the dominant objective.

This variant should only be promoted if its three-seed Melanoma_5K screen is
both above the original scMAE reference and more stable than standalone
NeighborMix.

## Melanoma_5K Screen Result

Completed on 2026-06-30 with `seeds = 42, 2024, 3407`, `epochs = 80`,
`n_top_genes = 1000`, and `batch_size = 128`.

```text
ARI mean = 0.654841
ARI std  = 0.016507
ACC mean = 0.724426
```

This is below the original scMAE Melanoma_5K ARI reference of about `0.668`.
The combination is therefore not promoted. The result suggests that simply
stacking adaptive masking and DEC prototypes on top of NeighborMix can dampen
the useful NeighborMix signal rather than stabilize it.
