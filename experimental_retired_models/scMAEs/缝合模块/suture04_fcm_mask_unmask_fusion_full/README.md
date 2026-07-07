# suture04_fcm_mask_unmask_fusion_full

This candidate adapts Feature Correction Module ideas into scMAE as masked/unmasked latent correction fusion.

The original FCM targets multi-source remote sensing features and uses spatial/channel fusion. Here it is rewritten as a latent channel correction module: the masked-cell latent is corrected toward a detached clean-view latent teacher with a learned per-dimension gate. It keeps scMAE mask prediction and masked expression reconstruction.

`scaled_expr` is only encoder input. `log_expr` is reconstruction target. No image reshape, no cell mixing, and no NeighborMix are used.

## Recommended screen setting

The first screen found that an unconstrained FCM gate can shut down to nearly zero. The current `gate06` setting keeps a weak but non-zero correction:

```bash
--module_weight 0.25 --min_gate 0.02 --target_gate 0.06 --gate_budget_weight 0.03
```

This passed the seed-42 screen gate on Melanoma_5K and Macosko by ARI, while Quake_10x_Spleen remained below baseline. Formal three-seed evaluation is still required before appending claims to the full benchmark.
