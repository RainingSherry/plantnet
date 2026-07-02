# suture12_attention_matching_teacher_consistency_full

This independent-full candidate rewrites Attention Matching for scRNA masked autoencoding.

## Mechanism

The original module matches support/query features using cosine attention. This version matches a masked latent view to a detached clean latent view from the same cell. The clean view is used only as a teacher context, not as labels.

The scMAE backbone remains intact:

- mask prediction;
- masked expression reconstruction.

The added branch is intentionally weak and can be disabled with `--match_weight 0` and `--match_loss_weight 0`.

## Gap Addressed

This candidate targets teacher/robust-loss consistency without graph propagation, prototype memory, or NeighborMix. It is intended to avoid the strong latent geometry shifts seen in previous failed candidates.

## NeighborMix

NeighborMix is not used. There is no cell mixing, and diagnostics report `mixed_cell_fraction=0.0`.

## Source

Reference module:

`/home/luolie/biopipeline/dimension-reduction/plantnet/缝合模块/即插即用/131 AttentionMacthcing(AAAI 2025).py`

Implementation source: rewritten from the mechanism description for clean/masked scRNA latent matching.
