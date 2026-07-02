# rank49_amp_drl_mask_policy_full

Independent-full scMAE candidate adapted from **Adaptive-Masking Policy with Deep Reinforcement Learning for Self-Supervised Medical Image Segmentation**.

## Method Basis

The local report recommends this paper as a policy search method for mask/mix, with the warning that RL should remain a small controller rather than the main contribution. The PDF models masking as a reinforcement-learning problem: a policy chooses mask position and size, and reconstruction feedback acts as reward.

No GitHub URL is listed in the index, so this implementation is reconstructed from the PDF and local report.

## scMAE Gap Addressed

This candidate targets the **mask** gap:

- scMAE mask prediction is retained.
- masked expression reconstruction is retained.
- A tiny dueling-Q controller chooses among a small action space.
- Actions are `zero`, `swap`, or variance-guided mask with mask ratios `0.25`, `0.40`, `0.55`.
- The reward is the batch-level improvement in scMAE reconstruction/mask loss.

## Data Semantics

- `scaled_expr` is used only as encoder input when `--scale_input true`.
- `log_expr` is used as reconstruction target.
- Variance-guided masks use log-expression gene variance, not scaled expression as a count target.
- No NB/ZINB or generated-cell evaluation is used.

## NeighborMix Relationship

NeighborMix is not used here. The method is independent and only controls masking. Therefore `mixed_cell_fraction=0.0`.

## Differences From Original AMP-DRL

- Medical-image mask boxes are replaced by gene-wise mask strategy and ratio actions.
- A compact dueling-Q bandit-style controller is used to avoid making RL the main engineering contribution.
- The main training loss remains scMAE reconstruction plus mask prediction.

## Screen Caveat

Smoke and screen results are candidate evidence only. They must not be appended to `全benchmark结果.csv` and are not formal performance claims.
