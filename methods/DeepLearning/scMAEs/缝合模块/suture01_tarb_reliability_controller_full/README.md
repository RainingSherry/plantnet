# suture01_tarb_reliability_controller_full

Independent-full scMAE candidate inspired by the read-only suture module:

`缝合模块/即插即用/225 Task Adaptive Restoration Block(TPAMI 2026).py`

The original module adaptively selects image restoration operations. This
candidate rewrites that idea for scRNA-seq as a latent operation controller.
The scMAE full-gene encoder, mask prediction, masked expression reconstruction,
and DEC-style cluster centers remain the main path.

## Mechanism

- `scaled_expr` is used only as encoder input when enabled.
- `log_expr` is used as reconstruction target.
- A small controller predicts weights for five latent operations.
- Risky smoothing-like operations are weakened on low-reliability cells.
- DEC KL is confidence-gated only; reliability does not gate DEC.

## Why This Design

Prior GatedNeighborMix and DualAxis experiments showed that the safe approach is
to keep rank13-style scMAE/DEC geometry intact and only gate risky auxiliary
operations. Replacing the encoder with compressed gene-module attention can lose
fine marker signals.

## Outputs

Each run saves `embedding_final.npy`, `labels.npy`, `training_history.json`,
`diagnostics.json`, `summary.json`, `operation_weights.npy`, `reliability.npy`,
and optional evaluation outputs.

Smoke and screen results are candidate evidence only and are not appended to
`全benchmark结果.csv`.
