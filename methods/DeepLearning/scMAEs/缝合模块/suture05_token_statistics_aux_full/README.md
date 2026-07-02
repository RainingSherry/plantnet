# suture05_token_statistics_aux_full

This candidate adapts Token Statistics Self-Attention into a gene-token statistics auxiliary branch for scMAE.

Each gene is treated as a token with a learned embedding scaled by expression. The branch computes linear token statistics and injects a small auxiliary latent into the full-gene scMAE encoder output. It keeps mask prediction and masked expression reconstruction.

No gene vector is reshaped into an image. `scaled_expr` is only encoder input and `log_expr` is the reconstruction/statistics target. NeighborMix is not used.
