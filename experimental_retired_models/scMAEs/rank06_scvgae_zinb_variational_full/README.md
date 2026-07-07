# rank06_scvgae_zinb_variational_full

Independent-full scMAE candidate inspired by scVGAE.

## Theory Basis

scVGAE combines variational graph-autoencoder ideas with ZINB loss for
single-cell count imputation. This candidate keeps scMAE as the main clustering
representation learner and adds a variational latent distribution plus a ZINB
decoder over a non-scaled count/count-like branch.

## scMAE Gap

This model addresses the distribution / VAE uncertainty gap. It preserves mask
prediction and masked log-expression reconstruction while adding `mu/logvar`,
ZINB mean, dispersion, and dropout heads.

## NeighborMix Relation

Independent and complementary. No cell mixing is used; `mixed_cell_fraction` is
always `0.0`.

