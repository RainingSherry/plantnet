# APA-scMAE

APA-scMAE is a prototype-aware adversarial masked autoencoder built on the scMAE training idea.

The implementation follows the attached model sketch:

- labels are loaded only for final evaluation;
- gene statistics are `[mean, variance, zero_rate, hvg_rank]`;
- prototypes are built by label-free PCA plus KMeans;
- the corruption module supplies replacement values `V`;
- the generator only selects mask positions and never generates replacement values;
- `M` is the selected hard mask, `E` marks positions where `V` differs from `X`, and `M_eff = M * E` is the effective-mask supervision target;
- embeddings are extracted from the clean Student path without Generator, mask, or donor replacement.

This package is not registered in the formal benchmark manifest yet.
