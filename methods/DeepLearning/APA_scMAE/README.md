# APA-scMAE

APA-scMAE is a prototype-aware adversarial masked autoencoder built on the scMAE training idea.

This is an experimental development method, not a verified formal-benchmark method.

The implementation follows the attached model sketch:

- labels are loaded only for final evaluation;
- gene statistics are `[mean, variance, zero_rate, hvg_rank]`;
- prototypes are built by label-free PCA plus KMeans;
- the corruption module supplies replacement values `V`;
- the generator only selects mask positions and never generates replacement values;
- `M` is the selected hard mask, `E` marks positions where `V` differs from `X`, and `M_eff = M * E` is the effective-mask supervision target;
- embeddings are extracted from the clean Student path without Generator, mask, or donor replacement.

This package is not registered in the formal benchmark manifest yet.

Notes:

- Do not add APA-scMAE to `methods/method_manifest.yaml` until a separate review approves formal benchmark inclusion.
- Prototypes are label-free; true labels may only be used by final evaluation code.
- `kmeans_known_k` uses the supplied `n_clusters` and must be treated as known-K / oracle-K evaluation, not fully unsupervised unknown-K clustering.
- When `--skip_eval true`, the runner may still save prediction labels from KMeans, but it does not report supervised ACC/NMI/ARI/F1 metrics and does not save `labels.npy`.
