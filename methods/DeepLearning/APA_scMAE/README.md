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

## v2 objective correction

APA-scMAE v2 keeps the same high-level Student/Generator structure, but changes the training objective so the cell embedding is directly shaped:

- Student warmup is enabled by default for 20 epochs; warmup uses random scMAE-style masks and does not update the Generator.
- Student training includes clean/masked representation consistency, variance/covariance anti-collapse penalties, EMA Teacher consistency, and embedding-prototype assignment consistency.
- The default decoder mode is `z_with_stopgrad_h`, so reconstruction may use token context but cannot route reconstruction gradients through token features to bypass the cell embedding.
- The Generator objective targets moderate corruption difficulty plus mask coverage/diversity/gene balance instead of rewarding the largest reconstruction deltas.

Useful one-command ablations:

```bash
--use_repr_loss false
--use_ema_teacher false
--use_proto_consistency false
--student_warmup_epochs 0
--decoder_mode current
```

Notes:

- Do not add APA-scMAE to `methods/method_manifest.yaml` until a separate review approves formal benchmark inclusion.
- Prototypes are label-free; true labels may only be used by final evaluation code.
- `kmeans_known_k` uses the supplied `n_clusters` and must be treated as known-K / oracle-K evaluation, not fully unsupervised unknown-K clustering.
- When `--skip_eval true`, the runner may still save prediction labels from KMeans, but it does not report supervised ACC/NMI/ARI/F1 metrics and does not save `labels.npy`.
