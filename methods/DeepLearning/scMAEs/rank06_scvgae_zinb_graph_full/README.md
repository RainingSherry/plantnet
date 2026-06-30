# Rank 06: scVGAE ZINB Graph scMAE

Source paper: **scVGAE: ZINB-Based Variational Graph Autoencoder for
Single-Cell RNA-Seq Imputation**.

The GitHub URL listed in the index could not be cloned anonymously in this
environment, so this implementation is based on the local PDF and the summary
report. The paper describes a GCN encoder, graph branches for ZINB mean/dropout
probability/dispersion, reconstruction loss, and joint optimization.

This independent scMAE-family variant implements:

- batch-local KNN cell graph with normalized adjacency;
- GCN encoder with GraphNorm;
- variational latent mean/log-variance and reparameterization;
- graph heads for ZINB `mu`, `theta`, and dropout logits `pi`;
- masked scMAE reconstruction branch;
- combined ZINB, masked reconstruction, KL, and mask-prediction losses.

Mask semantics: `1 = expression gene replaced for masked reconstruction branch`.
The masked reconstruction denominator is the number of effectively replaced gene
entries in the batch, clamped only for degenerate smoke tests.

ZINB target handling:

- the encoder input and evaluation embedding use the same fair scMAE protocol as
  the other variants (`n_top_genes=1000`, `target_sum=10000`,
  `scale_input=True` by default);
- because ZINB likelihood requires nonnegative observations, the ZINB branch uses
  the same selected HVG genes loaded with `scale_input=False` as a nonnegative
  log-normalized expression target;
- this auxiliary target is recorded in `preprocess_config.json` and is not used
  for KMeans evaluation.

Not reproduced:

- the original repository code could not be downloaded anonymously from the
  indexed GitHub URL;
- imputation-only postprocessing from the paper is not used because the target
  task here is fair clustering from learned embeddings.
