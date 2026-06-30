# Rank 03: Consistency Models scMAE

Source paper: **Improved Techniques for Training Consistency Models** (ICLR 2024).

This folder is an independent scMAE-family variant. It does not call the legacy
`scMAEs/common/model.py` switchboard. The implementation adapts the consistency
training objective to single-cell expression vectors:

- expression vectors are corrupted with masked Gaussian noise at Karras sigma levels;
- the online denoiser predicts clean expression at sigma `t`;
- an EMA target denoiser predicts the adjacent lower-noise state `t_next`;
- the consistency target follows the OpenAI consistency-model Euler step when no
  external score teacher is available;
- Pseudo-Huber distance is used for the consistency and reconstruction terms to
  reduce outlier/dropout domination;
- a mask prediction head preserves the scMAE masked-corruption signal.

Mask semantics: `1 = expression gene corrupted by masked Gaussian sigma noise`.
Loss denominators for consistency and reconstruction are the number of corrupted
gene entries in the batch, clamped only to avoid degenerate smoke-test division by
zero.

Fair protocol shared with the other independent variants:

- preprocessing uses the repository `scMAE_family` loader;
- default `n_top_genes=1000`, `target_sum=10000`, `scale_input=True`;
- evaluation writes KMeans known-k metrics through the same fixed protocol;
- no GPU 0 or 7 should be used by the benchmark runner.

Not reproduced:

- no image UNet, LPIPS metric, or external pretrained score model is used because
  those are not meaningful for scRNA expression clustering and would violate the
  no-large-external-weight assumption;
- the core consistency ingredients retained here are Karras sigma scheduling,
  boundary-condition denoising scalings, adjacent-time consistency, EMA target,
  and robust Pseudo-Huber regression.
