# rank38_latent_diffusion_denoise_full

Independent-full scMAE candidate adapted from **Latent Diffusion Models**.

## Method Basis

The local report recommends Latent Diffusion as a lightweight generative / denoising auxiliary rather than a heavy generator. The paper proposes training diffusion models in a learned latent space produced by an autoencoder, using a forward noising process and a denoising objective. The official GitHub implementation was inspected for `DDPM`, `LatentDiffusion`, `q_sample`, `predict_start_from_noise`, `p_losses`, and epsilon parameterization.

## scMAE Gap Addressed

This candidate targets the **robust loss / latent target** gap:

- scMAE mask prediction is retained.
- masked expression reconstruction is retained.
- the scMAE encoder provides a compact latent space;
- a lightweight timestep-conditioned denoiser predicts DDPM epsilon noise in latent space;
- an estimated `z0_hat` is matched to the clean scMAE latent as a consistency term.

Generated latent samples are not treated as real cells and are not used for benchmark evaluation.

## Data Semantics

- `scaled_expr` is used only as encoder input when `--scale_input true`.
- `log_expr` is used as masked expression reconstruction target.
- No NB/ZINB count likelihood is used.
- No generated expression or synthetic cell is evaluated as a real sample.

## NeighborMix Relationship

NeighborMix is not used. This candidate is independent and potentially complementary. `mixed_cell_fraction=0.0`.

## Screen Caveat

Smoke and screen results are candidate evidence only. They must not be appended to `全benchmark结果.csv` and are not formal performance claims.
