# rank03_consistency_ema_robust_full

Independent-full scMAE candidate inspired by improved consistency training.

It keeps scMAE mask prediction and masked expression reconstruction, then adds
an EMA teacher latent target. The student receives masked/noisy expression; the
teacher receives the clean view. Latent consistency uses a Pseudo-Huber metric.

This is not a full image consistency model or diffusion sampler. It is a
single-cell representation-learning adaptation of the robust denoising idea.

NeighborMix is not used; `mixed_cell_fraction` is always `0.0`.

