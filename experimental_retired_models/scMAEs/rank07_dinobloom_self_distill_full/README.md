# rank07_dinobloom_self_distill_full

Independent-full scMAE candidate inspired by DinoBloom/DINOv2.

The model keeps scMAE mask prediction and masked expression reconstruction. It
adds a DINO-style projection/prototype head, an EMA teacher, center/sharpened
teacher probabilities, and collapse diagnostics. It does not use external
DinoBloom image weights.

NeighborMix is not used; `mixed_cell_fraction` is always `0.0`.

