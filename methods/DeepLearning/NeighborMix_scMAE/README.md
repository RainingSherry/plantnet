# NeighborMix-scMAE

NeighborMix-scMAE keeps the original scMAE `AutoEncoder` and mask-prediction
loss, then adds one conservative NeighborMix branch during training.

Default variant: `nm_scmae_mid`

- `alpha = 0.9`
- `neighbor_k = 5`
- `mix_neighbors = 4`
- `mix_weight = 0.5`
- `consistency_weight = 0.02`
- `target_mode = original`
- `mask_ratio = 0.4`

