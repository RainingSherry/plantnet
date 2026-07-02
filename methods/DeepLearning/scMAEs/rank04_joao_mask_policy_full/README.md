# rank04_joao_mask_policy_full

Independent-full scMAE candidate inspired by JOAO.

JOAO's transferable idea is dynamic augmentation selection. Here the search
space is deliberately small and scMAE-specific: swap mask, absorbing mask,
module mask, and dropout-like mask. The controller updates policy probabilities
from per-policy training loss while retaining a small uniform prior.

NeighborMix is not used; `mixed_cell_fraction` is always `0.0`.

