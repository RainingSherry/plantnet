from __future__ import annotations

VARIANT_DEFAULTS = {
    "use_adaptive_mask": False,
    "use_prototype": False,
    "use_swav": False,
    "use_neighbormix": True,
    "neighbor_start_epoch": 20,
    "neighbor_mix_mode": "soft_first",
    "neighbor_soft_power": 1.0,
    "mix_alpha": 0.75,
    "mix_weight": 0.3,
}
