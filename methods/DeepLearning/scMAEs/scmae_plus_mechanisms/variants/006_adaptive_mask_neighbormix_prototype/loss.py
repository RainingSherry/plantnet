from __future__ import annotations

VARIANT_DEFAULTS = {
    "use_adaptive_mask": True,
    "use_prototype": True,
    "use_swav": False,
    "use_neighbormix": True,
    "mask_mode": "variance_adaptive",
    "prototype_start_epoch": 20,
    "neighbor_start_epoch": 20,
    "prototype_weight": 0.1,
    "mix_weight": 0.3,
    "neighbor_mix_mode": "first",
}

