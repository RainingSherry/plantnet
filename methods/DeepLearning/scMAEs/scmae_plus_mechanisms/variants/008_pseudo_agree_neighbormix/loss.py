from __future__ import annotations

VARIANT_DEFAULTS = {
    "use_adaptive_mask": False,
    "use_prototype": False,
    "use_swav": False,
    "use_neighbormix": True,
    "neighbor_start_epoch": 20,
    "neighbor_mix_mode": "first",
    "neighbor_pseudo_filter": "same_cluster",
    "neighbor_pseudo_confidence_quantile": 0.0,
    "mix_alpha": 0.9,
    "mix_weight": 0.3,
}
