from __future__ import annotations

VARIANT_DEFAULTS = {
    "use_adaptive_mask": False,
    "use_prototype": False,
    "use_swav": False,
    "use_neighbormix": True,
    "neighbor_start_epoch": 30,
    "neighbor_update_interval": 10,
    "neighbor_boundary_protect": True,
    "neighbor_boundary_confidence_quantile": 0.20,
    "neighbor_boundary_rare_quantile": 0.25,
    "neighbor_boundary_score_threshold": 0.84,
    "neighbor_mix_mode": "first",
    "mix_alpha": 0.9,
    "mix_weight": 0.3,
}
