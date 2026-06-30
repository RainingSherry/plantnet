from __future__ import annotations

VARIANT_DEFAULTS = {
    "use_adaptive_mask": False,
    "use_prototype": False,
    "use_swav": False,
    "use_neighbormix": True,
    "neighbor_start_epoch": 20,
    "neighbor_update_interval": 10,
    "neighbor_graph_embedding": "current",
    "neighbor_consensus_window": 3,
    "neighbor_consensus_min_hits": 2,
    "neighbor_mix_mode": "first",
    "mix_alpha": 0.9,
    "mix_weight": 0.3,
}
