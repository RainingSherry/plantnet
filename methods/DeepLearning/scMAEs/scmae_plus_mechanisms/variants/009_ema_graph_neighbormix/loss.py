from __future__ import annotations

VARIANT_DEFAULTS = {
    "use_adaptive_mask": False,
    "use_prototype": False,
    "use_swav": False,
    "use_neighbormix": True,
    "neighbor_start_epoch": 20,
    "neighbor_update_interval": 10,
    "neighbor_graph_embedding": "ema",
    "neighbor_embedding_ema_decay": 0.8,
    "neighbor_mix_mode": "first",
    "mix_alpha": 0.9,
    "mix_weight": 0.3,
}
