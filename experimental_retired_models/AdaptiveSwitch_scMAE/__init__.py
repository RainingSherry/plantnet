"""AdaptiveSwitch-scMAE — dataset-adaptive sharp<->soft clustering.

Assembles every positive signal from the whole search into ONE model:
  - SHARP head (rank13 DEC): trainable centers + confidence-gated KL(sharpen(q)||q).
    Wins Quake (0.91) / Melanoma (0.65). Fails Macosko.
  - SOFT regime (rank29 dissection): fuzzy-core-KL + boundary-entropy(maximize on
    low-clusterability cells) + variance(VICReg anti-collapse), NO balance loss.
    Wins Macosko (0.637). Catastrophic on Quake (0.17).
  - SWITCH signal: kl_ref = mean KL(sharpen(q)||q) measured on the full set at each
    refresh. Low (~0.01, Quake) => data is DEC-clusterable => gate->1 (sharp).
    High (~0.7, Macosko) => sharpening won't converge => gate->0 (soft).
      gate = 1 / (1 + (kl_ref / kappa)^2)

Total cluster objective (cluster_scale ramps after warmup):
  L_cluster = gate * L_sharp + (1-gate) * L_soft
  L_soft always carries variance anti-collapse (cheap, safe both regimes).

This is the first attempt at a cross-dataset winner: Quake/Melanoma auto-lean
sharp, Macosko auto-leans soft, driven by a data-derived signal, not a fixed knob.
"""
