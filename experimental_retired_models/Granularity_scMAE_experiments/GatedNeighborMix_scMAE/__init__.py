"""ReliabilityGatedNeighborMix-scMAE candidate line.

Reliability-gated fusion of three proven-effective mechanisms:
  - scMAE backbone (swap-noise + mask discriminator + weighted reconstruction)
  - NeighborMix expression-augment branch (per-cell gated)
  - DEC-style trainable cluster centers (per-cell gated KL)

A single per-cell reliability field r_i in [0,1] modulates BOTH the NeighborMix
alpha and the DEC KL weight, so every smoothing mechanism retreats together on
rare/boundary cells and hands them back to the self-anchored scMAE backbone.
"""
