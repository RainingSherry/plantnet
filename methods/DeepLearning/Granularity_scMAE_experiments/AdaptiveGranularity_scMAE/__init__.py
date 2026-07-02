"""Adaptive-Granularity Cluster Target (AGCT) scMAE.

Central hypothesis (first-principles synthesis of the whole search):
  Quake needs SHARP partitioning (DEC wins, KL converges to ~0.01);
  Macosko needs SOFT / entropy-tolerant partitioning (fuzzy wins, DEC KL stalls
  at ~0.68 and distorts). A FIXED sharpness is why every prior method fails on
  one dataset.

Mechanism: make the cluster TARGET's sharpness a per-cell function of local
clusterability c_i:
    p_i^adaptive = c_i * sharpen(q_i) + (1 - c_i) * q_i
  - core/clusterable cell (c_i->1): target = sharp DEC distribution -> full DEC.
  - boundary/rare cell (c_i->0): target = its own q_i -> KL(q||q)=0 -> DEC
    self-suspends, the cell is left fuzzy (fuzzy-rough boundary region).

Gating the TARGET (not the loss weight) is the key difference from the failed
GatedNeighborMix line: it keeps DEC full-strength where data is clusterable
(so Quake stays 0.91) while removing DEC pressure exactly on the cells that
break it (so Macosko's fine subtypes are not force-merged).

Backbone reuses rank13's proven DEC scMAE + rank29's SVD-anchor fusion for
manifold stability on fine-grained data (Macosko).
"""
