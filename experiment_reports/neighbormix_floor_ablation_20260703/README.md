# NeighborMix x std-floor ablation

This is an external experiment runner. It intentionally lives under
`experiment_reports/`, not under `methods/`, because it is an ablation scaffold
rather than formal benchmark code.

Question:

```text
Is NeighborMix compatible with the scMAE+DEC+std-floor diagnosis, or does local
neighbor smoothing conflict with fine-grained variance preservation?
```

Design for Macosko seeds 42, 43, and 44:

| Arm | `mix_mode` | `variance_weight` | Meaning |
|---|---|---:|---|
| DEC | `none` | `0.0` | sharp DEC control |
| DEC+floor | `none` | `0.02` | std-floor intervention |
| NeighborMix+DEC | `neighbor` | `0.0` | local denoising plus DEC |
| NeighborMix+DEC+floor | `neighbor` | `0.02` | compatibility test |

Implementation notes:

- The runner reuses `AdaptiveSwitch_scMAE` model/loss and `scMAE_family`
  loading/evaluation.
- NeighborMix is added only as an auxiliary denoising branch:
  original masked cell -> DEC/scMAE loss, mixed masked cell -> reconstruct the
  original log-expression target.
- Evaluation always extracts embeddings from original cells.

Summarize completed runs:

```bash
python experiment_reports/neighbormix_floor_ablation_20260703/summarize.py
```

Current Macosko result:

| Arm | ARI mean +/- sd | Effective dim PR mean +/- sd | Interpretation |
|---|---:|---:|---|
| DEC | 0.3298 +/- 0.0108 | 65.6 +/- 7.5 | collapsed low-variance control |
| DEC+floor | 0.7018 +/- 0.0051 | 105.7 +/- 12.4 | stable rescue by latent std-floor |
| NeighborMix+DEC | 0.3278 +/- 0.0057 | 40.9 +/- 5.6 | does not rescue collapse |
| NeighborMix+DEC+floor | 0.4322 +/- 0.0919 | 113.1 +/- 6.5 | geometrically compatible, but empirically hurts floor-only |

Conclusion: in this conservative compatibility test, NeighborMix should not be
used as the main explanatory mechanism for the Macosko rescue. The core
mechanism remains the latent std-floor intervention; NeighborMix is better
treated as a negative/boundary ablation showing that local smoothing can be
compatible with non-collapsed dimensions while still degrading fine-grained
cluster recovery.
