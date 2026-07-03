# NeighborMix x std-floor ablation

This is an external experiment runner. It intentionally lives under
`experiment_reports/`, not under `methods/`, because it is an ablation scaffold
rather than formal benchmark code.

Question:

```text
Is NeighborMix compatible with the scMAE+DEC+std-floor diagnosis, or does local
neighbor smoothing conflict with fine-grained variance preservation?
```

Design for Macosko seed 42:

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
