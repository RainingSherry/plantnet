# Mask-ratio Smoke Summary

Status: generated development-only mask-ratio sensitivity across three seeds. This is not validation evidence.

| mask ratio | runs | masked effective mean | global effective mean | known-K ARI mean | known-K ARI std | fixed-Leiden ARI mean | fixed-Leiden ARI std |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.2 | 3 | 0.191674 | 0.038335 | 0.462024 | 0.048981 | 0.408229 | 0.048769 |
| 0.4 | 3 | 0.191684 | 0.076674 | 0.490818 | 0.006156 | 0.430972 | 0.061284 |
| 0.6 | 3 | 0.191593 | 0.114956 | 0.470444 | 0.016125 | 0.480897 | 0.027309 |

## Sensitivity Notes

- Best mean known-K ARI in this smoke: mask ratio 0.4 (0.490818).
- Best mean fixed-Leiden ARI in this smoke: mask ratio 0.6 (0.480897).
- known-K ARI mean range: 0.028794.
- fixed-Leiden ARI mean range: 0.072668.
- maximum known-K seed standard deviation: 0.048981.
- maximum fixed-Leiden seed standard deviation: 0.061284.
- Masked-position effective corruption mean range: 0.000091.
- Global effective change estimate mean range: 0.076621.

## Claim Boundary

This supports reporting protocol sensitivity around nominal mask ratio. It must not be used to tune the frozen validation mask ratio or to claim that a different ratio is generally superior.
