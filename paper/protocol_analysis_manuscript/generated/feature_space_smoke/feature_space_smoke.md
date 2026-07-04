# Feature-space Smoke Summary

Status: generated development-only feature-space smoke. This is not validation evidence.

| feature space | genes | student params | known-K ARI | known-K NMI | fixed-Leiden ARI | fixed-Leiden F1 |
|---|---:|---:|---:|---:|---:|---:|
| HVG 2000 | 2000 | 9080688 | 0.485571 | 0.724467 | 0.463516 | 0.605773 |
| Full-gene stress | 23341 | 1101675865 | 0.443727 | 0.672599 | 0.389245 | 0.454948 |

## Paired Deltas

- Full-gene minus HVG known-K ARI: -0.041844.
- Full-gene minus HVG fixed-Leiden ARI: -0.074271.
- Full-gene minus HVG fixed-Leiden F1: -0.150825.
- Full-gene/HVG student parameter ratio: 121.32x.
- Full-gene/HVG gene-count ratio: 11.67x.

## Claim Boundary

This supports keeping HVG 2000 as the current dense-MLP protocol default under development evidence. It does not validate HVG 2000, reject all full-gene approaches, or evaluate sparse/gene-token full-gene architectures.
