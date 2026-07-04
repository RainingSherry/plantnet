# Feature-space stress smoke: HVG 2000 vs full-gene input

Date: 2026-06-26

Status: development-only smoke. This is not validation evidence and must not be used as a publication-level claim.

## Purpose

Test whether the corrected CAAM/scMAE protocol should continue to treat HVG 2000 as the main feature-space default, or whether full-gene input looks promising enough to justify a broader experiment.

## Frozen smoke setup

- Dataset: `Quake_Smart-seq2_Lung`
- Data path: `/data/luolie/biopipeline/dimension-reduction/plantnet/data/processed/Quake_Smart-seq2_Lung.h5ad`
- Variant: `control`
- Corruption: `scmae_shuffle`
- Mask selector: random
- Seed: `42`
- Epochs: `3`
- Input mode: `log1p`
- Scale input: `false`
- Benchmark mode: `true`
- Output root: `/tmp/caam_feature_space_smoke/dev_20260626_gpu`

Only `n_top_genes` was changed:

- HVG run: `n_top_genes=2000`
- Full-gene stress run: `n_top_genes=0`

## Results

| feature space | genes | student params | known-K ARI | known-K NMI | known-K F1 | fixed-Leiden ARI | fixed-Leiden NMI | fixed-Leiden F1 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| HVG 2000 | 2000 | 9,080,688 | 0.485571 | 0.724467 | 0.471594 | 0.463516 | 0.740922 | 0.605773 |
| full-gene stress | 23341 | 1,101,675,865 | 0.443727 | 0.672599 | 0.454177 | 0.389245 | 0.691040 | 0.454948 |

Full-gene minus HVG deltas:

- known-K ARI: `-0.041844`
- fixed-Leiden ARI: `-0.074271`
- student parameter ratio: `121.32x`

## Execution notes

An initial CPU full-gene attempt completed preprocessing and donor-candidate construction but was interrupted after exceeding the intended small-smoke cost boundary. The completed paired comparison was then rerun outside the sandbox on GPU 1, still writing only to `/tmp` and still using a single development dataset.

## Interpretation

This smoke supports keeping HVG 2000 as the current main protocol default. Full-gene input greatly increases parameter count under the current MLP encoder and does not improve either known-K or fixed-Leiden clustering in this single development smoke.

This result should not be generalized to all full-gene or gene-token architectures. It only says that full-gene dense MLP input is not an attractive default for the current CAAM/scMAE control path.

## Claim boundary

Safe:

```text
A single development smoke suggests that dense full-gene MLP input is much more expensive and weaker than HVG 2000 under the current corrected control protocol.
```

Unsafe:

```text
HVG 2000 is validated.
Full-gene input is always worse.
Gene-token or sparse full-gene models will fail.
This result supports a publication-level feature-space claim.
```

## Next decision

Do not expand full-gene dense MLP to a larger matrix by default. If feature-space remains central to the manuscript, the next safe step is either:

1. add this as a clearly scoped supplementary development diagnostic; or
2. design a separate sparse/gene-token feature-space experiment, treated as a new mechanism gate rather than a continuation of the dense MLP protocol.
