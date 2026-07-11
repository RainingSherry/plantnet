# scVICAR manuscript argument audit

## One-sentence argument

In masked single-cell representation learning, scVICAR uses bounded
graph-vicinal anchor recovery to exploit local cell geometry, achieves
competitive performance across development and frozen benchmarks, and exposes
when topology adaptation detects rather than necessarily overcomes graph risk.

## Evidence ladder

| Reader question | Evidence used in the main text | Defensible conclusion |
|---|---|---|
| Is the method competitive? | 15-dataset development benchmark; six-dataset external table | F and T are broadly competitive; T has the highest mean external ARI in the evaluated comparison. |
| Is the gain caused by the proposed corruption? | Six matched variants with identical backbone and budget | The incremental mean effect is small and dataset-dependent; capacity does not explain component differences. |
| Is local geometry necessary? | RandomMix negative control | Arbitrary convex smoothing does not reproduce the local method's behavior. |
| Which adaptive component matters? | Edge-only and gate-only ablations | Cell-wise perturbation control is the stronger descriptive clean-graph component. |
| Does topology adaptation measure graph risk? | Bad-edge injection, affinity AUROC, gate--purity association | Affinity and gating respond to injected contamination; ARI protection remains conditional. |
| Does the embedding retain downstream value? | Marker recovery, marker annotation, 10%/30% frozen probes | Low-label linear separability shows small positive trends; marker coherence is not uniformly improved. |
| Are conclusions sensitive to labels or clustering? | Full-label cohort and fixed-resolution Leiden | The direction varies with evaluation regime; full-label T is promising but not conclusive. |

## Claims retained

1. Graph-vicinal anchor recovery is a bounded local corruption framework.
2. scVICAR is competitive within the evaluated development and frozen
   benchmarks.
3. Topology-informed affinity and cell-wise gating detect graph contamination
   risk.
4. Low-label probes provide a small, directionally consistent downstream
   signal.

## Claims deliberately excluded

1. scVICAR-T is universally superior to scVICAR-F or NoMix.
2. The topology affinity is a calibrated probability of edge correctness.
3. The theoretical bounds prove improved clustering.
4. Low-label probes demonstrate cross-dataset transfer.
5. Clustering improvements guarantee better marker recovery.

## Main reviewer risks after revision

| Risk | Current response | Remaining action before submission |
|---|---|---|
| Matched effects are small | Lead with external competitiveness, then use matched analysis for attribution rather than overall ranking. | Preserve every dataset-level result and avoid significance language. |
| Development archive informed method choice | Separate it from frozen inference and disclose full history in the supplement. | Ensure cover letter and response use the same evidence boundary. |
| T is not uniformly robust | Frame T as measurable conditional risk control. | Do not promote 2/3 stress wins into a universal claim. |
| Gate-only has the best descriptive ARI | Interpret it as a design finding, not a post-hoc primary endpoint. | Consider a future gate-focused model only in later work. |
| Marker endpoints are weak | Separate predictive separability from marker coherence. | Retain negative marker results in the main paper. |
| Engineering detail dominates | Keep scientific provenance in the main paper and move storage mechanics to the supplement. | Maintain a concise availability section. |

## Q1 positioning

The paper's value is not a single large matched ARI delta. Its stronger
contribution is the combination of a bounded formulation, a competitive full
pipeline, component-level causal controls, direct graph-contamination tests,
and downstream evidence with explicit negative results. This package supports
a rigorous algorithmic paper when the prose keeps competitiveness, attribution,
mechanism, and biological utility as distinct claims.
