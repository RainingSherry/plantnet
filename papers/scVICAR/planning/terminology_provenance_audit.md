# Terminology provenance and style audit

This audit distinguishes established field terminology from names introduced
for the scVICAR manuscript.

| Term | Status | Manuscript decision |
|---|---|---|
| masked autoencoder | Established | Retain. |
| vicinal risk minimization | Established | Use only when discussing the learning principle. |
| mixup | Established | Retain with its original citation. |
| KNN, mutual KNN, shared-nearest-neighbor overlap | Established | Retain. |
| graph diffusion, random-walk Laplacian | Established | Retain in the theoretical interpretation. |
| Lipschitz perturbation bound | Established mathematical language | Retain with explicit assumptions. |
| scVICAR | Author-defined model name | Retain. |
| graph-vicinal anchor recovery | Author-defined method description | Retain and define once. It concisely states graph-local corruption and recovery of the original anchor. |
| topology-informed affinity | Author-defined descriptive name | Retain as an uncalibrated analytic edge weight; avoid calling it reliability or probability. |
| cell-wise perturbation budget | Author-defined metaphor | Replace with **cell-specific mixing coefficient**. The latter names the implemented scalar directly. |
| conditional risk control | Interpretive phrase | Remove. State the observed response to inaccurate edges directly. |
| diagnostically testable regularizer | Interpretive phrase | Remove. Report AUROC, correlation, and clustering response directly. |
| bounded trust region | Imported optimization metaphor | Remove from the main narrative. The convex-hull and displacement bounds carry the precise meaning. |
| falsifiable empirical prediction | Meta-scientific phrasing | Replace with the concrete experiment or measured quantity. |

## Style rule adopted

The revised manuscript advances the argument with affirmative statements.
Negative constructions remain only where they encode a necessary protocol or
scope constraint. The main text avoids the template `not X, but Y`, repeated
self-defense, and abstract labels for the evidence structure.
