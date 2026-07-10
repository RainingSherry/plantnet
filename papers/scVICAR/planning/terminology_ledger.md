# Terminology ledger

| Canonical term | First-use definition | Disallowed/drifting forms | Decision |
|---|---|---|---|
| scRNA-seq | single-cell RNA sequencing (scRNA-seq) | scRNAseq, single cell RNA-seq | Define once, then use scRNA-seq. |
| scVICAR | single-cell Vicinal Corruption with Anchor Recovery (scVICAR) | NeighborMix framework | Umbrella method name. |
| scVICAR-F | fixed graph-vicinal anchor recovery | NeighborMix-scMAE in main prose | Legacy name appears only when connecting to prior experiments. |
| scVICAR-T | topology-adaptive graph-vicinal anchor recovery | RG as reliability probability | Use topology-adaptive; the score is not calibrated. |
| topology-informed affinity | Analytic edge affinity from cosine similarity, mutual KNN and SNN | reliability probability, confidence probability | Canonical edge term. |
| cell-wise perturbation budget | Node-specific gate controlling vicinal displacement | learned gate, uncertainty gate | The current gate is analytic and static. |
| anchor recovery | Reconstruction of the original cell from a corrupted vicinal view | mixed-label reconstruction | Core training target. |
| development benchmark | Historical parameter-development datasets | test benchmark, held-out benchmark | Applies to the original 16-dataset table. |
| confirmatory evaluation | Frozen-protocol evaluation on six previously unused datasets | validation tuning set | No tuning is allowed. |
| ARI | adjusted Rand index (ARI) | adjusted rand score | Primary endpoint. |
| low-label transductive annotation | Frozen-embedding linear probe within a dataset | transfer learning, cross-dataset transfer | Exact downstream scope. |

