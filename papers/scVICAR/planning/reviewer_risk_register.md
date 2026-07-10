# Reviewer-risk register

| Risk | Current control | Acceptance condition |
|---|---|---|
| Sixteen datasets were used during model development | Treat them as development evidence and freeze a new protocol | No confirmatory wording is attached to the historical table. |
| F and T previously used different backbones | One RG runner and one immutable common argument dictionary | Automated test proves common fields match across all variants. |
| Topology score is called reliability without calibration | Rename to topology-informed affinity | Manuscript contains no probability/calibration claim. |
| T is not consistently better than F | Make bad-edge robustness the T hypothesis | Report T vs F on clean and contaminated graphs, including failures. |
| Labels could leak through graph construction | The clean path fixes stress contamination to zero; labels are consulted only inside an explicit positive-contamination branch, with regression tests for clean numerical identity | Every formal clean config records `stress_bad_edge_ratio=0`; stress runs are separately namespaced and disclosed. |
| Known-K evaluation is optimistic | Declare ARI/KMeans-known-K primary and fixed Leiden secondary | No resolution or K oracle sweep. |
| Marker annotation could reuse evaluation labels | Fixed reference/evaluation split | Reference labels never enter evaluation-cluster marker calculation. |
| scCluBench scripts select favorable seeds | Reimplement protocol and aggregate every seed | No best-seed helper in paper code. |
| Local disk exhaustion | Dataset-at-a-time cache and atomic remote promotion | Scheduler stops below 5 GiB and never deletes before verification. |
| Name collision | CrossRef/OpenAlex preliminary check only | Repeat Scholar/GitHub/PyPI/trademark checks before submission. |
| External baseline cell filtering or row reordering | Every adapter emits `cell_ids.npy` and the canonicalizer compares exact `obs_names` order before evaluation | A mismatch raises before `COMPLETED`; scVI and PCA contract smokes pass. |
| External methods use known K during training | Separate label-free representation baselines from known-K clustering-oriented baselines | scDCC/scDeepCluster/scRCL are explicitly disclosed; per-cell labels are withheld from optimization/checkpoint selection. |
