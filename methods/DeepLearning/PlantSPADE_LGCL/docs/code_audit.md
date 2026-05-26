# Code Audit: Runnable Entrypoints

This audit separates paper-facing entries from archive/exploratory code.

## Main Protocol Entrypoints

- `methods/DeepLearning/PlantSPADE_LGCL/run_plantspade.py`
- `methods/DeepLearning/PlantSPADE_LGCL/scripts/profile_datasets.py`
- `methods/DeepLearning/PlantSPADE_LGCL/scripts/run_single.py`
- `methods/DeepLearning/PlantSPADE_LGCL/scripts/run_suite.py`
- `methods/DeepLearning/PlantSPADE_LGCL/scripts/aggregate_results.py`

## Retained Main Baseline Entrypoints

Traditional:

- `methods/Traditional/ScanpyStandard/run.py`
- `methods/Traditional/Leiden/run.py`
- `methods/Traditional/Louvain/run.py`
- `methods/Traditional/sc3/run.py`

The new main protocol uses `traditional_pca` from `run_single.py` to avoid label-tuned Leiden/Louvain leakage. The older traditional scripts are retained for compatibility and supplementary checks.

Deep baselines:

- `methods/DeepLearning/PhytoCluster/run.py`
- `methods/DeepLearning/scVI/run.py`
- `methods/DeepLearning/scMAE/run.py`

The new runner re-evaluates their embeddings with the fixed protocol.

## Archive / Exploratory Entrypoints

These are not part of the main table:

- `methods/DeepLearning/PlantSPADE/run_plantspade.py`
- `methods/DeepLearning/_archive/**`
- `methods/GNN/scCDCG/run.py`
- `methods/GNN/scDSC/run.py`
- `methods/GNN/scGNN/run.py`
- `methods/GNN/AttentionAE-sc/run.py`
- `methods/Foundation/scGPT/run.py`
- `methods/Foundation/GeneFormer/run.py`
- `methods/Foundation/GeneCompass/run.py`
- root-level `infer_graphdiffusion.py`
- root-level `infer_graphdiffusion_ckpt.py`
- root-level `launch_doloris_eval.sh`
- root-level `run_infer.py`

Rationale: these entries are diffusion, graph-diffusion, unstable dependency, or foundation-model explorations. They can remain in the repository for history, but they should not be invoked by PlantSPADE-LGCL paper runs.
