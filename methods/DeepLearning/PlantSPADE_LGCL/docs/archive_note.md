# Archive / Exploratory Code

The main paper protocol is under `methods/DeepLearning/PlantSPADE_LGCL`.

Removed during cleanup:

- old non-LGCL `methods/DeepLearning/PlantSPADE`
- DOLORIS, GraphDiffusion, DiffusionBridge, and maskdiffusion archive code
- root-level GraphDiffusion/DOLORIS helper scripts
- old Diffusion/DOLORIS benchmark result folders

The following exploratory code remains excluded from main runners, main tables, and the README main flow:

- unstable or dependency-heavy GNN entries such as `scCDCG`, `scDSC`, `scGNN`
- foundation model experiments such as `scGPT`, `GeneFormer`, `GeneCompass`, and `scPlantLLM`

New paper-facing experiments should be added through the PlantSPADE-LGCL config and runner stack.
